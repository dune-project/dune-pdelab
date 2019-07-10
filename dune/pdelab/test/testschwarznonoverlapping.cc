// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// C++ includes
#include<math.h>
#include<iostream>

// dune-common includes
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/timer.hh>
#include<dune/common/enumset.hh>
#include<dune/common/parallel/indexset.hh>
#include<dune/common/parallel/plocalindex.hh>
#include<dune/common/parallel/interface.hh>
#include<dune/common/parallel/remoteindices.hh>
#include<dune/common/parallel/communicator.hh>
#include<dune/common/parallel/variablesizecommunicator.hh>
// dune-geometry includes
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
// dune-grid includes
#include<dune/grid/onedgrid.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#include<dune/alugrid/dgf.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
// dune-istl includes
#include<dune/istl/bvector.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/io.hh>

#include<dune/pdelab/backend/istl/geneo/schwarznonoverlapping.hh>

/**
 * \brief communication data handle for adding data
 *
 * V vector container. Elements of this vector are sent around. Shoul work for std::vector and ISTL::BlockVector
 */
template<typename GV, typename Vector>
class VectorAddDataHandle
  : public Dune::CommDataHandleIF<VectorAddDataHandle<GV,Vector>,typename Vector::block_type>
{
  const typename GV::IndexSet& indexset;
  Vector& v;

public:
  typedef typename Vector::value_type DataType;

  VectorAddDataHandle (const GV& gv, Vector& v_)
    : indexset(gv.indexSet()), v(v_)
  {}

  bool contains (int dim, int codim) const
  {
    return (codim==dim);
  }

  bool fixedSize (int dim, int codim) const
  {
    return true;
  }

  template<class EntityType>
  size_t size (const EntityType& e) const
  {
    return 1;
  }

  template<class MessageBufferImp, class EntityType>
  void gather (MessageBufferImp& buff, const EntityType& e) const
  {
    buff.write(v[indexset.index(e)]);
  }

  template<class MessageBufferImp, class EntityType>
  void scatter (MessageBufferImp& buff, const EntityType& e, size_t n)
  {
    DataType x;
    buff.read(x);
    v[indexset.index(e)] += x;
  }
};

//! Operator for the non-overlapping parallel case
/**
 * Calculate \f$y:=Ax\f$.
 *
 * \tparam GV  Grid View
 * \tparam Matrix   Type of the matrix.  Should be one of the ISTL matrix types.
 * \tparam Vector   Type of the vectors the matrix is applied to.
 */
template<typename GV, typename Matrix, typename Vector>
class NonoverlappingOperator
  : public Dune::LinearOperator<Vector,Vector>
{
  const GV& gv;

public:
  //! export type of matrix
  using matrix_type = Matrix;
  //! export type of vectors the matrix is applied to
  using domain_type = Vector;
  //! export type of result vectors
  using range_type = Vector;
  //! export type of the entries for x
  typedef typename Vector::field_type field_type;

  //! need to have this method now
  virtual Dune::SolverCategory::Category category() const
  {
    return Dune::SolverCategory::nonoverlapping;
  }

  //! Construct a non-overlapping operator
  /**
   * \param gfs_ GridFunctionsSpace for the vectors.
   * \param A    Matrix for this operator.  This should be the locally
   *             assembled matrix.
   *
   * \note The constructed object stores references to all the objects
   *       given as parameters here.  They should be valid for as long as
   *       the constructed object is used.  They are not needed to
   *       destruct the constructed object.
   */
  NonoverlappingOperator (const GV& gv_, const Matrix& A_)
    : gv(gv_), A(A_)
  { }

  //! apply operator
  /**
   * Compute \f$y:=A(x)\f$ on this process.
   * It is assumed that x is consistent, A is additive
   * then y is additive after multiplication and is immediately made consistent
   */
  virtual void apply (const Vector& x, Vector& y) const
  {
    A.mv(x,y); // A is additive and x is consistent; produces additive result
  }

  //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
  /**
   * Compute \f$y += \alpha A(x)\f$ on this process
   * It is assumed that x is consistent, A is additive and y is consistent
   */
  virtual void applyscaleadd (field_type alpha, const Vector& x, Vector& y) const
  {
    A.usmv(alpha,x,y); // A is additive and x is consistent; produces additive result
  }

private:
  const Matrix& A;
};


// parallel scalar product assuming no overlap
template<class GV, class Vector>
class NonoverlappingScalarProduct : public Dune::ScalarProduct<Vector>
{
public:
  //! export types
  typedef Vector domain_type;
  typedef typename Vector::field_type field_type;

  //! define the category
  virtual Dune::SolverCategory::Category category() const
  {
    return Dune::SolverCategory::nonoverlapping;
  }

  /*! \brief Constructor needs to know the grid view
   */
  NonoverlappingScalarProduct (const GV& gv_, const Vector& x)
    : gv(gv_)
  {
  }

  /*! \brief Dot product of two vectors.
   * It is assumed that one vector is consistent and the other is additive
   */
  virtual field_type dot (const Vector& x, const Vector& y)
  {
    field_type sum = 0.0;
    for (typename Vector::size_type i=0; i<x.N(); i++) sum += x[i]*y[i];
    auto sumsum = gv.comm().sum(sum);
    return sumsum;
  }

  /*! \brief Norm of a right-hand side vector.
   * It is assumed that x is additive.
   * This operation requires a local communication in addition
   * to make one argument consistent.
   */
  virtual double norm (const Vector& x)
  {
    // if (!check_vector_isfinite(x))
    //   std::cout << gv.comm().rank() << ": NaN in x detected in NonoverlappingScalarProduct.norm" << std::endl;
    Vector y(x);
    VectorAddDataHandle<GV,Vector> adddh(gv,y);
    if (gv.comm().size()>1)
      gv.communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);
    auto sp = static_cast<double>(this->dot(x,y));
    auto rv = std::sqrt(std::abs(sp)); // due to roundoff this may become negative for close to zero norms
    return rv;
  }

private:
  const GV gv;
};


int main(int argc, char **argv)
{
  try{
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper&
      helper = Dune::MPIHelper::instance(argc, argv);

#if HAVE_UG
    // load a grid file in parallel, loadbalance and refine
    using Grid = Dune::UGGrid<2>;
    std::shared_ptr<Grid> pgrid;
    const int dimension = Grid::dimension;

    // create and partition structured grid on a square process grid
    int m = std::sqrt(helper.size());
    if (m*m != helper.size())
      {
        if (helper.rank()==0) std::cout << helper.size() << " is not a square number of processors" << std::endl;
        return 1;
      }

    if (helper.rank()==0)
      {
        Dune::StructuredGridFactory<Grid> factory;
        Dune::FieldVector<double,dimension> lowerLeft; lowerLeft = -1.0;
        Dune::FieldVector<double,dimension> upperRight; upperRight = 1.0;
        auto cells = 5;
        std::array<unsigned int,dimension> cellsarray;
        for (int i=0; i<dimension; i++) cellsarray[i] = cells*m;
        pgrid = factory.createSimplexGrid(lowerLeft,upperRight,cellsarray);
        double HX = (upperRight[0]-lowerLeft[0])/m;
        double HY = (upperRight[1]-lowerLeft[1])/m;
        std::vector<Grid::Rank> target(pgrid->leafGridView().indexSet().size(0));
        for (const auto& e : elements(pgrid->leafGridView(),Dune::Partitions::all))
          {
            auto center = e.geometry().center();
            double X = (center[0]-lowerLeft[0])/HX;
            double Y = (center[1]-lowerLeft[1])/HY;
            int px = std::floor(X);
            int py = std::floor(Y);
            int p = py*m+px;
            //std::cout << "center=" << center << " px=" << px << " py=" << py << " p=" << p << std::endl;
            target[pgrid->leafGridView().indexSet().index(e)] = p;
          }
        pgrid->loadBalance(target,0);
      }
    else {
      Dune::GridFactory<Grid> factory;
      pgrid = std::shared_ptr<Grid>(factory.createGrid());
      std::vector<Grid::Rank> target(0);
      pgrid->loadBalance(target,0);
    }
    pgrid->globalRefine(2);

    // extract the leaf grid view
    auto gv = pgrid->leafGridView();
    using GV = decltype(gv); // and get the type

    // define the PDE problem
    auto f = [&](const auto& x){ // define the right hand side of the PDE
      return 0.0;};
    auto predicate = [&](const auto& x){ // return true if position x is on the Dirichlet boundary
      return (x[0]<-1.0+1e-6) || (x[0]>1.0-1e-6) || (x[1]<-1.0+1e-6) || (x[1]>1.0-1e-6);
    };
    auto g = [&](const auto& x){ // value at the Dirichlet boundary
      return x[0]*x[1];
    };

    // set up structure of the sparse matrix
    using K = double;
    const int components = 1;
    using ScalarVector = Dune::BlockVector<Dune::FieldVector<K,1>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<K,components>>;
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<K,components,components>>;
    auto& indexset = gv.indexSet();
    int avg = (dimension==1) ? 3 : ((dimension==2) ? 7 : 14);
    Matrix A(indexset.size(dimension),indexset.size(dimension),avg,0.05,Matrix::implicit);
    for (const auto& e : elements(gv,Dune::Partitions::all))
      for (int i=0; i<=dimension; i++)
        for (int j=0; j<=dimension; j++)
          A.entry(indexset.subIndex(e,i,dimension),indexset.subIndex(e,j,dimension)) = 0.0;
    auto stats = A.compress();
    //std::cout << "Proc " << helper.rank() << ": max_entries_per_row=" << stats.maximum << " mem_ratio=" << stats.mem_ratio << std::endl;

    // fill entries of local stiffness matrix
    // prepare basis functions on the reference element at reference element center
    enum {n=dimension+1};
    auto firstelement = *gv.template begin<0>();
    auto firstgeo = firstelement.geometry();
    auto rule1 = Dune::QuadratureRules<Grid::ctype,dimension>::rule(firstgeo.type(),1); // get lowest order quadrature rule
    if (rule1.size()>1) {
      std::cout << "Wrong quadrature rule!" << std::endl;
      exit(1);
    }
    auto weight = rule1[0].weight(); // quadrature weight on refelem
    auto qp = rule1[0].position();   // center of mass in reference element
    double gradphihat[dimension][n] = {{0.0}}; // gradients of basis functions on refelem
    for (int i=0; i<dimension; i++) {
      gradphihat[i][0] = -1.0;
      gradphihat[i][i+1] = 1.0;
    }
    for (const auto& e : elements(gv,Dune::Partitions::interior)) // nonoverlapping method: assemble only interior
      {
        // we use lowest order quadrature
        auto geo = e.geometry();
        auto S = geo.jacobianInverseTransposed(qp);
        auto factor = weight*geo.integrationElement(qp);

        // compute gradients of basis functions in transformed element
        double gradphi[dimension][n] = {{0.0}}; // coordinate x #basisfct
        for (int i=0; i<dimension; i++) // rows of S
          for (int k=0; k<dimension; k++) // columns of S
            for (int j=0; j<n; j++) // columns of gradhat
              gradphi[i][j] += S[i][k] * gradphihat[k][j];

        // compute gradphi^T * gradphi
        double A_e[n][n] = {{0.0}};
        for (int i=0; i<n; i++)
          for (int k=0; k<dimension; k++)
            for (int j=0; j<n; j++)
              A_e[i][j] += gradphi[k][i]*gradphi[k][j];

        // store in result
        for (int i=0; i<n; i++)
          for (int j=0; j<n; j++)
            for (int k=0; k<components; k++)
              A[indexset.subIndex(e,i,dimension)][indexset.subIndex(e,j,dimension)][k][k] += A_e[i][j]*factor;
      }

    // set up right hand side vector
    Vector b(indexset.size(dimension));
    b = 0.0;
    auto phihat = [&](int i, const auto& x){ // define basis functions on reference element
      if (i>0) return x[i-1];
      double s=0.0; for (int i=0; i<dimension; i++) s+=x[i];
      return 1.0-s;};
    for (const auto& c : elements(gv,Dune::Partitions::interior)) // nonoverlapping method: assemble only interior
      {
        auto geo = c.geometry();
        auto quadrature = Dune::QuadratureRules<Grid::ctype,dimension>::rule(geo.type(),4);
        for (const auto& q : quadrature) {
          auto factor = f(geo.global(q.position()))*geo.integrationElement(q.position())*q.weight();
          for (int i=0; i<n; i++)
            for (int k=0; k<components; k++)
              b[indexset.subIndex(c,i,dimension)][k] += phihat(i,q.position())*factor;
        }
      }

    // prepare for Dirichlet boundary conditions
    Vector x(indexset.size(dimension));
    x=0.0;
    std::vector<bool> dirichlet(indexset.size(dimension)); // store a flag for each vertex;
    bool floating = true;
    for (const auto& v : vertices(gv,Dune::Partitions::all))
      if (v.partitionType()==Dune::GhostEntity)
        {
          dirichlet[indexset.index(v)] = true;
          for (int k=0; k<components; k++)
            b[indexset.index(v)][k] = x[indexset.index(v)][k] = 0.0;
        }
      else if (predicate(v.geometry().corner(0)))
        {
          floating = false;
          dirichlet[indexset.index(v)] = true; // set the flag
          for (int k=0; k<components; k++)
            b[indexset.index(v)][k] = x[indexset.index(v)][k] = g(v.geometry().corner(0));
        }
    for (size_t i=0; i<A.N(); i++)
      if (dirichlet[i])
        {
          auto cIt = A[i].begin();
          auto cEndIt = A[i].end();
          for (; cIt!=cEndIt; ++cIt)
            if (i==cIt.index())
              {
                for (int k=0; k<components; k++)
                  for (int l=0; l<components; l++)
                    (*cIt)[k][l] = (k==l) ? 1.0 : 0.0;
              }
            else
              {
                for (int k=0; k<components; k++)
                  for (int l=0; l<components; l++)
                    (*cIt)[k][l] = 0.0;
              }
        }

    // provide a partitionType for each degree of freedom. Should be interior, border or ghost
    std::vector<Dune::PDELab::EPISAttribute> partitiontype(indexset.size(dimension));
    for (const auto& v : vertices(gv,Dune::Partitions::all))
      {
        if (v.partitionType()==Dune::InteriorEntity) partitiontype[indexset.index(v)] = Dune::PDELab::EPISAttribute::interior;
        if (v.partitionType()==Dune::BorderEntity) partitiontype[indexset.index(v)] = Dune::PDELab::EPISAttribute::border;
        if (v.partitionType()==Dune::OverlapEntity) partitiontype[indexset.index(v)] = Dune::PDELab::EPISAttribute::overlap;
        if (v.partitionType()==Dune::GhostEntity) partitiontype[indexset.index(v)] = Dune::PDELab::EPISAttribute::ghost;
      }

    // provide global id for each dof
    auto& globalidset = gv.grid().globalIdSet();
    using GlobalId = Grid::GlobalIdSet::IdType;
    std::vector<GlobalId> globalid(indexset.size(dimension));
    for (const auto& v : vertices(gv,Dune::Partitions::all))
      globalid[indexset.index(v)] = globalidset.id(v);

    // find ranks with which we might share data
    int overlap = 2;
    auto allmyneighborsvec = Dune::PDELab::findNeighboringRanks(gv,overlap);

    // solve the linear system iteratively
    // operators and scalar product
    NonoverlappingOperator<GV, Matrix,Vector> linearOperator(gv,A);
    NonoverlappingScalarProduct<GV,Vector> scalarproduct(gv,x);

    // preconditioner
    auto coarsespace = true;
    using CC = typename GV::CollectiveCommunication;
    Dune::PDELab::NonoverlappingSchwarzPreconditioner<CC,GlobalId,Matrix,Vector>
      preconditioner(gv.comm(),allmyneighborsvec,A,floating,partitiontype,globalid,avg,overlap,coarsespace);

    // solver
    int verbose=0;
    if (gv.comm().rank()==0) verbose=2;
    Dune::CGSolver<Vector> solver(linearOperator,scalarproduct,preconditioner,1e-6,2500,verbose);
    Dune::InverseOperatorResult stat;

    solver.apply(x,b,stat);
#endif

    return 0;
  }

  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }

}
