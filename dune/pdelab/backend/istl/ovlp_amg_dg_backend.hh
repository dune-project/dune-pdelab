#ifndef DUNE_PDELAB_OVLP_AMG_DG_BACKEND_HH
#define DUNE_PDELAB_OVLP_AMG_DG_BACKEND_HH

#include <dune/istl/matrixmatrix.hh>

#include <dune/grid/common/datahandleif.hh>

#include <dune/pdelab/backend/istl/vector.hh>
#include <dune/pdelab/backend/istl/bcrsmatrix.hh>
#include <dune/pdelab/backend/istl/ovlpistlsolverbackend.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>

namespace Dune {
  namespace PDELab {
    //***********************************************************
    // A data handle / function to communicate matrix entries
    // in the overlap. It turned out that it is actually not
    // required to do that, but anyway it is there now.
    //***********************************************************

    /** Data handle to build up local index to global id and global id to local index map for codim 0 in overlap
     */
    template<class GFS>
    class LocalGlobalMapDataHandle
      : public Dune::CommDataHandleIF<LocalGlobalMapDataHandle<GFS>,int>
    {
    public:
      // some types from the grid view
      typedef typename GFS::Traits::GridView GV;
      typedef typename GV::IndexSet IndexSet;
      typedef typename IndexSet::IndexType IndexType;
      typedef typename GV::Grid Grid;
      typedef typename Grid::Traits::GlobalIdSet GlobalIdSet;
      typedef typename GlobalIdSet::IdType IdType;

      //! export type of data for message buffer
      typedef int DataType;

      // map types
      typedef std::map<IndexType,IdType> LocalToGlobalMap;
      typedef std::map<IdType,IndexType> GlobalToLocalMap;

      //! returns true if data for this codim should be communicated
      bool contains (int dim, int codim) const
      {
        return (codim==0);
      }

      //! returns true if size per entity of given dim and codim is a constant
      bool fixedsize (int dim, int codim) const
      {
        return true;
      }

      /*! how many objects of type DataType have to be sent for a given entity

        Note: Only the sender side needs to know this size.
      */
      template<class EntityType>
      size_t size (EntityType& e) const
      {
        return 1;
      }

      //! pack data from user to message buffer
      template<class MessageBuffer, class EntityType>
      void gather (MessageBuffer& buff, const EntityType& e) const
      {
        // fill map
        IndexType myindex = indexSet.index(e);
        IdType myid = globalIdSet.id(e);
        l2g[myindex] = myid;
        g2l[myid] = myindex;
        //std::cout << gv.comm().rank() << ": gather local=" << myindex << " global=" << myid << std::endl;

        buff.write(0); // this is a dummy, we are not interested in the data
      }

      /*! unpack data from message buffer to user

        n is the number of objects sent by the sender
      */
      template<class MessageBuffer, class EntityType>
      void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
      {
        DataType x;
        buff.read(x); // this is a dummy, we are not interested in the data

        IndexType myindex = indexSet.index(e);
        IdType myid = globalIdSet.id(e);
        l2g[myindex] = myid;
        g2l[myid] = myindex;
        //std::cout << gv.comm().rank() << ": scatter local=" << myindex << " global=" << myid << std::endl;
      }

      //! constructor
      LocalGlobalMapDataHandle (const GFS& gfs_, LocalToGlobalMap& l2g_, GlobalToLocalMap& g2l_)
        : gfs(gfs_), gv(gfs.gridView()), indexSet(gv.indexSet()),
          grid(gv.grid()), globalIdSet(grid.globalIdSet()),
          l2g(l2g_), g2l(g2l_)
      {
      }

    private:
      const GFS& gfs;
      GV gv;
      const IndexSet& indexSet;
      const Grid& grid;
      const GlobalIdSet& globalIdSet;
      LocalToGlobalMap& l2g;
      GlobalToLocalMap& g2l;
    };


    // A DataHandle class to exchange rows of a sparse matrix
    template<class GFS, class M> // mapper type and vector type
    class MatrixExchangeDataHandle
      : public Dune::CommDataHandleIF<MatrixExchangeDataHandle<GFS,M>,
                                      std::pair<typename GFS::Traits::GridView::Grid::Traits::GlobalIdSet::IdType,
                                                typename M::block_type> >
    {
    public:
      // some types from the grid view
      typedef typename GFS::Traits::GridView GV;
      typedef typename GV::IndexSet IndexSet;
      typedef typename IndexSet::IndexType IndexType;
      typedef typename GV::Grid Grid;
      typedef typename Grid::Traits::GlobalIdSet GlobalIdSet;
      typedef typename GlobalIdSet::IdType IdType;

      // some types from the matrix
      typedef typename M::block_type B;

      //! export type of data for message buffer
      typedef std::pair<IdType,B> DataType;

      // map types
      typedef std::map<IndexType,IdType> LocalToGlobalMap;
      typedef std::map<IdType,IndexType> GlobalToLocalMap;

      //! returns true if data for this codim should be communicated
      bool contains (int dim, int codim) const
      {
        return (codim==0);
      }

      //! returns true if size per entity of given dim and codim is a constant
      bool fixedsize (int dim, int codim) const
      {
        return false;
      }

      /*! how many objects of type DataType have to be sent for a given entity

        Note: Only the sender side needs to know this size.
      */
      template<class EntityType>
      size_t size (EntityType& e) const
      {
        IndexType myindex = indexSet.index(e);
        typename M::row_type myrow = m[myindex];
        typename M::row_type::iterator endit=myrow.end();
        size_t count=0;
        for (typename M::row_type::iterator it=myrow.begin(); it!=endit; ++it)
          {
            typename LocalToGlobalMap::const_iterator find=l2g.find(it.index());
            if (find!=l2g.end())
              {
                count++;
              }
          }
        //std::cout << gv.comm().rank() << ": row=" << myindex << " size=" << count << std::endl;
        return count;
      }

      //! pack data from user to message buffer
      template<class MessageBuffer, class EntityType>
      void gather (MessageBuffer& buff, const EntityType& e) const
      {
        IndexType myindex = indexSet.index(e);
        //std::cout << gv.comm().rank() << ": begin gather local=" << myindex << std::endl;
        typename M::row_type myrow = m[myindex];
        typename M::row_type::iterator endit=myrow.end();
        size_t count = 0;
        for (typename M::row_type::iterator it=myrow.begin(); it!=endit; ++it)
          {
            typename LocalToGlobalMap::const_iterator find=l2g.find(it.index());
            if (find!=l2g.end())
              {
                buff.write(std::make_pair(find->second,*it));
                //std::cout << gv.comm().rank() << ":   ==> local=" << find->first << " global=" << find->second << std::endl;
                count++;
              }
          }
        //std::cout << gv.comm().rank() << ": end gather row=" << myindex << " size=" << count << std::endl;
      }

      /*! unpack data from message buffer to user

        n is the number of objects sent by the sender
      */
      template<class MessageBuffer, class EntityType>
      void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
      {
        IndexType myindex = indexSet.index(e);
        std::cout << gv.comm().rank() << ": begin scatter local=" << myindex << " size=" << n << std::endl;
        DataType x;
        size_t count = 0;
        for (size_t i=0; i<n; ++i)
          {
            buff.read(x);
            std::cout << gv.comm().rank() << ":   --> received global=" << x.first << std::endl;
            typename GlobalToLocalMap::const_iterator find=g2l.find(x.first);
            if (find!=g2l.end())
              {
                IndexType nbindex=find->second;
                if (m.exists(myindex,nbindex))
                  {
                    m[myindex][nbindex] = x.second;
                    B block(x.second);
                    block -= m2[myindex][nbindex];
                    std::cout << gv.comm().rank() << ":   compare i=" << myindex << " j=" << nbindex
                              << " norm=" << block.infinity_norm() << std::endl;

                    count++;
                    //std::cout << gv.comm().rank() << ":   --> local=" << find->first << " global=" << find->second << std::endl;
                  }
              }
          }
        //std::cout << gv.comm().rank() << ": end scatter row=" << myindex << " size=" << count << std::endl;
      }

      //! constructor
      MatrixExchangeDataHandle (const GFS& gfs_, M& m_, const LocalToGlobalMap& l2g_, const GlobalToLocalMap& g2l_,M& m2_)
        : gfs(gfs_), m(m_), gv(gfs.gridView()), indexSet(gv.indexSet()),
          grid(gv.grid()), globalIdSet(grid.globalIdSet()),
          l2g(l2g_), g2l(g2l_), m2(m2_)
      {
      }

    private:
      const GFS& gfs;
      M& m;
      GV gv;
      const IndexSet& indexSet;
      const Grid& grid;
      const GlobalIdSet& globalIdSet;
      const LocalToGlobalMap& l2g;
      const GlobalToLocalMap& g2l;
      M& m2;
    };


  /** A function to communicate matrix entries
   */
  template<class GFS, class T, class A, int n, int m>
  void restore_overlap_entries (const GFS& gfs, Dune::BCRSMatrix<Dune::FieldMatrix<T,n,m>,A>& matrix,
                                Dune::BCRSMatrix<Dune::FieldMatrix<T,n,m>,A>& matrix2)
  {
    typedef Dune::FieldMatrix<T,n,m> B;
    typedef Dune::BCRSMatrix<B,A> M; // m is of type M
    typedef typename LocalGlobalMapDataHandle<GFS>::LocalToGlobalMap LocalToGlobalMap;
    typedef typename LocalGlobalMapDataHandle<GFS>::GlobalToLocalMap GlobalToLocalMap;

    // build up two maps to associate local indices and global ids
    LocalToGlobalMap l2g;
    GlobalToLocalMap g2l;
    LocalGlobalMapDataHandle<GFS> lgdh(gfs,l2g,g2l);
    if (gfs.gridView().comm().size()>1)
      gfs.gridView().communicate(lgdh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);

    // exchange matrix entries
    MatrixExchangeDataHandle<GFS,M> mexdh(gfs,matrix,l2g,g2l,matrix2);
    if (gfs.gridView().comm().size()>1)
      gfs.gridView().communicate(mexdh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
  }

  //***********************************************************
  // The DG/AMG preconditioner in the overlapping case
  //***********************************************************

  /** An ISTL preconditioner for DG based on AMG applied to CG subspace

      The template parameters are:
      DGGFS       DG space
      DGMatrix    BCRSMatrix assembled with DG
      DGPrec      preconditioner to be used for DG
      CGPrec      preconditioner to be used on CG subspace
      P           BCRSMatrix for grid transfer
  */
  template<class DGGFS, class DGMatrix, class DGPrec, class DGCC,
           class CGGFS, class CGPrec, class CGCC,
           class P, class DGHelper, class Comm>
  class OvlpDGAMGPrec
    : public Dune::Preconditioner<Dune::PDELab::Backend::Vector<DGGFS,typename DGPrec::domain_type::field_type>,
                                  Dune::PDELab::Backend::Vector<DGGFS,typename DGPrec::range_type::field_type>>
  {
    const DGGFS& dggfs;
    DGMatrix& dgmatrix;
    DGPrec& dgprec;
    const DGCC& dgcc;
    const CGGFS& cggfs;
    CGPrec& cgprec;
    const CGCC& cgcc;
    P& p;
    const DGHelper& dghelper;
    const Comm& comm;
    int n1,n2;

  public:
    using V = Dune::PDELab::Backend::Vector<DGGFS,typename DGPrec::domain_type::field_type>;
    using W = Dune::PDELab::Backend::Vector<DGGFS,typename DGPrec::range_type::field_type>;
    typedef typename V::BaseT X;
    typedef typename W::BaseT Y;
    using CGV = Dune::PDELab::Backend::Vector<CGGFS,typename CGPrec::domain_type::field_type>;
    using CGW = Dune::PDELab::Backend::Vector<CGGFS,typename CGPrec::range_type::field_type>;

    // define the category
    enum {
      //! \brief The category the preconditioner is part of.
      category=Dune::SolverCategory::overlapping
    };

    /*! \brief Constructor.

      Constructor gets all parameters to operate the prec.
      \param A The matrix to operate on.
      \param n The number of iterations to perform.
      \param w The relaxation factor.
    */
    OvlpDGAMGPrec (const DGGFS& dggfs_, DGMatrix& dgmatrix_, DGPrec& dgprec_, const DGCC& dgcc_,
                   const CGGFS& cggfs_, CGPrec& cgprec_, const CGCC& cgcc_, P& p_,
                   const DGHelper& dghelper_, const Comm& comm_, int n1_, int n2_)
      : dggfs(dggfs_), dgmatrix(dgmatrix_), dgprec(dgprec_), dgcc(dgcc_),
        cggfs(cggfs_), cgprec(cgprec_), cgcc(cgcc_), p(p_), dghelper(dghelper_),
        comm(comm_), n1(n1_), n2(n2_)
    {
    }

    /*!
      \brief Prepare the preconditioner.

      \copydoc Preconditioner::pre(X&,Y&)
    */
    virtual void pre (V& x, W& b)
    {
      using Backend::native;
      dgprec.pre(native(x),native(b));
      CGW cgd(cggfs,0.0);
      CGV cgv(cggfs,0.0);
      cgprec.pre(native(cgv),native(cgd));
    }

    /*!
      \brief Apply the precondioner.

      \copydoc Preconditioner::apply(X&,const Y&)
    */
    virtual void apply (V& x, const W& b)
    {
      using Backend::native;
      // need local copies to store defect and solution
      W d(b);
      Dune::PDELab::set_constrained_dofs(dgcc,0.0,d);
      V v(x); // only to get correct size

      // pre-smoothing on DG matrix
      for (int i=0; i<n1; i++)
        {
          using Backend::native;
          v = 0.0;
          dgprec.apply(native(v),native(d));
          Dune::PDELab::AddDataHandle<DGGFS,V> adddh(dggfs,v);
          if (dggfs.gridView().comm().size()>1)
            dggfs.gridView().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);
          dgmatrix.mmv(native(v),native(d));
          Dune::PDELab::set_constrained_dofs(dgcc,0.0,d);
          x += v;
        }

      // restrict defect to CG subspace
      dghelper.maskForeignDOFs(d); // DG defect is additive for overlap 1, but in case we use more
      CGW cgd(cggfs,0.0);
      p.mtv(native(d),native(cgd));
      Dune::PDELab::AddDataHandle<CGGFS,CGW> adddh(cggfs,cgd);
      if (cggfs.gridView().comm().size()>1)
        cggfs.gridView().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication); // now we have consistent defect on coarse grid
      Dune::PDELab::set_constrained_dofs(cgcc,0.0,cgd);
      comm.project(native(cgd));
      CGV cgv(cggfs,0.0);


      // call preconditioner
      cgprec.apply(native(cgv),native(cgd));

      // prolongate correction
      p.mv(native(cgv),native(v));
      dgmatrix.mmv(native(v),native(d));
      Dune::PDELab::set_constrained_dofs(dgcc,0.0,d);
      x += v;

      // post-smoothing on DG matrix
      for (int i=0; i<n2; i++)
        {
          v = 0.0;
          dgprec.apply(native(v),native(d));
          Dune::PDELab::AddDataHandle<DGGFS,V> adddh(dggfs,v);
          if (dggfs.gridView().comm().size()>1)
            dggfs.gridView().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);
          dgmatrix.mmv(native(v),native(d));
          Dune::PDELab::set_constrained_dofs(dgcc,0.0,d);
          x += v;
        }
    }

    /*!
      \brief Clean up.

      \copydoc Preconditioner::post(X&)
    */
    virtual void post (V& x)
    {
      dgprec.post(Backend::native(x));
      CGV cgv(cggfs,0.0);
      cgprec.post(Backend::native(cgv));
    }
  };

//***********************************************************
// The DG/AMG linear solver backend in the overlapping case
//***********************************************************

/** Overlapping solver backend for using AMG for DG in PDELab

    The template parameters are:
    DGGO         GridOperator for DG discretization, allows access to matrix, vector and grid function space
    DGCC         constraints container for DG problem
    CGGFS        grid function space for CG subspace
    CGCC         constraints container for CG problem
    TransferLOP  local operator to assemble prolongation from CGGFS to DGGFS
    DGPrec       preconditioner for DG problem
    Solver       solver to be used on the complete problem
    int s        size of global index to be used in AMG
*/
template<class DGGO, class DGCC, class CGGFS, class CGCC, class TransferLOP,
         template<class,class,class,int> class DGPrec, template<class> class Solver, int s=96>
class ISTLBackend_OVLP_AMG_4_DG :
  public Dune::PDELab::OVLPScalarProductImplementation<typename DGGO::Traits::TrialGridFunctionSpace>,
  public Dune::PDELab::LinearResultStorage
{
public:
  // DG grid function space
  typedef typename DGGO::Traits::TrialGridFunctionSpace GFS;

  // vectors and matrices on DG level
  typedef typename DGGO::Traits::Jacobian M; // wrapped istl DG matrix
  typedef typename DGGO::Traits::Domain V;   // wrapped istl DG vector
  typedef typename M::BaseT Matrix;          // istl DG matrix
  typedef typename V::BaseT Vector;          // istl DG vector
  typedef typename Vector::field_type field_type;

  // vectors and matrices on CG level
  using CGV = Dune::PDELab::Backend::Vector<CGGFS,field_type>; // wrapped istl CG vector
  typedef typename CGV::BaseT CGVector;                               // istl CG vector

  // prolongation matrix
  typedef Dune::PDELab::ISTLMatrixBackend MBE;
  typedef Dune::PDELab::EmptyTransformation CC;
  typedef TransferLOP CGTODGLOP; // local operator
  typedef Dune::PDELab::GridOperator<CGGFS,GFS,CGTODGLOP,MBE,field_type,field_type,field_type,CC,CC> PGO;
  typedef typename PGO::Jacobian PMatrix; // wrapped ISTL prolongation matrix
  typedef typename PMatrix::BaseT P;      // ISTL prolongation matrix

  // CG subspace matrix
  typedef typename Dune::TransposedMatMultMatResult<P,Matrix>::type PTADG;
  typedef typename Dune::MatMultMatResult<PTADG,P>::type CGMatrix; // istl coarse space matrix

private:

  const GFS& gfs;
  DGGO& dggo;
  const DGCC& dgcc;
  CGGFS& cggfs;
  const CGCC& cgcc;
  unsigned maxiter;
  int verbose;
  bool usesuperlu;

  CGTODGLOP cgtodglop;  // local operator to assemble prolongation matrix
  PGO pgo;              // grid operator to assemble prolongation matrix
  PMatrix pmatrix;      // wrapped prolongation matrix

  /** an empty local operator to assemble processor boundary constraints
   */
  class EmptyLop : public Dune::PDELab::NumericalJacobianApplyVolume<EmptyLop>,
                   public Dune::PDELab::FullVolumePattern,
                   public Dune::PDELab::LocalOperatorDefaultFlags,
                   public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  {
  };

public:

  // access to prolongation matrix
  PMatrix& prolongation_matrix ()
  {
    return pmatrix;
  }

  /** make backend object
   */
  ISTLBackend_OVLP_AMG_4_DG(DGGO& dggo_, const DGCC& dgcc_, CGGFS& cggfs_, const CGCC& cgcc_,
                            unsigned maxiter_=5000, int verbose_=1, bool usesuperlu_=true) :
    Dune::PDELab::OVLPScalarProductImplementation<typename DGGO::Traits::TrialGridFunctionSpace>(dggo_.trialGridFunctionSpace()),
    gfs(dggo_.trialGridFunctionSpace()), dggo(dggo_), dgcc(dgcc_), cggfs(cggfs_), cgcc(cgcc_), maxiter(maxiter_), verbose(verbose_), usesuperlu(usesuperlu_),
    cgtodglop(), pgo(cggfs,dggo.trialGridFunctionSpace(),cgtodglop), pmatrix(pgo)
  {
#if !HAVE_SUPERLU
    if (usesuperlu == true)
      {
        if (gfs.gridView().comm().rank()==0)
          std::cout << "WARNING: You are using AMG without SuperLU!"
                    << " Please consider installing SuperLU,"
                    << " or set the usesuperlu flag to false"
                    << " to suppress this warning." << std::endl;
      }
#endif

    // assemble prolongation matrix; this will not change from one apply to the next
    pmatrix = 0.0;
    if (verbose>0 && gfs.gridView().comm().rank()==0) std::cout << "allocated prolongation matrix of size " << pmatrix.N() << " x " << pmatrix.M() << std::endl;
    CGV cgx(cggfs,0.0);         // need vector to call jacobian
    pgo.jacobian(cgx,pmatrix);
  }

  /*! \brief solve the given linear system

    \param[in] A the given matrix
    \param[out] z the solution vector to be computed
    \param[in] r right hand side
    \param[in] reduction to be achieved
  */
  void apply (M& A, V& z, V& r, typename V::ElementType reduction)
  {
    using Backend::native;
    // make operator and scalar product for overlapping solver
    typedef Dune::PDELab::OverlappingOperator<DGCC,M,V,V> POP;
    POP pop(dgcc,A);
    typedef Dune::PDELab::OVLPScalarProduct<GFS,V> PSP;
    PSP psp(*this);

    // compute CG matrix
    // make grid operator with empty local operator => matrix data type and constraints assembly
    EmptyLop emptylop;
    typedef Dune::PDELab::GridOperator<CGGFS,CGGFS,EmptyLop,MBE,field_type,field_type,field_type,CGCC,CGCC> CGGO;
    CGGO cggo(cggfs,cgcc,cggfs,cgcc,emptylop);
    typedef typename CGGO::Jacobian CGM;

    // do triple matrix product ACG = P^T ADG P; this is purely local
    Dune::Timer watch;
    watch.reset();
    tags::attached_container attached_container;
    CGM acg(attached_container);
    {
      PTADG ptadg;
      Dune::transposeMatMultMat(ptadg,native(pmatrix),native(A)); // 1a
      //Dune::transposeMatMultMat(ptadg,native(pmatrix),native(A2));   // 1b
      Dune::matMultMat(native(acg),ptadg,native(pmatrix));
    }
    double triple_product_time = watch.elapsed();
    if (verbose>0 && gfs.gridView().comm().rank()==0) std::cout << "=== triple matrix product " << triple_product_time << " s" << std::endl;
    //Dune::printmatrix(std::cout,native(acg),"triple product matrix","row",10,2);
    CGV cgx(cggfs,0.0);     // need vector to call jacobian
    cggo.jacobian(cgx,acg); // insert trivial rows at processor boundaries
    //std::cout << "CG constraints: " << cgcc.size() << " out of " << cggfs.globalSize() << std::endl;

    // NOW we need to insert the processor boundary conditions in DG matrix
    typedef Dune::PDELab::GridOperator<GFS,GFS,EmptyLop,MBE,field_type,field_type,field_type,DGCC,DGCC> DGGOEmpty;
    DGGOEmpty dggoempty(gfs,dgcc,gfs,dgcc,emptylop);
    dggoempty.jacobian(z,A);

    // and in the residual
    Dune::PDELab::set_constrained_dofs(dgcc,0.0,r);

    // now set up parallel AMG solver for the CG subspace
    typedef typename Dune::PDELab::istl::CommSelector<s,Dune::MPIHelper::isFake>::type Comm;
    Comm oocc(gfs.gridView().comm());
    typedef Dune::PDELab::istl::ParallelHelper<CGGFS> CGHELPER;
    CGHELPER cghelper(cggfs,2);
    cghelper.createIndexSetAndProjectForAMG(acg,oocc);
    typedef Dune::OverlappingSchwarzOperator<CGMatrix,CGVector,CGVector,Comm> ParCGOperator;
    ParCGOperator paroop(native(acg),oocc);
    Dune::OverlappingSchwarzScalarProduct<CGVector,Comm> sp(oocc);
    typedef Dune::Amg::Parameters Parameters; // AMG parameters (might be nice to change from outside)
    Parameters params(15,2000);
    params.setDefaultValuesIsotropic(CGGFS::Traits::GridViewType::Traits::Grid::dimension);
    params.setDebugLevel(verbose);
    params.setCoarsenTarget(2000);
    params.setMaxLevel(20);
    params.setProlongationDampingFactor(1.6);
    params.setNoPreSmoothSteps(3);
    params.setNoPostSmoothSteps(3);
    params.setGamma(1);
    params.setAdditive(false);
    //params.setAccumulate(Dune::Amg::AccumulationMode::noAccu); // atOnceAccu results in deadlock
    typedef Dune::SeqSSOR<CGMatrix,CGVector,CGVector,1> Smoother;
    typedef Dune::BlockPreconditioner<CGVector,CGVector,Comm,Smoother> ParSmoother;
    typedef typename Dune::Amg::SmootherTraits<ParSmoother>::Arguments SmootherArgs;
    SmootherArgs smootherArgs;
    smootherArgs.iterations = 2;
    smootherArgs.relaxationFactor = 0.92;
    typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<CGMatrix,Dune::Amg::FirstDiagonal> > Criterion;
    Criterion criterion(params);
    typedef Dune::Amg::AMG<ParCGOperator,CGVector,ParSmoother,Comm> AMG;
    watch.reset();
    AMG amg(paroop,criterion,smootherArgs,oocc);
    double amg_setup_time = watch.elapsed();
    if (verbose>0 && gfs.gridView().comm().rank()==0) std::cout << "=== AMG setup " <<amg_setup_time << " s" << std::endl;

    // set up hybrid DG/CG preconditioner
    typedef DGPrec<Matrix,Vector,Vector,1> DGPrecType;
    DGPrecType dgprec(native(A),1,0.92);
    //DGPrecType dgprec(native(A),0.92);
    typedef Dune::PDELab::istl::ParallelHelper<GFS> DGHELPER;
    typedef OvlpDGAMGPrec<GFS,Matrix,DGPrecType,DGCC,CGGFS,AMG,CGCC,P,DGHELPER,Comm> HybridPrec;
    HybridPrec hybridprec(gfs,native(A),dgprec,dgcc,cggfs,amg,cgcc,native(pmatrix),
                          this->parallelHelper(),oocc,3,3);

    // /********/
    // /* Test */
    // /********/
    // Solver<CGVector> testsolver(paroop,sp,amg,1e-8,100,2);
    // CGV cgxx(cggfs,0.0);
    // CGV cgdd(cggfs,1.0);
    // Dune::InverseOperatorResult statstat;
    // testsolver.apply(native(cgxx),native(cgdd),statstat);
    // /********/

    // set up solver
    int verb=verbose;
    if (gfs.gridView().comm().rank()>0) verb=0;
    Solver<V> solver(pop,psp,hybridprec,reduction,maxiter,verb);

    // solve
    Dune::InverseOperatorResult stat;
    watch.reset();
    solver.apply(z,r,stat);
    double amg_solve_time = watch.elapsed();
    if (verbose>0 && gfs.gridView().comm().rank()==0) std::cout << "=== Hybrid total solve time " << amg_solve_time+amg_setup_time+triple_product_time << " s" << std::endl;
    res.converged  = stat.converged;
    res.iterations = stat.iterations;
    res.elapsed    = amg_solve_time+amg_setup_time+triple_product_time;
    res.reduction  = stat.reduction;
    res.conv_rate  = stat.conv_rate;
  }

};
}
}
#endif
