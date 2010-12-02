// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_NOVLPISTLSOLVERBACKEND_HH
#define DUNE_NOVLPISTLSOLVERBACKEND_HH

#include <cstddef>

#include <dune/common/deprecated.hh>
#include <dune/common/mpihelper.hh>

#include <dune/grid/common/gridenums.hh>

#include <dune/istl/io.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>

#include "../gridfunctionspace/constraints.hh"
#include "../gridfunctionspace/genericdatahandle.hh"
#include "../newton/newton.hh"
#include "istlvectorbackend.hh"
#include "parallelistlhelper.hh"

namespace Dune {
  namespace PDELab {

    //========================================================
    // Generic support for nonoverlapping grids
    //========================================================

    //! Operator for the non-overlapping parallel case
    /**
     * Calculate \f$y:=Ax\f$.
     *
     * \tparam GFS The GridFunctionSpace the vectors apply to.
     * \tparam M   Type of the matrix.  Should be one of the ISTL matrix types.
     * \tparam X   Type of the vectors the matrix is applied to.
     * \tparam Y   Type of the result vectors.
     */
    template<class GFS, class M, class X, class Y>
    class NonoverlappingOperator : public Dune::AssembledLinearOperator<M,X,Y>
    {
    public:
      //! export type of matrix
      typedef M matrix_type;
      //! export type of vectors the matrix is applied to
      typedef X domain_type;
      //! export type of result vectors
      typedef Y range_type;
      //! export type of the entries for x
      typedef typename X::field_type field_type;

      //redefine the category, that is the only difference
      enum {category=Dune::SolverCategory::nonoverlapping};

      //! Construct a non-overlapping operator
      /**
       * \param gfs_    GridFunctionsSpace for the vectors.
       * \param A       Matrix for this operator.  This should be the locally
       *                assembled matrix.
       * \param helper_ Helper for parallel communication (not used).
       *
       * \note The constructed object stores references to all the objects
       *       given as parameters here.  They should be valid for as long as
       *       the constructed object is used.  They are not needed to
       *       destruct the constructed object.
       *
       * \deprecated The helper_ parameter is unused.  Use the constructor
       *             without the helper_ parameter instead.
       */
      NonoverlappingOperator (const GFS& gfs_, const M& A,
                              const ParallelISTLHelper<GFS>& helper_)
        DUNE_DEPRECATED
        : gfs(gfs_), _A_(A)
      {
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
      NonoverlappingOperator (const GFS& gfs_, const M& A)
        : gfs(gfs_), _A_(A)
      { }

      //! apply operator
      /**
       * Compute \f$y:=A(x)\f$ on this process, then make y consistent (sum up
       * corresponding entries of y on the different processes and store the
       * result back in y on each process).
       */
      virtual void apply (const X& x, Y& y) const
      {
        // apply local operator; now we have sum y_p = sequential y
        _A_.mv(x,y);

        // accumulate y on border
        Dune::PDELab::AddDataHandle<GFS,Y> adddh(gfs,y);
        if (gfs.gridview().comm().size()>1)
          gfs.gridview().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
      }

      //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
      /**
       * Compute \f$y:=\alpha A(x)\f$ on this process, then make y consistent
       * (sum up corresponding entries of y on the different processes and
       * store the result back in y on each process).
       */
      virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
      {
        // apply local operator; now we have sum y_p = sequential y
        _A_.usmv(alpha,x,y);

        // accumulate y on border
        Dune::PDELab::AddDataHandle<GFS,Y> adddh(gfs,y);
        if (gfs.gridview().comm().size()>1)
          gfs.gridview().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
      }

      //! extract the matrix
      virtual const M& getmat () const
      {
        return _A_;
      }

    private:
      const GFS& gfs;
      const M& _A_;
    };


    // parallel scalar product assuming no overlap
    template<class GFS, class X>
    class NonoverlappingScalarProduct : public Dune::ScalarProduct<X>
    {
    public:
      //! export types
      typedef X domain_type;
      typedef typename X::ElementType field_type;

      //! define the category
      enum {category=Dune::SolverCategory::nonoverlapping};

      /*! \brief Constructor needs to know the grid function space
       */
      NonoverlappingScalarProduct (const GFS& gfs_, const ParallelISTLHelper<GFS>& helper_)
        : gfs(gfs_), helper(helper_)
      {}

      /*! \brief Dot product of two vectors.
        It is assumed that the vectors are consistent on the interior+border
        partition.
      */
      virtual field_type dot (const X& x, const X& y)
      {
        // do local scalar product on unique partition
        field_type sum = 0;
        for (typename X::size_type i=0; i<x.N(); ++i)
          for (typename X::size_type j=0; j<x[i].N(); ++j)
            sum += (x[i][j]*y[i][j])*helper.mask(i,j);

        // do global communication
        return gfs.gridview().comm().sum(sum);
      }

      /*! \brief Norm of a right-hand side vector.
        The vector must be consistent on the interior+border partition
      */
      virtual double norm (const X& x)
      {
        return sqrt(static_cast<double>(this->dot(x,x)));
      }

      /*! \brief make additive vector consistent
       */
      void make_consistent (X& x) const
      {
        Dune::PDELab::AddDataHandle<GFS,X> adddh(gfs,x);
        if (gfs.gridview().comm().size()>1)
          gfs.gridview().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
      }

    private:
      const GFS& gfs;
      const ParallelISTLHelper<GFS>& helper;
    };

    // parallel Richardson preconditioner
    template<class GFS, class X, class Y>
    class NonoverlappingRichardson : public Dune::Preconditioner<X,Y>
    {
    public:
      //! \brief The domain type of the preconditioner.
      typedef X domain_type;
      //! \brief The range type of the preconditioner.
      typedef Y range_type;
      //! \brief The field type of the preconditioner.
      typedef typename X::ElementType field_type;

      // define the category
      enum {
        //! \brief The category the preconditioner is part of.
        category=Dune::SolverCategory::nonoverlapping
      };

      //! \brief Constructor.
      NonoverlappingRichardson (const GFS& gfs_, const ParallelISTLHelper<GFS>& helper_)
        : gfs(gfs_), helper(helper_)
      {
      }

      /*!
        \brief Prepare the preconditioner.

        \copydoc Preconditioner::pre(X&,Y&)
      */
      virtual void pre (X& x, Y& b) {}

      /*!
        \brief Apply the precondioner.

        \copydoc Preconditioner::apply(X&,const Y&)
      */
      virtual void apply (X& v, const Y& d)
      {
        v = d;
        // no communication is necessary here because defect is already consistent!
//         helper.mask(v);
//         Dune::PDELab::AddDataHandle<GFS,X> adddh(gfs,v);
//         if (gfs.gridview().comm().size()>1)
//           gfs.gridview().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
      }

      /*!
        \brief Clean up.

        \copydoc Preconditioner::post(X&)
      */
      virtual void post (X& x) {}

    private:
      const GFS& gfs;
      const ParallelISTLHelper<GFS>& helper;
    };

    //! parallel non-overlapping Jacobi preconditioner
    /**
     * \tparam Diagonal Vector type used to store the diagonal of the matrix
     * \tparam X        Vector type used to store the result of applying the
     *                  preconditioner.
     * \tparam Y        Vector type used to store the defect.
     *
     * The Jacobi preconditioner approximates the inverse of a matrix M by
     * taking the diagonal diag(M) and inverting that.  In the parallel case
     * the matrix M is assumed to be inconsistent, so diagonal entries for
     * dofs on the border are summed up over all relevant processes by this
     * precoditioner before the inverse is computed.
     */
    template<class Diagonal, class X, class Y>
    class NonoverlappingJacobi : public Dune::Preconditioner<X,Y>
    {
      typedef typename Diagonal::Backend DBackend;

      std::size_t gfsSize;
      Diagonal diagonal;

    public:
      //! The domain type of the operator.
      /**
       * The preconditioner is an inverse operator, so this is the output type
       * of the preconditioner.
       */
      typedef X domain_type;
      //! \brief The range type of the operator.
      /**
       * The preconditioner is an inverse operator, so this is the input type
       * of the preconditioner.
       */
      typedef Y range_type;
      //! \brief The field type of the preconditioner.
      typedef typename X::ElementType field_type;

      enum {
        //! \brief The category the preconditioner is part of.
        category=Dune::SolverCategory::nonoverlapping
      };

      //! \brief Constructor.
      /**
       * \param gfs The GridFunctionSpace the matrix and the vectors live on.
       * \param m   The matrix whose inverse the preconditioner should
       *            estimate.  m is assumed to be inconsistent (i.e. rows for
       *            dofs on the border only contain the contribution of the
       *            local process).
       *
       * The preconditioner does not store any reference to the gfs or the
       * matrix m.  The diagonal of m is copied, since it has to be made
       * consistent.
       */
      template<class GFS, class Matrix>
      NonoverlappingJacobi(const GFS& gfs, const Matrix &m) :
        gfsSize(gfs.size()), diagonal(gfs)
      {
        typedef typename Matrix::Backend MBackend;

        for(std::size_t i = 0; i < gfsSize; ++i)
          DBackend::access(diagonal, i) = MBackend::access(m, i, i);

        AddDataHandle<GFS, Diagonal> addDH(gfs, diagonal);
        gfs.gridview().communicate(addDH,
                                   InteriorBorder_InteriorBorder_Interface,
                                   ForwardCommunication);
      }

      //! Prepare the preconditioner.
      /**
       * \copydoc Preconditioner::pre(X&,Y&)
       */
      virtual void pre (X& x, Y& b) {}

      //! Apply the precondioner.
      /*
       * \copydoc Preconditioner::apply(X&,const Y&)
       *
       * For this preconditioner, this method works with both consistent and
       * inconsistent vectors: if d is consistent, v will be consistent, if d
       * is inconsistent, v will be inconsistent.
       */
      virtual void apply (X& v, const Y& d)
      {
        typedef typename X::Backend XBackend;
        typedef typename Y::Backend YBackend;

        for(std::size_t i = 0; i < gfsSize; ++i)
          XBackend::access(v, i) =
            YBackend::access(d, i) / DBackend::access(diagonal, i);
      }

      //! Clean up.
      /*
       * \copydoc Preconditioner::post(X&)
       */
      virtual void post (X& x) {}
    };

    template<class GFS>
    class ISTLBackend_NOVLP_CG_NOPREC
    {
      typedef Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;

    public:
      /*! \brief make a linear solver object

        \param[in] gfs a grid function space
        \param[in] maxiter maximum number of iterations to do
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_NOVLP_CG_NOPREC (const GFS& gfs_,
                                            unsigned maxiter_=5000,
                                            int verbose_=1)
        : gfs(gfs_), phelper(gfs), maxiter(maxiter_), verbose(verbose_)
      {}

      /*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      template<class V>
      typename V::ElementType norm (const V& v) const
      {
        V x(v); // make a copy because it has to be made consistent
        typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        psp.make_consistent(x);
        return psp.norm(x);
      }

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename V::ElementType reduction)
      {
        typedef Dune::PDELab::NonoverlappingOperator<GFS,M,V,W> POP;
        POP pop(gfs,A);
        typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        typedef Dune::PDELab::NonoverlappingRichardson<GFS,V,W> PRICH;
        PRICH prich(gfs,phelper);
        int verb=0;
        if (gfs.gridview().comm().rank()==0) verb=verbose;
        Dune::CGSolver<V> solver(pop,psp,prich,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
      }

      /*! \brief Return access to result data */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res;
      }

    private:
      const GFS& gfs;
      PHELPER phelper;
      Dune::PDELab::LinearSolverResult<double> res;
      unsigned maxiter;
      int verbose;
    };

    template<class GFS>
    class ISTLBackend_NOVLP_BCGS_NOPREC
    {
      typedef Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;

    public:
      /*! \brief make a linear solver object

        \param[in] gfs a grid function space
        \param[in] maxiter maximum number of iterations to do
        \param[in] verbose print messages if true
      */
      explicit ISTLBackend_NOVLP_BCGS_NOPREC (const GFS& gfs_, unsigned maxiter_=5000, int verbose_=1)
        : gfs(gfs_), phelper(gfs), maxiter(maxiter_), verbose(verbose_)
      {}

      /*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      template<class V>
      typename V::ElementType norm (const V& v) const
      {
        V x(v); // make a copy because it has to be made consistent
        typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        psp.make_consistent(x);
        return psp.norm(x);
      }

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename V::ElementType reduction)
      {
        typedef Dune::PDELab::NonoverlappingOperator<GFS,M,V,W> POP;
        POP pop(gfs,A);
        typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        typedef Dune::PDELab::NonoverlappingRichardson<GFS,V,W> PRICH;
        PRICH prich(gfs,phelper);
        int verb=0;
        if (gfs.gridview().comm().rank()==0) verb=verbose;
        Dune::BiCGSTABSolver<V> solver(pop,psp,prich,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
      }

      /*! \brief Return access to result data */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res;
      }

    private:
      const GFS& gfs;
      PHELPER phelper;
      Dune::PDELab::LinearSolverResult<double> res;
      unsigned maxiter;
      int verbose;
    };

    //! Solver to be used for explicit time-steppers with (block-)diagonal mass matrix
    template<typename GFS>
    class ISTLBackend_NOVLP_ExplicitDiagonal
    {
      typedef Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;

      const GFS& gfs;
      PHELPER phelper;
      Dune::PDELab::LinearSolverResult<double> res;

    public:
      /*! \brief make a linear solver object

        \param[in] gfs GridFunctionSpace, used to identify DoFs for parallel
        communication
      */
      explicit ISTLBackend_NOVLP_ExplicitDiagonal(const GFS& gfs_)
        : gfs(gfs_), phelper(gfs)
      {}

      /*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      template<class V>
      typename V::ElementType norm (const V& v) const
      {
        typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,V> PSP;
        V x(v); // make a copy because it has to be made consistent
        PSP psp(gfs,phelper);
        psp.make_consistent(x);
        return psp.norm(x);
      }

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename W::ElementType reduction)
      {
        Dune::SeqJac<M,V,W> jac(A,1,1.0);
        jac.pre(z,r);
        jac.apply(z,r);
        jac.post(z);
        if (gfs.gridview().comm().size()>1)
        {
          Dune::PDELab::AddDataHandle<GFS,V> adddh(gfs,z);
          gfs.gridview().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
        }
        res.converged  = true;
        res.iterations = 1;
        res.elapsed    = 0.0;
        res.reduction  = reduction;
      }

      /*! \brief Return access to result data */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res;
      }
    };


    /**
    * @brief Helper class for adding up matrix entries on border.
    * @tparam GridOperatorSpace The grid operator space to work on.
    * @tparam Scalar The type of the scalar matrix entries.
    */
    template<class GridOperatorSpace, class Scalar>
    class VertexExchanger
    {
      typedef typename GridOperatorSpace::Traits GridOperatorSpaceTraits;
      typedef typename GridOperatorSpaceTraits::GridViewType GridView;
      enum {dim = GridView::dimension};
      typedef typename GridView::Traits::Grid Grid;
      typedef typename GridOperatorSpace::template MatrixContainer<Scalar>::Type Matrix;
      typedef typename Matrix::block_type BlockType;
      typedef typename GridView::template Codim<dim>::Iterator  VertexIterator;
      typedef typename Grid::Traits::GlobalIdSet IDS;
      typedef typename IDS::IdType IdType;

    public:
      /*! \brief Constructor. Sets up the local to global relations.

      \param[in] gridView The grid view to operate on.
      */
      VertexExchanger(const GridView& gridView)
        : gridView_(gridView)
      {
        gid2Index_.clear();
        index2GID_.clear();


        VertexIterator vertexEndIt = gridView_.template end<dim>();
        for (VertexIterator vertexIt = gridView_.template begin<dim>(); vertexIt != vertexEndIt; ++vertexIt)
        {
          if (vertexIt->partitionType() == BorderEntity)
          {
            int localIdx = gridView_.indexSet().index(*vertexIt);
            IdType globalIdx = gridView_.grid().globalIdSet().id(*vertexIt);

            std::pair<IdType,int> g2iPair(globalIdx, localIdx);
            gid2Index_.insert(g2iPair);

            std::pair<int,IdType> i2gPair(localIdx, globalIdx);
            index2GID_.insert(i2gPair);

          }
        }
      }

      //! Local matrix blocks associated with the global id set
      struct MatEntry
      {
        IdType first;
        BlockType second;
        MatEntry (const IdType& f, const BlockType& s) : first(f),second(s) {}
        MatEntry () {}
      };

      //! A DataHandle class to exchange matrix entries
      class MatEntryExchange
      : public CommDataHandleIF<MatEntryExchange,MatEntry> {
        typedef typename Matrix::RowIterator RowIterator;
        typedef typename Matrix::ColIterator ColIterator;
      public:
        //! Export type of data for message buffer
        typedef MatEntry DataType;

        /** @brief Returns true if data for given valid codim should be communicated
        @copydoc CommDataHandleIF::contains(int, int)
        */
        bool contains (int dim, int codim) const
        {
          return (codim==dim);
        }

        /** @brief Returns true if size of data per entity of given dim and codim is a constant
        @copydoc CommDataHandleIF::fixedsize(int, int)
        */
        bool fixedsize (int dim, int codim) const
        {
          return false;
        }

        /** @brief How many objects of type DataType have to be sent for a given entity
        @copydoc CommDataHandleIF::size(EntityType&)
        */
        template<class EntityType>
        size_t size (EntityType& e) const
        {
          int i = gridView_.indexSet().index(e);
          int n = 0;
          for (ColIterator j = A_[i].begin(); j != A_[i].end(); ++j)
            {
            typename std::map<int,IdType>::const_iterator it = index2GID_.find(j.index());
            if (it != index2GID_.end())
              n++;
            }
          
            return n;
        }

        /** @brief Pack data from user to message buffer
        @copydoc CommDataHandleIF::gather(MessageBuffer&, const EntityType&)
        */
        template<class MessageBuffer, class EntityType>
        void gather (MessageBuffer& buff, const EntityType& e) const
        {
          int i = gridView_.indexSet().index(e);
          for (ColIterator j = A_[i].begin(); j != A_[i].end(); ++j)
            {
              typename std::map<int,IdType>::const_iterator it=index2GID_.find(j.index());
              if (it != index2GID_.end())
                buff.write(MatEntry(it->second,*j));
            }
        
        }

        /** @brief Unpack data from message buffer to user
        @copydoc CommDataHandleIF::scatter(MessageBuffer&, const EntityType&, size_t)
        */
        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
        {
          int i = gridView_.indexSet().index(e);
          for (size_t k = 0; k < n; k++)
          {
            MatEntry m;
            buff.read(m);
            // only add entries corresponding to border entities
            typename std::map<IdType,int>::const_iterator it = gid2Index_.find(m.first);
            if (it != gid2Index_.end())
              if (A_[i].find(it->second) != A_[i].end())
                A_[i][it->second] += m.second;
          }
        }

        /** @brief Constructor
        @param[in] gridView Grid view.
        @param[in] g2i Global to local index map.
        @param[in] i2g Local to global index map.
        @param[in] A Matrix to operate on.
        */
        MatEntryExchange (const GridView& gridView, const std::map<IdType,int>& g2i,
            const std::map<int,IdType>& i2g,
                          Matrix& A)
          : gridView_(gridView), gid2Index_(g2i), index2GID_(i2g), A_(A)
        {}

      private:
        const GridView& gridView_;
        const std::map<IdType,int>& gid2Index_;
        const std::map<int,IdType>& index2GID_;
        Matrix& A_;
      };

      /** @brief Sums up the entries corresponding to border vertices.
      @param A Matrix to operate on.
      */
      void sumEntries (Matrix& A)
      {
        if (gridView_.comm().size() > 1)
        {
          MatEntryExchange datahandle(gridView_, gid2Index_, index2GID_, A);
          gridView_.communicate(datahandle, InteriorBorder_InteriorBorder_Interface, ForwardCommunication);
        }
      }

      private:
      const GridView& gridView_;
      std::map<IdType,int> gid2Index_;
      std::map<int,IdType> index2GID_;
    };
      
    /**
    * @brief Nonoverlapping parallel BiCGSTAB solver preconditioned by block SSOR.
    * @tparam GridOperatorSpace The grid operator space to work on.
    * @tparam Scalar The type of the scalar matrix entries.
    * @tparam C The communication object.
    *
    * The solver uses a NonoverlappingWrappedPreconditioner with underlying
    * sequential SSOR preconditioner. The crucial step is to add up the matrix entries
    * corresponding to the border vertices on each process. This is achieved by
    * performing a VertexExchanger::sumEntries(Matrix&) before constructing the
    * sequential SSOR.
    *
    * Currently the solver is only working for nonoverlapping grids without ghosts 
    * (YaspGrid, UG with option setGhosts(false) (not yet checked in))
    */
    template<class GridOperatorSpace, class Scalar, class C>
    class ISTLBackend_NOVLP_BCGS_SSORk
    {
      typedef typename GridOperatorSpace::Traits::TrialGridFunctionSpace GridFunctionSpace;
      typedef Dune::PDELab::ParallelISTLHelper<GridFunctionSpace> PHELPER;
      typedef typename GridFunctionSpace::template VectorContainer<Scalar>::Type U;
      
    public:
      /*! \brief Constructor.

      \param[in] gfs a grid function space
      \param[in] maxiter maximum number of iterations to do
      \param[in] kssor iteration count for the SSOR preconditioner
      \param[in] verbose print messages if true
      */
      explicit ISTLBackend_NOVLP_BCGS_SSORk (const GridFunctionSpace& gfs, const C& c_, unsigned maxiter = 5000, unsigned kssor = 5, int verbose = 1)
        : gfs_(gfs), c(c_), phelper_(gfs_),
          maxiter_(maxiter), kssor_(kssor), verbose_(verbose), exchanger_(gfs_.gridview())
      {
      }

      /*! \brief Compute global norm of a vector.

      \param[in] v the given vector
      */
      template<class Vector>
     typename Vector::ElementType norm (const Vector& v) const
      {
        Vector x(v); // make a copy because it has to be made consistent
        typedef Dune::PDELab::NonoverlappingScalarProduct<GridFunctionSpace,Vector> PSP;
        PSP psp(gfs_,phelper_);
        psp.make_consistent(x);
        return psp.norm(x);
      }

      /*! \brief Solve the given linear system.

      \param[in] A the given matrix
      \param[out] z the solution vector to be computed
      \param[in] r right hand side
      \param[in] reduction to be achieved
      */
      template<class Matrix, class SolVector, class RhsVector>
      void apply(Matrix& A, SolVector& z, RhsVector& r, typename SolVector::ElementType reduction)
      { typedef typename Matrix::BaseT MatrixType;
        typedef typename BlockProcessor<GridFunctionSpace>::template AMGVectorTypeSelector<SolVector>::Type
          VectorType;
        typedef typename SolVector::ElementType Scalar;
        typedef typename GridOperatorSpace::template MatrixContainer<Scalar>::Type M;
        typedef Dune::SeqSSOR<MatrixType,VectorType,VectorType> SeqPreCond;
        exchanger_.sumEntries(A);
        MatrixType& mat=A.base();
        SeqPreCond seqPreCond(mat, kssor_, 1.0);
        typedef typename CommSelector<96,Dune::MPIHelper::isFake>::type Comm;
        Comm oocc(gfs_.gridview().comm());
        phelper_.createIndexSetAndProjectForAMG(mat, oocc);
        typedef Dune::NonoverlappingSchwarzScalarProduct<VectorType,Comm> PSP;
        PSP psp(oocc);      
        typedef Dune::NonoverlappingSchwarzOperator<MatrixType,VectorType,VectorType,Comm> POP;
        POP pop(mat,oocc);
        typedef Dune::NonoverlappingBlockPreconditioner<Comm, SeqPreCond> ParPreCond;
        ParPreCond parPreCond(seqPreCond, oocc);
        int verb=0;
        if (gfs_.gridview().comm().rank()==0) verb=verbose_;
        Dune::BiCGSTABSolver<VectorType> solver(pop,psp,parPreCond,reduction,maxiter_,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res_.converged  = stat.converged;
        res_.iterations = stat.iterations;
        res_.elapsed    = stat.elapsed;
        res_.reduction  = stat.reduction;
      }

      /*! \brief Return access to result data. */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res_;
      }

    private:
      const GridFunctionSpace& gfs_;
      const C& c;
      PHELPER phelper_;
      Dune::PDELab::LinearSolverResult<double> res_;
      unsigned maxiter_;
      unsigned kssor_;
      int verbose_;
      VertexExchanger<GridOperatorSpace, Scalar> exchanger_;
    };

  } // namespace PDELab
} // namespace Dune

#endif
