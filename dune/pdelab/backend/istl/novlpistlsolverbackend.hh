// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_NOVLPISTLSOLVERBACKEND_HH
#define DUNE_PDELAB_BACKEND_ISTL_NOVLPISTLSOLVERBACKEND_HH

// this is here for backwards compatibility and deprecation warnings, remove after 2.5.0
#include "ensureistlinclude.hh"

#include <cstddef>

#include <dune/common/deprecated.hh>
#include <dune/common/parallel/mpihelper.hh>

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

#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/backend/istl/vector.hh>
#include <dune/pdelab/backend/istl/bcrsmatrix.hh>
#include <dune/pdelab/backend/istl/blockmatrixdiagonal.hh>
#include <dune/pdelab/backend/istl/parallelhelper.hh>
#include <dune/pdelab/backend/istl/seqistlsolverbackend.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Backend
    //! \ingroup PDELab
    //! \{

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
    template<typename GFS, typename M, typename X, typename Y>
    class NonoverlappingOperator
      : public Dune::AssembledLinearOperator<M,X,Y>
    {
    public:
      //! export type of matrix
      using matrix_type = Backend::Native<M>;
      //! export type of vectors the matrix is applied to
      using domain_type = Backend::Native<X>;
      //! export type of result vectors
      using range_type = Backend::Native<Y>;
      //! export type of the entries for x
      typedef typename X::field_type field_type;

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
        using Backend::native;
        // apply local operator; now we have sum y_p = sequential y
        native(_A_).mv(native(x),native(y));

        // accumulate y on border
        Dune::PDELab::AddDataHandle<GFS,Y> adddh(gfs,y);
        if (gfs.gridView().comm().size()>1)
          gfs.gridView().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
      }

      //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
      /**
       * Compute \f$y:=\alpha A(x)\f$ on this process, then make y consistent
       * (sum up corresponding entries of y on the different processes and
       * store the result back in y on each process).
       */
      virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
      {
        using Backend::native;
        // apply local operator; now we have sum y_p = sequential y
        native(_A_).usmv(alpha,native(x),native(y));

        // accumulate y on border
        Dune::PDELab::AddDataHandle<GFS,Y> adddh(gfs,y);
        if (gfs.gridView().comm().size()>1)
          gfs.gridView().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
      }

      SolverCategory::Category category() const override
      {
        return SolverCategory::nonoverlapping;
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

      SolverCategory::Category category() const override
      {
        return SolverCategory::nonoverlapping;
      }

      /*! \brief Constructor needs to know the grid function space
       */
      NonoverlappingScalarProduct (const GFS& gfs_, const ISTL::ParallelHelper<GFS>& helper_)
        : gfs(gfs_), helper(helper_)
      {}

      /*! \brief Dot product of two vectors.
        It is assumed that the vectors are consistent on the interior+border
        partition.
      */
      virtual field_type dot (const X& x, const X& y)
      {
        // do local scalar product on unique partition
        field_type sum = helper.disjointDot(x,y);

        // do global communication
        return gfs.gridView().comm().sum(sum);
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
        if (gfs.gridView().comm().size()>1)
          gfs.gridView().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
      }

    private:
      const GFS& gfs;
      const ISTL::ParallelHelper<GFS>& helper;
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
      SolverCategory::Category category() const override
      {
        return SolverCategory::nonoverlapping;
      }

      //! \brief Constructor.
      NonoverlappingRichardson (const GFS& gfs_, const ISTL::ParallelHelper<GFS>& helper_)
        : gfs(gfs_), helper(helper_)
      {
      }

      /*!
        \brief Prepare the preconditioner.
      */
      virtual void pre (X& x, Y& b) {}

      /*!
        \brief Apply the precondioner.
      */
      virtual void apply (X& v, const Y& d)
      {
        v = d;
      }

      /*!
        \brief Clean up.
      */
      virtual void post (X& x) {}

    private:
      const GFS& gfs;
      const ISTL::ParallelHelper<GFS>& helper;
    };

    //! parallel non-overlapping Jacobi preconditioner
    /**
     * \tparam Diagonal Vector type used to store the diagonal of the matrix
     * \tparam X        Vector type used to store the result of applying the
     *                  preconditioner.
     * \tparam Y        Vector type used to store the defect.
     * \tparam A        The matrix type to be used.
     *
     * The Jacobi preconditioner approximates the inverse of a matrix M by
     * taking the diagonal diag(M) and inverting that.  In the parallel case
     * the matrix M is assumed to be inconsistent, so diagonal entries for
     * dofs on the border are summed up over all relevant processes by this
     * precoditioner before the inverse is computed.
     */
    template<typename A, typename X, typename Y>
    class NonoverlappingJacobi
      : public Dune::Preconditioner<X,Y>
    {

      typedef typename ISTL::BlockMatrixDiagonal<A>::MatrixElementVector Diagonal;

      Diagonal _inverse_diagonal;

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

      SolverCategory::Category category() const override
      {
        return SolverCategory::nonoverlapping;
      }

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
      template<typename GFS>
      NonoverlappingJacobi(const GFS& gfs, const A &m)
        : _inverse_diagonal(m)
      {
        // make the diagonal consistent...
        typename ISTL::BlockMatrixDiagonal<A>::template AddMatrixElementVectorDataHandle<GFS> addDH(gfs, _inverse_diagonal);
        gfs.gridView().communicate(addDH,
                                   InteriorBorder_InteriorBorder_Interface,
                                   ForwardCommunication);

        // ... and then invert it
        _inverse_diagonal.invert();

      }

      //! Prepare the preconditioner.
      virtual void pre (X& x, Y& b) {}

      //! Apply the precondioner.
      /*
       * For this preconditioner, this method works with both consistent and
       * inconsistent vectors: if d is consistent, v will be consistent, if d
       * is inconsistent, v will be inconsistent.
       */
      virtual void apply (X& v, const Y& d)
      {
        _inverse_diagonal.mv(d,v);
      }

      //! Clean up.
      virtual void post (X& x) {}
    };

    //! \addtogroup PDELab_novlpsolvers Nonoverlapping Solvers
    //! \{

    //! \brief Nonoverlapping parallel CG solver without preconditioner
    template<class GFS>
    class ISTLBackend_NOVLP_CG_NOPREC
    {
      typedef ISTL::ParallelHelper<GFS> PHELPER;

    public:
      /*! \brief make a linear solver object

        \param[in] gfs_ a grid function space
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_NOVLP_CG_NOPREC (const GFS& gfs_,
                                            unsigned maxiter_=5000,
                                            int verbose_=1)
        : gfs(gfs_), phelper(gfs,verbose_), maxiter(maxiter_), verbose(verbose_)
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
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename V::ElementType >::real_type reduction)
      {
        typedef Dune::PDELab::NonoverlappingOperator<GFS,M,V,W> POP;
        POP pop(gfs,A);
        typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        typedef Dune::PDELab::NonoverlappingRichardson<GFS,V,W> PRICH;
        PRICH prich(gfs,phelper);
        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        Dune::CGSolver<V> solver(pop,psp,prich,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
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

    //! \brief Nonoverlapping parallel CG solver with Jacobi preconditioner
    template<class GFS>
    class ISTLBackend_NOVLP_CG_Jacobi
    {
      typedef ISTL::ParallelHelper<GFS> PHELPER;

      const GFS& gfs;
      PHELPER phelper;
      LinearSolverResult<double> res;
      unsigned maxiter;
      int verbose;

    public:
      //! make a linear solver object
      /**
       * \param gfs_     A grid function space
       * \param maxiter_ Maximum number of iterations to do.
       * \param verbose_ Verbosity level, directly handed to the CGSolver.
       */
      explicit ISTLBackend_NOVLP_CG_Jacobi(const GFS& gfs_,
                                           unsigned maxiter_ = 5000,
                                           int verbose_ = 1) :
        gfs(gfs_), phelper(gfs,verbose_), maxiter(maxiter_), verbose(verbose_)
      {}

      //! compute global norm of a vector
      /**
       * \param v The vector to compute the norm of.  Should be an
       *          inconsistent vector (i.e. the entries corresponding a DoF on
       *          the border should only contain the summand of this process).
       */
      template<class V>
      typename V::ElementType norm (const V& v) const
      {
        V x(v); // make a copy because it has to be made consistent
        typedef NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        psp.make_consistent(x);
        return psp.norm(x);
      }

      //! solve the given linear system
      /**
       * \param A         The matrix to solve.  Should be a matrix from one of
       *                  PDELabs ISTL backends (only ISTLBCRSMatrixBackend at
       *                  the moment).
       * \param z         The solution vector to be computed
       * \param r         Right hand side
       * \param reduction to be achieved
       *
       * Solve the linear system A*z=r such that
       * norm(A*z0-r)/norm(A*z-r) < reduction where z0 is the initial value of
       * z.
       */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename V::ElementType >::real_type reduction)
      {
        typedef NonoverlappingOperator<GFS,M,V,W> POP;
        POP pop(gfs,A);
        typedef NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);

        typedef NonoverlappingJacobi<M,V,W> PPre;
        PPre ppre(gfs,Backend::native(A));

        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        CGSolver<V> solver(pop,psp,ppre,reduction,maxiter,verb);
        InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }

      //! Return access to result data
      const LinearSolverResult<double>& result() const
      { return res; }
    };

    //! \brief Nonoverlapping parallel BiCGStab solver without preconditioner
    template<class GFS>
    class ISTLBackend_NOVLP_BCGS_NOPREC
    {
      typedef ISTL::ParallelHelper<GFS> PHELPER;

    public:
      /*! \brief make a linear solver object

        \param[in] gfs_ a grid function space
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_NOVLP_BCGS_NOPREC (const GFS& gfs_, unsigned maxiter_=5000, int verbose_=1)
        : gfs(gfs_), phelper(gfs,verbose_), maxiter(maxiter_), verbose(verbose_)
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
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename V::ElementType >::real_type reduction)
      {
        typedef Dune::PDELab::NonoverlappingOperator<GFS,M,V,W> POP;
        POP pop(gfs,A);
        typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        typedef Dune::PDELab::NonoverlappingRichardson<GFS,V,W> PRICH;
        PRICH prich(gfs,phelper);
        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        Dune::BiCGSTABSolver<V> solver(pop,psp,prich,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
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


    //! \brief Nonoverlapping parallel BiCGStab solver with Jacobi preconditioner
    template<class GFS>
    class ISTLBackend_NOVLP_BCGS_Jacobi
    {
      typedef ISTL::ParallelHelper<GFS> PHELPER;

    public:
      /*! \brief make a linear solver object

        \param[in] gfs_ a grid function space
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_NOVLP_BCGS_Jacobi (const GFS& gfs_, unsigned maxiter_=5000, int verbose_=1)
        : gfs(gfs_), phelper(gfs,verbose_), maxiter(maxiter_), verbose(verbose_)
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
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename V::ElementType >::real_type reduction)
      {
        typedef Dune::PDELab::NonoverlappingOperator<GFS,M,V,W> POP;
        POP pop(gfs,A);
        typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);

        typedef NonoverlappingJacobi<M,V,W> PPre;
        PPre ppre(gfs,A);

        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        Dune::BiCGSTABSolver<V> solver(pop,psp,ppre,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
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
      typedef ISTL::ParallelHelper<GFS> PHELPER;

      const GFS& gfs;
      PHELPER phelper;
      Dune::PDELab::LinearSolverResult<double> res;

    public:
      /*! \brief make a linear solver object

        \param[in] gfs_ GridFunctionSpace, used to identify DoFs for parallel
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
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename W::ElementType >::real_type reduction)
      {
        Dune::SeqJac<M,V,W> jac(A,1,1.0);
        jac.pre(z,r);
        jac.apply(z,r);
        jac.post(z);
        if (gfs.gridView().comm().size()>1)
        {
          Dune::PDELab::AddDataHandle<GFS,V> adddh(gfs,z);
          gfs.gridView().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
        }
        res.converged  = true;
        res.iterations = 1;
        res.elapsed    = 0.0;
        res.reduction  = reduction;
        res.conv_rate  = reduction; // pow(reduction,1.0/1)
      }

      /*! \brief Return access to result data */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res;
      }
    };
    //! \} Nonoverlapping Solvers


    /*! \brief Utility base class for preconditioned novlp backends.
     * \tparam GO The type of the grid operator for the spatial discretization.
     *            This class will be used to adjust the discretization matrix.
     *            and extract the trial grid function space.
     * \tparam Preconditioner The type of preconditioner to use.
     * \tparam Solver The type of solver to use.
     */
    template<class GO,
             template<class,class,class,int> class Preconditioner,
             template<class> class Solver>
    class ISTLBackend_NOVLP_BASE_PREC
    {
      typedef typename GO::Traits::TrialGridFunctionSpace GFS;
      typedef ISTL::ParallelHelper<GFS> PHELPER;

    public:
      /*! \brief Constructor.

        \param[in] gfs_ a grid function space
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] steps_ number of preconditioner steps to apply as inner iteration
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_NOVLP_BASE_PREC (const GO& grid_operator, unsigned maxiter_ = 5000, unsigned steps_ = 5, int verbose_ = 1)
        : _grid_operator(grid_operator)
        , gfs(grid_operator.trialGridFunctionSpace())
        , phelper(gfs,verbose_)
        , maxiter(maxiter_)
        , steps(steps_)
        , verbose(verbose_)
      {}

      /*! \brief Compute global norm of a vector.

        \param[in] v the given vector
      */
      template<class Vector>
      typename Vector::ElementType norm (const Vector& v) const
      {
        Vector x(v); // make a copy because it has to be made consistent
        typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,Vector> PSP;
        PSP psp(gfs,phelper);
        psp.make_consistent(x);
        return psp.norm(x);
      }

      /*! \brief Solve the given linear system.

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename V::ElementType reduction)
      {
        using MatrixType = Backend::Native<M>;
        MatrixType& mat = Backend::native(A);
        using VectorType = Backend::Native<W>;
#if HAVE_MPI
        typedef typename ISTL::CommSelector<96,Dune::MPIHelper::isFake>::type Comm;
        _grid_operator.make_consistent(A);
        ISTL::assertParallelUG(gfs.gridView().comm());
        Comm oocc(gfs.gridView().comm(),Dune::SolverCategory::nonoverlapping);
        phelper.createIndexSetAndProjectForAMG(mat, oocc);
        typedef Preconditioner<MatrixType,VectorType,VectorType,1> Smoother;
        Smoother smoother(mat, steps, 1.0);
        typedef Dune::NonoverlappingSchwarzScalarProduct<VectorType,Comm> PSP;
        PSP psp(oocc);
        typedef Dune::NonoverlappingSchwarzOperator<MatrixType,VectorType,VectorType,Comm> Operator;
        Operator oop(mat,oocc);
        typedef Dune::NonoverlappingBlockPreconditioner<Comm, Smoother> ParSmoother;
        ParSmoother parsmoother(smoother, oocc);
#else
        typedef Preconditioner<MatrixType,VectorType,VectorType,1> ParSmoother;
        ParSmoother parsmoother(mat, steps, 1.0);
        typedef Dune::SeqScalarProduct<VectorType> PSP;
        PSP psp;
        typedef Dune::MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
        Operator oop(mat);
#endif
        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        Solver<VectorType> solver(oop,psp,parsmoother,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        //make r consistent
        if (gfs.gridView().comm().size()>1){
          Dune::PDELab::AddDataHandle<GFS,V> adddh(gfs,r);
          gfs.gridView().communicate(adddh,
                                     Dune::InteriorBorder_InteriorBorder_Interface,
                                     Dune::ForwardCommunication);
        }

        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }

      /*! \brief Return access to result data. */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res;
      }

    private:
      const GO& _grid_operator;
      const GFS& gfs;
      PHELPER phelper;
      Dune::PDELab::LinearSolverResult<double> res;
      unsigned maxiter;
      unsigned steps;
      int verbose;
    };

    //! \addtogroup PDELab_novlpsolvers Nonoverlapping Solvers
    //! \{

  /**
   * @brief Nonoverlapping parallel BiCGSTAB solver preconditioned by block SSOR.
   * @tparam GO The type of the grid operator used for the spatial discretization
   * (or the fakeGOTraits class for the old grid operator space). It is used
   * to adjust the discretization matrix and extract the trial grid function space.
   *
   * The solver uses a NonoverlappingBlockPreconditioner with underlying
   * sequential SSOR preconditioner. The crucial step is to add up the matrix entries
   * corresponding to the border vertices on each process. This is achieved by
   * performing a VertexExchanger::sumEntries(Matrix&) before constructing the
   * sequential SSOR.
   */
  template<class GO>
  class ISTLBackend_NOVLP_BCGS_SSORk
    : public ISTLBackend_NOVLP_BASE_PREC<GO,Dune::SeqSSOR, Dune::BiCGSTABSolver>
  {

  public:
    /*! \brief make a linear solver object

      \param[in] gfs_ a grid function space
      \param[in] maxiter_ maximum number of iterations to do
      \param[in] steps_ number of SSOR steps to apply as inner iteration
      \param[in] verbose_ print messages if true
    */
    explicit ISTLBackend_NOVLP_BCGS_SSORk (const GO& grid_operator, unsigned maxiter_=5000,
                                           int steps_=5, int verbose_=1)
      : ISTLBackend_NOVLP_BASE_PREC<GO,Dune::SeqSSOR, Dune::BiCGSTABSolver>(grid_operator, maxiter_, steps_, verbose_)
    {}
  };

  /**
   * @brief Nonoverlapping parallel CG solver preconditioned by block SSOR.
   * @tparam GO The type of the grid operator used for the spatial discretization.
   *            This class will be used to adjust the discretization matrix.
   *             and extract the trial grid function space.
   */
  template<class GO>
  class ISTLBackend_NOVLP_CG_SSORk
    : public ISTLBackend_NOVLP_BASE_PREC<GO,Dune::SeqSSOR, Dune::CGSolver>
  {

  public:
    /*! \brief make a linear solver object

      \param[in] gfs_ a grid function space
      \param[in] maxiter_ maximum number of iterations to do
      \param[in] steps_ number of SSOR steps to apply as inner iteration
      \param[in] verbose_ print messages if true
    */
    explicit ISTLBackend_NOVLP_CG_SSORk (const GO& grid_operator, unsigned maxiter_=5000,
                                         int steps_=5, int verbose_=1)
      : ISTLBackend_NOVLP_BASE_PREC<GO,Dune::SeqSSOR, Dune::CGSolver>(grid_operator, maxiter_, steps_, verbose_)
    {}
  };
    //! \} Nonoverlapping Solvers
    //! \} group Backend

    template<class GO,int s, template<class,class,class,int> class Preconditioner,
             template<class> class Solver>
    class ISTLBackend_AMG_NOVLP : public LinearResultStorage
    {
      typedef typename GO::Traits::TrialGridFunctionSpace GFS;
      typedef typename ISTL::ParallelHelper<GFS> PHELPER;
      typedef typename GO::Traits::Jacobian M;
      typedef Backend::Native<M> MatrixType;
      typedef typename GO::Traits::Domain V;
      typedef Backend::Native<V> VectorType;
      typedef typename ISTL::CommSelector<s,Dune::MPIHelper::isFake>::type Comm;
#if HAVE_MPI
      typedef Preconditioner<MatrixType,VectorType,VectorType,1> Smoother;
      typedef Dune::NonoverlappingBlockPreconditioner<Comm,Smoother> ParSmoother;
      typedef Dune::NonoverlappingSchwarzOperator<MatrixType,VectorType,VectorType,Comm> Operator;
#else
      typedef Preconditioner<MatrixType,VectorType,VectorType,1> ParSmoother;
      typedef Dune::MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
#endif
      typedef typename Dune::Amg::SmootherTraits<ParSmoother>::Arguments SmootherArgs;
      typedef Dune::Amg::AMG<Operator,VectorType,ParSmoother,Comm> AMG;
      typedef Dune::Amg::Parameters Parameters;

    public:
      ISTLBackend_AMG_NOVLP(const GO& grid_operator, unsigned maxiter_=5000,
                            int verbose_=1, bool reuse_=false,
                            bool usesuperlu_=true)
        : _grid_operator(grid_operator)
        , gfs(grid_operator.trialGridFunctionSpace())
        , phelper(gfs,verbose_)
        , maxiter(maxiter_)
        , params(15,2000,1.2,1.6,Dune::Amg::atOnceAccu)
        , verbose(verbose_)
        , reuse(reuse_)
        , firstapply(true)
        , usesuperlu(usesuperlu_)
      {
        params.setDefaultValuesIsotropic(GFS::Traits::GridViewType::Traits::Grid::dimension);
        params.setDebugLevel(verbose_);
#if !HAVE_SUPERLU
        if (phelper.rank() == 0 && usesuperlu == true)
          {
            std::cout << "WARNING: You are using AMG without SuperLU!"
                      << " Please consider installing SuperLU,"
                      << " or set the usesuperlu flag to false"
                      << " to suppress this warning." << std::endl;
          }
#endif
      }

       /*! \brief set AMG parameters

        \param[in] params_ a parameter object of Type Dune::Amg::Parameters
      */
      void setParameters(const Parameters& params_)
      {
        params = params_;
      }

      /**
       * @brief Get the parameters describing the behaviuour of AMG.
       *
       * The returned object can be adjusted to ones needs and then can be
       * reset using setParameters.
       * @return The object holding the parameters of AMG.
       */
      const Parameters& parameters() const
      {
        return params;
      }

      //! Set whether the AMG should be reused again during call to apply().
      void setReuse(bool reuse_)
      {
        reuse = reuse_;
      }

      //! Return whether the AMG is reused during call to apply()
      bool getReuse() const
      {
        return reuse;
      }


      /*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      typename V::ElementType norm (const V& v) const
      {
        V x(v); // make a copy because it has to be made consistent
        typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        psp.make_consistent(x);
        return psp.norm(x);
      }

      void apply(M& A, V& z, V& r, typename Dune::template FieldTraits<typename V::ElementType >::real_type reduction)
      {
        Timer watch;
        MatrixType& mat = Backend::native(A);
        typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<MatrixType,
          Dune::Amg::FirstDiagonal> > Criterion;
#if HAVE_MPI
        Comm oocc(gfs.gridView().comm(),Dune::SolverCategory::nonoverlapping);
        _grid_operator.make_consistent(A);
        phelper.createIndexSetAndProjectForAMG(A, oocc);
        Dune::NonoverlappingSchwarzScalarProduct<VectorType,Comm> sp(oocc);
        Operator oop(mat, oocc);
#else
        Comm oocc(gfs.gridView().comm());
        Operator oop(mat);
        Dune::SeqScalarProduct<VectorType> sp;
#endif
        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1;
        //use noAccu or atOnceAccu
        Criterion criterion(params);
        stats.tprepare=watch.elapsed();
        watch.reset();

        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        //only construct a new AMG if the matrix changes
        if (reuse==false || firstapply==true){
          amg.reset(new AMG(oop, criterion, smootherArgs, oocc));
          firstapply = false;
          stats.tsetup = watch.elapsed();
          stats.levels = amg->maxlevels();
          stats.directCoarseLevelSolver=amg->usesDirectCoarseLevelSolver();
        }

        Dune::InverseOperatorResult stat;
        // make r consistent
        if (gfs.gridView().comm().size()>1) {
          Dune::PDELab::AddDataHandle<GFS,V> adddh(gfs,r);
          gfs.gridView().communicate(adddh,
                                     Dune::InteriorBorder_InteriorBorder_Interface,
                                     Dune::ForwardCommunication);
        }
        watch.reset();
        Solver<VectorType> solver(oop,sp,*amg,reduction,maxiter,verb);
        solver.apply(Backend::native(z),Backend::native(r),stat);
        stats.tsolve= watch.elapsed();
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }

      /**
       * @brief Get statistics of the AMG solver (no of levels, timings).
       * @return statistis of the AMG solver.
       */
      const ISTLAMGStatistics& statistics() const
      {
        return stats;
      }

    private:
      const GO& _grid_operator;
      const GFS& gfs;
      PHELPER phelper;
      unsigned maxiter;
      Parameters params;
      int verbose;
      bool reuse;
      bool firstapply;
      bool usesuperlu;
      std::shared_ptr<AMG> amg;
      ISTLAMGStatistics stats;
    };

    //! \addtogroup PDELab_novlpsolvers Nonoverlapping Solvers
    //! \{

  /**
   * @brief Nonoverlapping parallel CG solver preconditioned with AMG smoothed by SSOR.
   * @tparam GO The type of the grid operator
   * (or the fakeGOTraits class for the old grid operator space).
   * This class will be used to adjust the discretization matrix.
   * and extract the trial grid function space.
   * @tparam s The bits to use for the global index.
   *
   * The solver uses AMG with underlying
   * sequential SSOR preconditioner. The crucial step is to add up the matrix entries
   * corresponding to the border vertices on each process. This is achieved by
   * performing a VertexExchanger::sumEntries(Matrix&).
   */

    template<class GO, int s=96>
    class ISTLBackend_NOVLP_CG_AMG_SSOR
      : public ISTLBackend_AMG_NOVLP<GO, s, Dune::SeqSSOR, Dune::CGSolver>
    {

    public:
      ISTLBackend_NOVLP_CG_AMG_SSOR(const GO& grid_operator, unsigned maxiter_=5000,
                                    int verbose_=1, bool reuse_=false,
                                    bool usesuperlu_=true)
        : ISTLBackend_AMG_NOVLP<GO, s, Dune::SeqSSOR, Dune::CGSolver>(grid_operator, maxiter_,verbose_,reuse_,usesuperlu_)
      {}
    };

  /**
   * @brief Nonoverlapping parallel BiCGStab solver preconditioned with AMG smoothed by SSOR.
   * @tparam GO The type of the grid operator
   * (or the fakeGOTraits class for the old grid operator space).
   * This class will be used to adjust the discretization matrix.
   * and extract the trial grid function space.
   * @tparam s The bits to use for the global index.
   *
   * The solver uses AMG with underlying
   * sequential SSOR preconditioner. The crucial step is to add up the matrix entries
   * corresponding to the border vertices on each process. This is achieved by
   * performing a VertexExchanger::sumEntries(Matrix&).
   */

    template<class GO, int s=96>
    class ISTLBackend_NOVLP_BCGS_AMG_SSOR
      : public ISTLBackend_AMG_NOVLP<GO, s, Dune::SeqSSOR, Dune::BiCGSTABSolver>
    {

    public:
      ISTLBackend_NOVLP_BCGS_AMG_SSOR(const GO& grid_operator, unsigned maxiter_=5000,
                                      int verbose_=1, bool reuse_=false,
                                      bool usesuperlu_=true)
        : ISTLBackend_AMG_NOVLP<GO, s, Dune::SeqSSOR, Dune::BiCGSTABSolver>(grid_operator, maxiter_,verbose_,reuse_,usesuperlu_)
      {}
    };

  /**
   * @brief Nonoverlapping parallel LoopSolver preconditioned with AMG smoothed by SSOR.
   * @tparam GO The type of the grid operator used for the spatial discretization
   * (or the fakeGOTraits class for the old grid operator space).
   * This class will be used to adjust the discretization matrix.
   * and extract the trial grid function space.
   * @tparam s The bits to use for the global index.
   *
   * The solver uses AMG with underlying
   * sequential SSOR preconditioner. The crucial step is to add up the matrix entries
   * corresponding to the border vertices on each process. This is achieved by
   * performing a VertexExchanger::sumEntries(Matrix&).
   */

    template<class GO, int s=96>
    class ISTLBackend_NOVLP_LS_AMG_SSOR
      : public ISTLBackend_AMG_NOVLP<GO, s, Dune::SeqSSOR, Dune::LoopSolver>
    {

    public:
      ISTLBackend_NOVLP_LS_AMG_SSOR(const GO& grid_operator, unsigned maxiter_=5000,
                                    int verbose_=1, bool reuse_=false,
                                    bool usesuperlu_=true)
        : ISTLBackend_AMG_NOVLP<GO, s, Dune::SeqSSOR, Dune::LoopSolver>(grid_operator, maxiter_,verbose_,reuse_,usesuperlu_)
      {}
    };
    //! \} Nonoverlapping Solvers
    //! \} group Backend

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_NOVLPISTLSOLVERBACKEND_HH
