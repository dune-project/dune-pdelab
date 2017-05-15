// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_OVLPISTLSOLVERBACKEND_HH
#define DUNE_PDELAB_BACKEND_ISTL_OVLPISTLSOLVERBACKEND_HH

// this is here for backwards compatibility and deprecation warnings, remove after 2.5.0
#include "ensureistlinclude.hh"

#include <dune/common/deprecated.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/io.hh>
#include <dune/istl/superlu.hh>

#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/backend/istl/vector.hh>
#include <dune/pdelab/backend/istl/bcrsmatrix.hh>
#include <dune/pdelab/backend/istl/parallelhelper.hh>
#include <dune/pdelab/backend/istl/seqistlsolverbackend.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Backend
    //! \ingroup PDELab
    //! \{

    //========================================================
    // Generic support for overlapping grids
    // (need to be used with appropriate constraints)
    //========================================================

    // operator that resets result to zero at constrained DOFS
    template<class CC, class M, class X, class Y>
    class OverlappingOperator
      : public Dune::AssembledLinearOperator<M,X,Y>
    {
    public:
      //! export types
      typedef M matrix_type;
      typedef X domain_type;
      typedef Y range_type;
      typedef typename X::ElementType field_type;

      OverlappingOperator (const CC& cc_, const M& A)
        : cc(cc_), _A_(A)
      {}

      //! apply operator to x:  \f$ y = A(x) \f$
      virtual void apply (const domain_type& x, range_type& y) const
      {
        using Backend::native;
        native(_A_).mv(native(x),native(y));
        Dune::PDELab::set_constrained_dofs(cc,0.0,y);
      }

      //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
      virtual void applyscaleadd (field_type alpha, const domain_type& x, range_type& y) const
      {
        using Backend::native;
        native(_A_).usmv(alpha,native(x),native(y));
        Dune::PDELab::set_constrained_dofs(cc,0.0,y);
      }

      SolverCategory::Category category() const override
      {
        return SolverCategory::overlapping;
      }

      //! get matrix via *
      virtual const M& getmat () const
      {
        return _A_;
      }

    private:
      const CC& cc;
      const M& _A_;
    };

    // new scalar product assuming at least overlap 1
    // uses unique partitioning of nodes for parallelization
    template<class GFS, class X>
    class OverlappingScalarProduct
      : public Dune::ScalarProduct<X>
    {
    public:
      //! export types
      typedef X domain_type;
      typedef typename X::ElementType field_type;

      /*! \brief Constructor needs to know the grid function space
       */
      OverlappingScalarProduct (const GFS& gfs_, const ISTL::ParallelHelper<GFS>& helper_)
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

      SolverCategory::Category category() const override
      {
        return SolverCategory::overlapping;
      }

    private:
      const GFS& gfs;
      const ISTL::ParallelHelper<GFS>& helper;
    };

    // wrapped sequential preconditioner
    template<class CC, class GFS, class P>
    class OverlappingWrappedPreconditioner
      : public Dune::Preconditioner<Dune::PDELab::Backend::Vector<GFS,typename P::domain_type::field_type>,
                                    Dune::PDELab::Backend::Vector<GFS,typename P::range_type::field_type>>
    {
    public:
      //! \brief The domain type of the preconditioner.
      using domain_type = Dune::PDELab::Backend::Vector<GFS,typename P::domain_type::field_type>;
      //! \brief The range type of the preconditioner.
      using range_type = Dune::PDELab::Backend::Vector<GFS,typename P::range_type::field_type>;

      //! Constructor.
      OverlappingWrappedPreconditioner (const GFS& gfs_, P& prec_, const CC& cc_,
                                        const ISTL::ParallelHelper<GFS>& helper_)
        : gfs(gfs_), prec(prec_), cc(cc_), helper(helper_)
      {}

      /*!
        \brief Prepare the preconditioner.
      */
      virtual void pre (domain_type& x, range_type& b)
      {
        prec.pre(x,b);
      }

      /*!
        \brief Apply the preconditioner.
      */
      virtual void apply (domain_type& v, const range_type& d)
      {
        range_type dd(d);
        set_constrained_dofs(cc,0.0,dd);
        prec.apply(Backend::native(v),Backend::native(dd));
        Dune::PDELab::AddDataHandle<GFS,domain_type> adddh(gfs,v);
        if (gfs.gridView().comm().size()>1)
          gfs.gridView().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);
      }

      SolverCategory::Category category() const override
      {
        return SolverCategory::overlapping;
      }

      /*!
        \brief Clean up.
      */
      virtual void post (domain_type& x)
      {
        prec.post(Backend::native(x));
      }

    private:
      const GFS& gfs;
      P& prec;
      const CC& cc;
      const ISTL::ParallelHelper<GFS>& helper;
    };


#if HAVE_SUITESPARSE_UMFPACK || DOXYGEN
    // exact subdomain solves with UMFPack as preconditioner
    template<class GFS, class M, class X, class Y>
    class UMFPackSubdomainSolver : public Dune::Preconditioner<X,Y>
    {
      typedef Backend::Native<M> ISTLM;

    public:
      //! \brief The domain type of the preconditioner.
      typedef X domain_type;
      //! \brief The range type of the preconditioner.
      typedef Y range_type;
      //! \brief The field type of the preconditioner.
      typedef typename X::ElementType field_type;

      /*! \brief Constructor.

        Constructor gets all parameters to operate the prec.
        \param gfs_ The grid function space.
        \param A_ The matrix to operate on.
      */
      UMFPackSubdomainSolver (const GFS& gfs_, const M& A_)
        : gfs(gfs_), solver(Backend::native(A_),false) // this does the decomposition
      {}

      /*!
        \brief Prepare the preconditioner.
      */
      virtual void pre (X& x, Y& b) {}

      /*!
        \brief Apply the precondioner.
      */
      virtual void apply (X& v, const Y& d)
      {
        Dune::InverseOperatorResult stat;
        Y b(d); // need copy, since solver overwrites right hand side
        solver.apply(Backend::native(v),Backend::native(b),stat);
        if (gfs.gridView().comm().size()>1)
          {
            AddDataHandle<GFS,X> adddh(gfs,v);
            gfs.gridView().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);
          }
      }

      SolverCategory::Category category() const override
      {
        return SolverCategory::overlapping;
      }

      /*!
        \brief Clean up.
      */
      virtual void post (X& x) {}

    private:
      const GFS& gfs;
      Dune::UMFPack<ISTLM> solver;
    };
#endif

#if HAVE_SUPERLU
    // exact subdomain solves with SuperLU as preconditioner
    template<class GFS, class M, class X, class Y>
    class SuperLUSubdomainSolver : public Dune::Preconditioner<X,Y>
    {
      typedef Backend::Native<M> ISTLM;

    public:
      //! \brief The domain type of the preconditioner.
      typedef X domain_type;
      //! \brief The range type of the preconditioner.
      typedef Y range_type;
      //! \brief The field type of the preconditioner.
      typedef typename X::ElementType field_type;

      /*! \brief Constructor.

        Constructor gets all parameters to operate the prec.
        \param gfs_ The grid function space.
        \param A_ The matrix to operate on.
      */
      SuperLUSubdomainSolver (const GFS& gfs_, const M& A_)
        : gfs(gfs_), solver(Backend::native(A_),false) // this does the decomposition
      {}

      /*!
        \brief Prepare the preconditioner.
      */
      virtual void pre (X& x, Y& b) {}

      /*!
        \brief Apply the precondioner.
      */
      virtual void apply (X& v, const Y& d)
      {
        Dune::InverseOperatorResult stat;
        Y b(d); // need copy, since solver overwrites right hand side
        solver.apply(Backend::native(v),Backend::native(b),stat);
        if (gfs.gridView().comm().size()>1)
          {
            AddDataHandle<GFS,X> adddh(gfs,v);
            gfs.gridView().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);
          }
      }

      SolverCategory::Category category() const override
      {
        return SolverCategory::overlapping;
      }

      /*!
        \brief Clean up.
      */
      virtual void post (X& x) {}

    private:
      const GFS& gfs;
      Dune::SuperLU<ISTLM> solver;
    };

    // exact subdomain solves with SuperLU as preconditioner
    template<class GFS, class M, class X, class Y>
    class RestrictedSuperLUSubdomainSolver : public Dune::Preconditioner<X,Y>
    {
      typedef typename M::Container ISTLM;

    public:
      //! \brief The domain type of the preconditioner.
      typedef X domain_type;
      //! \brief The range type of the preconditioner.
      typedef Y range_type;
      //! \brief The field type of the preconditioner.
      typedef typename X::ElementType field_type;

      /*! \brief Constructor.

        Constructor gets all parameters to operate the prec.
        \param gfs_ The grid function space.
        \param A_ The matrix to operate on.
        \param helper_ The parallel istl helper.
      */
      RestrictedSuperLUSubdomainSolver (const GFS& gfs_, const M& A_,
                                        const ISTL::ParallelHelper<GFS>& helper_)
        : gfs(gfs_), solver(Backend::native(A_),false), helper(helper_) // this does the decomposition
      {}

      /*!
        \brief Prepare the preconditioner.
      */
      virtual void pre (X& x, Y& b) {}

      /*!
        \brief Apply the precondioner.
      */
      virtual void apply (X& v, const Y& d)
      {
        using Backend::native;
        Dune::InverseOperatorResult stat;
        Y b(d); // need copy, since solver overwrites right hand side
        solver.apply(native(v),native(b),stat);
        if (gfs.gridView().comm().size()>1)
          {
            helper.maskForeignDOFs(native(v));
            AddDataHandle<GFS,X> adddh(gfs,v);
            gfs.gridView().communicate(adddh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
          }
      }

      SolverCategory::Category category() const override
      {
        return SolverCategory::overlapping;
      }

      /*!
        \brief Clean up.
      */
      virtual void post (X& x) {}

    private:
      const GFS& gfs;
      Dune::SuperLU<ISTLM> solver;
      const ISTL::ParallelHelper<GFS>& helper;
    };
#endif

    template<typename GFS>
    class OVLPScalarProductImplementation
    {
    public:
      OVLPScalarProductImplementation(const GFS& gfs_)
        : gfs(gfs_), helper(gfs_)
      {}

      /*! \brief Dot product of two vectors.
        It is assumed that the vectors are consistent on the interior+border
        partition.
      */
      template<typename X>
      typename X::ElementType dot (const X& x, const X& y) const
      {
        // do local scalar product on unique partition
        typename X::ElementType sum = helper.disjointDot(x,y);

        // do global communication
        return gfs.gridView().comm().sum(sum);
      }

      /*! \brief Norm of a right-hand side vector.
        The vector must be consistent on the interior+border partition
      */
      template<typename X>
      typename Dune::template FieldTraits<typename X::ElementType >::real_type norm (const X& x) const
      {
        using namespace std;
        return sqrt(static_cast<double>(this->dot(x,x)));
      }

      const ISTL::ParallelHelper<GFS>& parallelHelper() const
      {
        return helper;
      }

      // need also non-const version;
      ISTL::ParallelHelper<GFS>& parallelHelper() // P.B.: needed for createIndexSetAndProjectForAMG
      {
        return helper;
      }

    private:
      const GFS& gfs;
      ISTL::ParallelHelper<GFS> helper;
    };


    template<typename GFS, typename X>
    class OVLPScalarProduct
      : public ScalarProduct<X>
    {
    public:
      SolverCategory::Category category() const override
      {
        return SolverCategory::overlapping;
      }

      OVLPScalarProduct(const OVLPScalarProductImplementation<GFS>& implementation_)
        : implementation(implementation_)
      {}

      virtual typename X::Container::field_type dot(const X& x, const X& y)
      {
        return implementation.dot(x,y);
      }

      virtual typename X::Container::field_type norm (const X& x)
      {
        using namespace std;
        return sqrt(static_cast<double>(this->dot(x,x)));
      }

    private:
      const OVLPScalarProductImplementation<GFS>& implementation;
    };

    template<class GFS, class C,
             template<class,class,class,int> class Preconditioner,
             template<class> class Solver>
    class ISTLBackend_OVLP_Base
      : public OVLPScalarProductImplementation<GFS>, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] gfs_ a grid function space
        \param[in] c_ a constraints object
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] steps_ number of SSOR steps to apply as inner iteration
        \param[in] verbose_ print messages if true
      */
      ISTLBackend_OVLP_Base (const GFS& gfs_, const C& c_, unsigned maxiter_=5000,
                                            int steps_=5, int verbose_=1)
        : OVLPScalarProductImplementation<GFS>(gfs_), gfs(gfs_), c(c_), maxiter(maxiter_), steps(steps_), verbose(verbose_)
      {}

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename V::ElementType >::real_type reduction)
      {
        using Backend::Native;
        using Backend::native;
        typedef OverlappingOperator<C,M,V,W> POP;
        POP pop(c,A);
        typedef OVLPScalarProduct<GFS,V> PSP;
        PSP psp(*this);
        typedef Preconditioner<
          Native<M>,
          Native<V>,
          Native<W>,
          1
          > SeqPrec;
        SeqPrec seqprec(native(A),steps,1.0);
        typedef OverlappingWrappedPreconditioner<C,GFS,SeqPrec> WPREC;
        WPREC wprec(gfs,seqprec,c,this->parallelHelper());
        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        Solver<V> solver(pop,psp,wprec,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }
    private:
      const GFS& gfs;
      const C& c;
      unsigned maxiter;
      int steps;
      int verbose;
    };

    // Base class for ILU0 as preconditioner
    template<class GFS, class C,
             template<class> class Solver>
    class ISTLBackend_OVLP_ILU0_Base
      : public OVLPScalarProductImplementation<GFS>, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] gfs_ a grid function space
        \param[in] c_ a constraints object
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ verbosity level (0=silent)
      */
      ISTLBackend_OVLP_ILU0_Base (const GFS& gfs_, const C& c_, unsigned maxiter_=5000, int verbose_=1)
        : OVLPScalarProductImplementation<GFS>(gfs_), gfs(gfs_), c(c_), maxiter(maxiter_), verbose(verbose_)
      {}

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename V::ElementType >::real_type reduction)
      {
        using Backend::Native;
        using Backend::native;
        typedef OverlappingOperator<C,M,V,W> POP;
        POP pop(c,A);
        typedef OVLPScalarProduct<GFS,V> PSP;
        PSP psp(*this);
        typedef SeqILU0<
          Native<M>,
          Native<V>,
          Native<W>,
          1
          > SeqPrec;
        SeqPrec seqprec(native(A),1.0);
        typedef OverlappingWrappedPreconditioner<C,GFS,SeqPrec> WPREC;
        WPREC wprec(gfs,seqprec,c,this->parallelHelper());
        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        Solver<V> solver(pop,psp,wprec,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }
    private:
      const GFS& gfs;
      const C& c;
      unsigned maxiter;
      int steps;
      int verbose;
    };

    // Base class for ILUn as preconditioner
    template<class GFS, class C,
             template<class> class Solver>
    class ISTLBackend_OVLP_ILUn_Base
      : public OVLPScalarProductImplementation<GFS>, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] gfs_ a grid function space
        \param[in] c_ a constraints object
        \param[in] n_ level for ILUn
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ verbosity level (0=silent)
      */
      ISTLBackend_OVLP_ILUn_Base (const GFS& gfs_, const C& c_, int n_=1, unsigned maxiter_=5000, int verbose_=1)
        : OVLPScalarProductImplementation<GFS>(gfs_), gfs(gfs_), c(c_), n(n_), maxiter(maxiter_), verbose(verbose_)
      {}

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename V::ElementType >::real_type reduction)
      {
        using Backend::Native;
        using Backend::native;
        typedef OverlappingOperator<C,M,V,W> POP;
        POP pop(c,A);
        typedef OVLPScalarProduct<GFS,V> PSP;
        PSP psp(*this);
        typedef SeqILUn<
          Native<M>,
          Native<V>,
          Native<W>,
          1
          > SeqPrec;
        SeqPrec seqprec(native(A),n,1.0);
        typedef OverlappingWrappedPreconditioner<C,GFS,SeqPrec> WPREC;
        WPREC wprec(gfs,seqprec,c,this->parallelHelper());
        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        Solver<V> solver(pop,psp,wprec,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }
    private:
      const GFS& gfs;
      const C& c;
      int n;
      unsigned maxiter;
      int steps;
      int verbose;
    };

    //! \addtogroup PDELab_ovlpsolvers Overlapping Solvers
    //! \{

    /**
     * @brief Overlapping parallel BiCGStab solver with SSOR preconditioner
     * @tparam GFS The Type of the GridFunctionSpace.
     * @tparam CC The Type of the Constraints Container.
     */
    template<class GFS, class CC>
    class ISTLBackend_OVLP_BCGS_SSORk
      : public ISTLBackend_OVLP_Base<GFS,CC,Dune::SeqSSOR, Dune::BiCGSTABSolver>
    {
    public:
      /*! \brief make a linear solver object

        \param[in] gfs a grid function space
        \param[in] cc a constraints container object
        \param[in] maxiter maximum number of iterations to do
        \param[in] steps number of SSOR steps to apply as inner iteration
        \param[in] verbose print messages if true
      */
      ISTLBackend_OVLP_BCGS_SSORk (const GFS& gfs, const CC& cc, unsigned maxiter=5000,
                                            int steps=5, int verbose=1)
        : ISTLBackend_OVLP_Base<GFS,CC,Dune::SeqSSOR, Dune::BiCGSTABSolver>(gfs, cc, maxiter, steps, verbose)
      {}
    };
    /**
     * @brief Overlapping parallel BiCGStab solver with ILU0 preconditioner
     * @tparam GFS The Type of the GridFunctionSpace.
     * @tparam CC The Type of the Constraints Container.
     */
    template<class GFS, class CC>
    class ISTLBackend_OVLP_BCGS_ILU0
      : public ISTLBackend_OVLP_ILU0_Base<GFS,CC,Dune::BiCGSTABSolver>
    {
    public:
      /*! \brief make a linear solver object

        \param[in] gfs a grid function space
        \param[in] cc a constraints container object
        \param[in] maxiter maximum number of iterations to do
        \param[in] verbose print messages if true
      */
      ISTLBackend_OVLP_BCGS_ILU0 (const GFS& gfs, const CC& cc, unsigned maxiter=5000, int verbose=1)
        : ISTLBackend_OVLP_ILU0_Base<GFS,CC,Dune::BiCGSTABSolver>(gfs, cc, maxiter, verbose)
      {}
    };
    /**
     * @brief Overlapping parallel BiCGStab solver with ILU0 preconditioner
     * @tparam GFS The Type of the GridFunctionSpace.
     * @tparam CC The Type of the Constraints Container.
     */
    template<class GFS, class CC>
    class ISTLBackend_OVLP_BCGS_ILUn
      : public ISTLBackend_OVLP_ILUn_Base<GFS,CC,Dune::BiCGSTABSolver>
    {
    public:
      /*! \brief make a linear solver object

        \param[in] gfs a grid function space
        \param[in] cc a constraints container object
        \param[in] n level for ILUn
        \param[in] maxiter maximum number of iterations to do
        \param[in] verbose print messages if true
      */
      ISTLBackend_OVLP_BCGS_ILUn (const GFS& gfs, const CC& cc, int n=1, unsigned maxiter=5000, int verbose=1)
        : ISTLBackend_OVLP_ILUn_Base<GFS,CC,Dune::BiCGSTABSolver>(gfs, cc, n, maxiter, verbose)
      {}
    };
    /**
     * @brief Overlapping parallel CGS solver with SSOR preconditioner
     * @tparam GFS The Type of the GridFunctionSpace.
     * @tparam CC The Type of the Constraints Container.
     */
    template<class GFS, class CC>
    class ISTLBackend_OVLP_CG_SSORk
      : public ISTLBackend_OVLP_Base<GFS,CC,Dune::SeqSSOR, Dune::CGSolver>
    {
    public:
      /*! \brief make a linear solver object

        \param[in] gfs a grid function space
        \param[in] cc a constraints container object
        \param[in] maxiter maximum number of iterations to do
        \param[in] steps number of SSOR steps to apply as inner iteration
        \param[in] verbose print messages if true
      */
      ISTLBackend_OVLP_CG_SSORk (const GFS& gfs, const CC& cc, unsigned maxiter=5000,
                                            int steps=5, int verbose=1)
        : ISTLBackend_OVLP_Base<GFS,CC,Dune::SeqSSOR, Dune::CGSolver>(gfs, cc, maxiter, steps, verbose)
      {}
    };

    /**
     * @brief Overlapping parallel restarted GMRes solver with ILU0 preconditioner
     * @tparam GFS The Type of the GridFunctionSpace.
     * @tparam CC The Type of the Constraints Container.
     */
    template<class GFS, class CC>
    class ISTLBackend_OVLP_GMRES_ILU0
      : public OVLPScalarProductImplementation<GFS>, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] gfs a grid function space
        \param[in] cc a constraints container object
        \param[in] maxiter maximum number of iterations to do
        \param[in] verbose print messages if true
      */
        ISTLBackend_OVLP_GMRES_ILU0 (const GFS& gfs_, const CC& cc_, unsigned maxiter_=5000, int verbose_=1,
                                     int restart_ = 20)
        : OVLPScalarProductImplementation<GFS>(gfs_), gfs(gfs_), cc(cc_), maxiter(maxiter_), verbose(verbose_),
          restart(restart_)
      {}

      /*! \brief solve the given linear system
        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename V::ElementType >::real_type reduction)
      {
        using Backend::Native;
        using Backend::native;
        typedef OverlappingOperator<CC,M,V,W> POP;
        POP pop(cc,A);
        typedef OVLPScalarProduct<GFS,V> PSP;
        PSP psp(*this);
        typedef SeqILU0<
          Native<M>,
          Native<V>,
          Native<W>,
          1
          > SeqPrec;
        SeqPrec seqprec(native(A),1.0);
        typedef OverlappingWrappedPreconditioner<CC,GFS,SeqPrec> WPREC;
        WPREC wprec(gfs,seqprec,cc,this->parallelHelper());
        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        RestartedGMResSolver<V> solver(pop,psp,wprec,reduction,restart,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }

    private:
      const GFS& gfs;
      const CC& cc;
      unsigned maxiter;
      int steps;
      int verbose;
      int restart;
    };

    //! \} Solver

    template<class GFS, class C, template<typename> class Solver>
    class ISTLBackend_OVLP_SuperLU_Base
      : public OVLPScalarProductImplementation<GFS>, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] gfs_ a grid function space
        \param[in] c_ a constraints object
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      ISTLBackend_OVLP_SuperLU_Base (const GFS& gfs_, const C& c_, unsigned maxiter_=5000,
                                              int verbose_=1)
        : OVLPScalarProductImplementation<GFS>(gfs_), gfs(gfs_), c(c_), maxiter(maxiter_), verbose(verbose_)
      {}

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename V::ElementType >::real_type reduction)
      {
        typedef OverlappingOperator<C,M,V,W> POP;
        POP pop(c,A);
        typedef OVLPScalarProduct<GFS,V> PSP;
        PSP psp(*this);
#if HAVE_SUPERLU
        typedef SuperLUSubdomainSolver<GFS,M,V,W> PREC;
        PREC prec(gfs,A);
        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        Solver<V> solver(pop,psp,prec,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
#else
        std::cout << "No superLU support, please install and configure it." << std::endl;
#endif
      }

    private:
      const GFS& gfs;
      const C& c;
      unsigned maxiter;
      int verbose;
    };

        //! \} Solver

    template<class GFS, class C, template<typename> class Solver>
    class ISTLBackend_OVLP_UMFPack_Base
      : public OVLPScalarProductImplementation<GFS>, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] gfs_ a grid function space
        \param[in] c_ a constraints object
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      ISTLBackend_OVLP_UMFPack_Base (const GFS& gfs_, const C& c_, unsigned maxiter_=5000,
                                              int verbose_=1)
        : OVLPScalarProductImplementation<GFS>(gfs_), gfs(gfs_), c(c_), maxiter(maxiter_), verbose(verbose_)
      {}

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename V::ElementType >::real_type reduction)
      {
        typedef OverlappingOperator<C,M,V,W> POP;
        POP pop(c,A);
        typedef OVLPScalarProduct<GFS,V> PSP;
        PSP psp(*this);
#if HAVE_SUITESPARSE_UMFPACK || DOXYGEN
        typedef UMFPackSubdomainSolver<GFS,M,V,W> PREC;
        PREC prec(gfs,A);
        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        Solver<V> solver(pop,psp,prec,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
#else
        std::cout << "No UMFPack support, please install and configure it." << std::endl;
#endif
      }

    private:
      const GFS& gfs;
      const C& c;
      unsigned maxiter;
      int verbose;
    };

    //! \addtogroup PDELab_ovlpsolvers Overlapping Solvers
    //! \{
    /**
     * @brief Overlapping parallel BiCGStab solver with SuperLU preconditioner
     * @tparam GFS The Type of the GridFunctionSpace.
     * @tparam CC The Type of the Constraints Container.
     */
    template<class GFS, class CC>
    class ISTLBackend_OVLP_BCGS_SuperLU
      : public ISTLBackend_OVLP_SuperLU_Base<GFS,CC,Dune::BiCGSTABSolver>
    {
    public:

      /*! \brief make a linear solver object

        \param[in] gfs_ a grid function space
        \param[in] cc_ a constraints container object
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      ISTLBackend_OVLP_BCGS_SuperLU (const GFS& gfs_, const CC& cc_, unsigned maxiter_=5000,
                                              int verbose_=1)
        : ISTLBackend_OVLP_SuperLU_Base<GFS,CC,Dune::BiCGSTABSolver>(gfs_,cc_,maxiter_,verbose_)
      {}
    };

    /**
     * @brief Overlapping parallel CG solver with SuperLU preconditioner
     * @tparam GFS The Type of the GridFunctionSpace.
     * @tparam CC The Type of the Constraints Container.
     */
    template<class GFS, class CC>
    class ISTLBackend_OVLP_CG_SuperLU
      : public ISTLBackend_OVLP_SuperLU_Base<GFS,CC,Dune::CGSolver>
    {
    public:

      /*! \brief make a linear solver object

        \param[in] gfs_ a grid function space
        \param[in] cc_ a constraints object
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      ISTLBackend_OVLP_CG_SuperLU (const GFS& gfs_, const CC& cc_,
                                              unsigned maxiter_=5000,
                                              int verbose_=1)
        : ISTLBackend_OVLP_SuperLU_Base<GFS,CC,Dune::CGSolver>(gfs_,cc_,maxiter_,verbose_)
      {}
    };

    /**
     * @brief Overlapping parallel CG solver with UMFPack preconditioner
     * @tparam GFS The Type of the GridFunctionSpace.
     * @tparam CC The Type of the Constraints Container.
     */
    template<class GFS, class CC>
    class ISTLBackend_OVLP_CG_UMFPack
      : public ISTLBackend_OVLP_UMFPack_Base<GFS,CC,Dune::CGSolver>
    {
    public:

      /*! \brief make a linear solver object

        \param[in] gfs_ a grid function space
        \param[in] cc_ a constraints object
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      ISTLBackend_OVLP_CG_UMFPack (const GFS& gfs_, const CC& cc_,
                                              unsigned maxiter_=5000,
                                              int verbose_=1)
        : ISTLBackend_OVLP_UMFPack_Base<GFS,CC,Dune::CGSolver>(gfs_,cc_,maxiter_,verbose_)
      {}
    };


    /** @brief Solver to be used for explicit time-steppers with (block-)diagonal mass matrix
     * @tparam GFS The Type of the GridFunctionSpace.
     */
    template<class GFS>
    class ISTLBackend_OVLP_ExplicitDiagonal
      : public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] gfs_ a grid function space
      */
      explicit ISTLBackend_OVLP_ExplicitDiagonal (const GFS& gfs_)
        : gfs(gfs_)
      {}

      explicit ISTLBackend_OVLP_ExplicitDiagonal (const ISTLBackend_OVLP_ExplicitDiagonal& other_)
        : gfs(other_.gfs)
      {}

      /*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      template<class V>
      typename V::ElementType norm(const V& v) const
      {
        static_assert
          (AlwaysFalse<V>::value,
           "ISTLBackend_OVLP_ExplicitDiagonal::norm() should not be "
           "neccessary, so we skipped the implementation.  If you have a "
           "scenario where you need it, please implement it or report back to "
           "us.");
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
        using Backend::Native;
        using Backend::native;
        Dune::SeqJac<
          Native<M>,
          Native<V>,
          Native<W>
          > jac(native(A),1,1.0);
        jac.pre(native(z),native(r));
        jac.apply(native(z),native(r));
        jac.post(native(z));
        if (gfs.gridView().comm().size()>1)
        {
          CopyDataHandle<GFS,V> copydh(gfs,z);
          gfs.gridView().communicate(copydh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
        }
        res.converged  = true;
        res.iterations = 1;
        res.elapsed    = 0.0;
        res.reduction  = static_cast<double>(reduction);
        res.conv_rate  = static_cast<double>(reduction); // pow(reduction,1.0/1)
      }

    private:
      const GFS& gfs;
    };
    //! \} Overlapping Solvers

    template<class GO, int s, template<class,class,class,int> class Preconditioner,
             template<class> class Solver>
    class ISTLBackend_AMG : public LinearResultStorage
    {
      typedef typename GO::Traits::TrialGridFunctionSpace GFS;
      typedef ISTL::ParallelHelper<GFS> PHELPER;
      typedef typename GO::Traits::Jacobian M;
      typedef Backend::Native<M> MatrixType;
      typedef typename GO::Traits::Domain V;
      typedef Backend::Native<V> VectorType;
      typedef typename ISTL::CommSelector<s,Dune::MPIHelper::isFake>::type Comm;
#if HAVE_MPI
      typedef Preconditioner<MatrixType,VectorType,VectorType,1> Smoother;
      typedef Dune::BlockPreconditioner<VectorType,VectorType,Comm,Smoother> ParSmoother;
      typedef Dune::OverlappingSchwarzOperator<MatrixType,VectorType,VectorType,Comm> Operator;
#else
      typedef Preconditioner<MatrixType,VectorType,VectorType,1> ParSmoother;
      typedef Dune::MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
#endif
      typedef typename Dune::Amg::SmootherTraits<ParSmoother>::Arguments SmootherArgs;
      typedef Dune::Amg::AMG<Operator,VectorType,ParSmoother,Comm> AMG;

      typedef typename V::ElementType RF;

    public:

      /**
       * @brief Parameters object to customize matrix hierachy building.
       */
      typedef Dune::Amg::Parameters Parameters;

    public:
      ISTLBackend_AMG(const GFS& gfs_, unsigned maxiter_=5000,
                      int verbose_=1, bool reuse_=false,
                      bool usesuperlu_=true)
        : gfs(gfs_), phelper(gfs,verbose_), maxiter(maxiter_), params(15,2000),
          verbose(verbose_), reuse(reuse_), firstapply(true),
          usesuperlu(usesuperlu_)
      {
        params.setDefaultValuesIsotropic(GFS::Traits::GridViewType::Traits::Grid::dimension);
        params.setDebugLevel(verbose_);
#if !HAVE_SUPERLU
        if (gfs.gridView().comm().rank() == 0 && usesuperlu == true)
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
        typedef OverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        return psp.norm(v);
      }

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      void apply(M& A, V& z, V& r, typename Dune::template FieldTraits<typename V::ElementType >::real_type reduction)
      {
        Timer watch;
        Comm oocc(gfs.gridView().comm());
        MatrixType& mat=Backend::native(A);
        typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<MatrixType,
          Dune::Amg::FirstDiagonal> > Criterion;
#if HAVE_MPI
        phelper.createIndexSetAndProjectForAMG(A, oocc);
        Operator oop(mat, oocc);
        Dune::OverlappingSchwarzScalarProduct<VectorType,Comm> sp(oocc);
#else
        Operator oop(mat);
        Dune::SeqScalarProduct<VectorType> sp;
#endif
        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1;
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
        watch.reset();
        Solver<VectorType> solver(oop,sp,*amg,RF(reduction),maxiter,verb);
        Dune::InverseOperatorResult stat;

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
      const GFS& gfs;
      PHELPER phelper;
      unsigned maxiter;
      Parameters params;
      int verbose;
      bool reuse;
      bool firstapply;
      bool usesuperlu;
      shared_ptr<AMG> amg;
      ISTLAMGStatistics stats;
    };

    //! \addtogroup PDELab_ovlpsolvers Overlapping Solvers
    //! \{

    /**
     * @brief Overlapping parallel conjugate gradient solver preconditioned with AMG smoothed by SSOR
     * @tparam GO The type of the grid operator
     * (or the fakeGOTraits class for the old grid operator space).
     * @tparam s The bits to use for the global index.
     */
    template<class GO, int s=96>
    class ISTLBackend_CG_AMG_SSOR
      : public ISTLBackend_AMG<GO, s, Dune::SeqSSOR, Dune::CGSolver>
    {
      typedef typename GO::Traits::TrialGridFunctionSpace GFS;
    public:
      /**
       * @brief Constructor
       * @param gfs_ The grid function space used.
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       * @param reuse_ Set true, if the Matrix to be used is always identical
       * (AMG aggregation is then only performed once).
       * @param usesuperlu_ Set false, to suppress the no SuperLU warning
       */
      ISTLBackend_CG_AMG_SSOR(const GFS& gfs_, unsigned maxiter_=5000,
                              int verbose_=1, bool reuse_=false,
                              bool usesuperlu_=true)
        : ISTLBackend_AMG<GO, s, Dune::SeqSSOR, Dune::CGSolver>
          (gfs_, maxiter_, verbose_, reuse_, usesuperlu_)
      {}
    };

    /**
     * @brief Overlapping parallel BiCGStab solver preconditioned with AMG smoothed by SSOR.
     * @tparam GO The type of the grid operator
     * (or the fakeGOTraits class for the old grid operator space).
     * @tparam s The bits to use for the globale index.
     */
    template<class GO, int s=96>
    class ISTLBackend_BCGS_AMG_SSOR
      : public ISTLBackend_AMG<GO, s, Dune::SeqSSOR, Dune::BiCGSTABSolver>
    {
      typedef typename GO::Traits::TrialGridFunctionSpace GFS;
    public:
      /**
       * @brief Constructor
       * @param gfs_ The grid function space used.
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       * @param reuse_ Set true, if the Matrix to be used is always identical
       * (AMG aggregation is then only performed once).
       * @param usesuperlu_ Set false, to suppress the no SuperLU warning
       */
      ISTLBackend_BCGS_AMG_SSOR(const GFS& gfs_, unsigned maxiter_=5000,
                                int verbose_=1, bool reuse_=false,
                                bool usesuperlu_=true)
        : ISTLBackend_AMG<GO, s, Dune::SeqSSOR, Dune::BiCGSTABSolver>
          (gfs_, maxiter_, verbose_, reuse_, usesuperlu_)
      {}
    };

    /**
     * @brief Overlapping parallel BiCGStab solver preconditioned with AMG smoothed by ILU0.
     * @tparam GO The type of the grid operator
     * (or the fakeGOTraits class for the old grid operator space).
     * @tparam s The bits to use for the globale index.
     */
    template<class GO, int s=96>
    class ISTLBackend_BCGS_AMG_ILU0
      : public ISTLBackend_AMG<GO, s, Dune::SeqILU0, Dune::BiCGSTABSolver>
    {
      typedef typename GO::Traits::TrialGridFunctionSpace GFS;
    public:
      /**
       * @brief Constructor
       * @param gfs_ The grid function space used.
       * @param maxiter_ The maximum number of iterations allowed.
       * @param verbose_ The verbosity level to use.
       * @param reuse_ Set true, if the Matrix to be used is always identical
       * (AMG aggregation is then only performed once).
       * @param usesuperlu_ Set false, to suppress the no SuperLU warning
       */
      ISTLBackend_BCGS_AMG_ILU0(const GFS& gfs_, unsigned maxiter_=5000,
                                int verbose_=1, bool reuse_=false,
                                bool usesuperlu_=true)
        : ISTLBackend_AMG<GO, s, Dune::SeqILU0, Dune::BiCGSTABSolver>
          (gfs_, maxiter_, verbose_, reuse_, usesuperlu_)
      {}
    };

    //! \} Overlapping Solvers

    //! \} group Backend

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_OVLPISTLSOLVERBACKEND_HH
