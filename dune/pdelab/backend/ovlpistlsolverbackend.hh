// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_OVLPISTLSOLVERBACKEND_HH
#define DUNE_OVLPISTLSOLVERBACKEND_HH

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
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/backend/istl/parallelhelper.hh>
#include <dune/pdelab/backend/seqistlsolverbackend.hh>

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

      //redefine the category, that is the only difference
      enum {category=Dune::SolverCategory::overlapping};

      OverlappingOperator (const CC& cc_, const M& A)
        : cc(cc_), _A_(A)
      {}

      //! apply operator to x:  \f$ y = A(x) \f$
      virtual void apply (const domain_type& x, range_type& y) const
      {
        istl::raw(_A_).mv(istl::raw(x),istl::raw(y));
        Dune::PDELab::set_constrained_dofs(cc,0.0,y);
      }

      //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
      virtual void applyscaleadd (field_type alpha, const domain_type& x, range_type& y) const
      {
        istl::raw(_A_).usmv(alpha,istl::raw(x),istl::raw(y));
        Dune::PDELab::set_constrained_dofs(cc,0.0,y);
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

      //! define the category
      enum {category=Dune::SolverCategory::overlapping};

      /*! \brief Constructor needs to know the grid function space
       */
      OverlappingScalarProduct (const GFS& gfs_, const istl::ParallelHelper<GFS>& helper_)
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

    private:
      const GFS& gfs;
      const istl::ParallelHelper<GFS>& helper;
    };

    // wrapped sequential preconditioner
    template<class CC, class GFS, class P>
    class OverlappingWrappedPreconditioner
      : public Dune::Preconditioner<typename Dune::PDELab::BackendVectorSelector<GFS,typename P::domain_type::field_type>::Type,
                                    typename Dune::PDELab::BackendVectorSelector<GFS,typename P::range_type::field_type>::Type>
    {
    public:
      //! \brief The domain type of the preconditioner.
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,typename P::domain_type::field_type>::Type
      domain_type;
      //! \brief The range type of the preconditioner.
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,typename P::range_type::field_type>::Type
      range_type;

      // define the category
      enum {
        //! \brief The category the preconditioner is part of.
        category=Dune::SolverCategory::overlapping
      };

      //! Constructor.
      OverlappingWrappedPreconditioner (const GFS& gfs_, P& prec_, const CC& cc_,
                                        const istl::ParallelHelper<GFS>& helper_)
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
        prec.apply(istl::raw(v),istl::raw(dd));
        Dune::PDELab::AddDataHandle<GFS,domain_type> adddh(gfs,v);
        if (gfs.gridView().comm().size()>1)
          gfs.gridView().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);
      }

      /*!
        \brief Clean up.
      */
      virtual void post (domain_type& x)
      {
        prec.post(istl::raw(x));
      }

    private:
      const GFS& gfs;
      P& prec;
      const CC& cc;
      const istl::ParallelHelper<GFS>& helper;
    };


#if HAVE_SUPERLU
    // exact subdomain solves with SuperLU as preconditioner
    template<class GFS, class M, class X, class Y>
    class SuperLUSubdomainSolver : public Dune::Preconditioner<X,Y>
    {
      typedef typename M::BaseT ISTLM;

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
        category=Dune::SolverCategory::overlapping
      };

      /*! \brief Constructor.

        Constructor gets all parameters to operate the prec.
        \param gfs_ The grid function space.
        \param A_ The matrix to operate on.
      */
      SuperLUSubdomainSolver (const GFS& gfs_, const M& A_)
        : gfs(gfs_), solver(istl::raw(A_),false) // this does the decomposition
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
        solver.apply(istl::raw(v),istl::raw(b),stat);
        if (gfs.gridView().comm().size()>1)
          {
            AddDataHandle<GFS,X> adddh(gfs,v);
            gfs.gridView().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);
          }
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
      typedef typename M::BaseT ISTLM;

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
        category=Dune::SolverCategory::overlapping
      };

      /*! \brief Constructor.

        Constructor gets all parameters to operate the prec.
        \param gfs_ The grid function space.
        \param A_ The matrix to operate on.
        \param helper_ The parallel istl helper.
      */
      RestrictedSuperLUSubdomainSolver (const GFS& gfs_, const M& A_,
                                        const istl::ParallelHelper<GFS>& helper_)
        : gfs(gfs_), solver(istl::raw(A_),false), helper(helper_) // this does the decomposition
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
        solver.apply(istl::raw(v),istl::raw(b),stat);
        if (gfs.gridView().comm().size()>1)
          {
            helper.maskForeignDOFs(istl::raw(v));
            AddDataHandle<GFS,X> adddh(gfs,v);
            gfs.gridView().communicate(adddh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
          }
      }

      /*!
        \brief Clean up.
      */
      virtual void post (X& x) {}

    private:
      const GFS& gfs;
      Dune::SuperLU<ISTLM> solver;
      const istl::ParallelHelper<GFS>& helper;
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
      typename X::ElementType norm (const X& x) const
      {
        return sqrt(static_cast<double>(this->dot(x,x)));
      }

      const istl::ParallelHelper<GFS>& parallelHelper() const
      {
        return helper;
      }

      // need also non-const version;
      istl::ParallelHelper<GFS>& parallelHelper() // P.B.: needed for createIndexSetAndProjectForAMG
      {
        return helper;
      }

    private:
      const GFS& gfs;
      istl::ParallelHelper<GFS> helper;
    };


    template<typename GFS, typename X>
    class OVLPScalarProduct
      : public ScalarProduct<X>
    {
    public:
      enum {category=Dune::SolverCategory::overlapping};
      OVLPScalarProduct(const OVLPScalarProductImplementation<GFS>& implementation_)
        : implementation(implementation_)
      {}

      virtual typename X::BaseT::field_type dot(const X& x, const X& y)
      {
        return implementation.dot(x,y);
      }

      virtual typename X::BaseT::field_type norm (const X& x)
      {
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
      void apply(M& A, V& z, W& r, typename V::ElementType reduction)
      {
        typedef OverlappingOperator<C,M,V,W> POP;
        POP pop(c,A);
        typedef OVLPScalarProduct<GFS,V> PSP;
        PSP psp(*this);
        typedef Preconditioner<typename M::BaseT,typename V::BaseT,typename W::BaseT,1> SeqPrec;
        SeqPrec seqprec(istl::raw(A),steps,1.0);
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
      void apply(M& A, V& z, W& r, typename V::ElementType reduction)
      {
        typedef OverlappingOperator<C,M,V,W> POP;
        POP pop(c,A);
        typedef OVLPScalarProduct<GFS,V> PSP;
        PSP psp(*this);
        typedef SeqILU0<typename M::BaseT,typename V::BaseT,typename W::BaseT,1> SeqPrec;
        SeqPrec seqprec(istl::raw(A),1.0);
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
      void apply(M& A, V& z, W& r, typename V::ElementType reduction)
      {
        typedef OverlappingOperator<C,M,V,W> POP;
        POP pop(c,A);
        typedef OVLPScalarProduct<GFS,V> PSP;
        PSP psp(*this);
        typedef SeqILUn<typename M::BaseT,typename V::BaseT,typename W::BaseT,1> SeqPrec;
        SeqPrec seqprec(istl::raw(A),n,1.0);
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
            int restart_ = 20, bool recalc_defect_ = false)
        : OVLPScalarProductImplementation<GFS>(gfs_), gfs(gfs_), cc(cc_), maxiter(maxiter_), verbose(verbose_),
          restart(restart_), recalc_defect(recalc_defect_)
      {}

      /*! \brief solve the given linear system
        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename V::ElementType reduction)
      {
        typedef OverlappingOperator<CC,M,V,W> POP;
        POP pop(cc,A);
        typedef OVLPScalarProduct<GFS,V> PSP;
        PSP psp(*this);
        typedef SeqILU0<typename M::BaseT,typename V::BaseT,typename W::BaseT,1> SeqPrec;
        SeqPrec seqprec(istl::raw(A),1.0);
        typedef OverlappingWrappedPreconditioner<CC,GFS,SeqPrec> WPREC;
        WPREC wprec(gfs,seqprec,cc,this->parallelHelper());
        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        RestartedGMResSolver<V> solver(pop,psp,wprec,reduction,restart,maxiter,verb,recalc_defect);
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
      bool recalc_defect;
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
      void apply(M& A, V& z, W& r, typename V::ElementType reduction)
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
      void apply(M& A, V& z, W& r, typename W::ElementType reduction)
      {
        Dune::SeqJac<typename M::BaseT,typename V::BaseT,typename W::BaseT> jac(istl::raw(A),1,1.0);
        jac.pre(istl::raw(z),istl::raw(r));
        jac.apply(istl::raw(z),istl::raw(r));
        jac.post(istl::raw(z));
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
      typedef istl::ParallelHelper<GFS> PHELPER;
      typedef typename GO::Traits::Jacobian M;
      typedef typename M::BaseT MatrixType;
      typedef typename GO::Traits::Domain V;
      typedef typename V::BaseT VectorType;
      typedef typename istl::CommSelector<s,Dune::MPIHelper::isFake>::type Comm;
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

      void setparams(Parameters params_) DUNE_DEPRECATED_MSG("setparams() is deprecated, use setParameters() instead")
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
      void apply(M& A, V& z, V& r, typename V::ElementType reduction)
      {
        Timer watch;
        Comm oocc(gfs.gridView().comm());
        MatrixType& mat=istl::raw(A);
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

        solver.apply(istl::raw(z),istl::raw(r),stat);
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

#endif
