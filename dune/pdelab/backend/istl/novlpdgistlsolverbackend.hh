#ifndef DUNE_NOVLDGISTLSOLVERBACKEND_HH
#define DUNE_NOVLDGISTLSOLVERBACKEND_HH

#include "novlpistlsolverbackend.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup Backend
    //! \ingroup PDELab
    //! \{

    template<typename GO, typename X, typename Y>
    class NonOverlappingDGOnTheFlyOperator : public Dune::LinearOperator<X,Y>
    {
    public:
      typedef X domain_type;
      typedef Y range_type;
      typedef typename X::field_type field_type;

      enum {category=Dune::SolverCategory::nonoverlapping};

      NonOverlappingDGOnTheFlyOperator (const GO& go_)
        : go(go_)
      {}

      virtual void apply (const X& x, Y& y) const
      {
        y = 0.0;
        go.jacobian_apply(x,y);
      }

      virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
      {
        Y temp(y);
        temp = 0.0;
        go.jacobian_apply(x,temp);
        y.axpy(alpha,temp);
      }

    private:
      const GO& go;
    };


    // wrapped sequential preconditioner
    template<class GFS, class P>
    class DGNonOverlappingWrappedPreconditioner
      : public Dune::Preconditioner<Dune::PDELab::Backend::Vector<GFS,typename P::domain_type::field_type>,
                                    Dune::PDELab::Backend::Vector<GFS,typename P::range_type::field_type>>
    {
    public:
      //! \brief The domain type of the preconditioner.
      using domain_type = Dune::PDELab::Backend::Vector<GFS,typename P::domain_type::field_type>;
      //! \brief The range type of the preconditioner.
      using range_type = Dune::PDELab::Backend::Vector<GFS,typename P::range_type::field_type>;

      // define the category
      enum {
        //! \brief The category the preconditioner is part of.
        category=Dune::SolverCategory::nonoverlapping
      };

      //! Constructor.
      DGNonOverlappingWrappedPreconditioner (const GFS& gfs_, P& prec_)
        : gfs(gfs_), prec(prec_)
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
        prec.apply(Backend::native(v),Backend::native(dd));
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
    };



    template<class GO,
             template<class,class,class,int> class Preconditioner,
             template<class> class Solver>
    class ISTLBackend_DGNOVLP_Base
      : public LinearResultStorage
    {
      using GFS = typename GO::Traits::TrialGridFunctionSpace;
      using PHELPER = istl::ParallelHelper<GFS>;

    public:
      /*! \brief make a linear solver object

        \param[in] go_ a grid operator
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] steps_ number of SSOR steps to apply as inner iteration
        \param[in] verbose_ print messages if true
      */
      ISTLBackend_DGNOVLP_Base (const GO& go_, unsigned maxiter_=5000,
                                int steps_=5, int verbose_=1)
        : go(go_), gfs(go_.trialGridFunctionSpace()), phelper(gfs,verbose_), maxiter(maxiter_), steps(steps_), verbose(verbose_)
    {}

      /*! \brief Compute global norm of a vector.

        \param[in] v the given vector
      */
      template<class Vector>
      typename Vector::ElementType norm (const Vector& v) const
      {
        typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,Vector> PSP;
        PSP psp(gfs,phelper);
        return psp.norm(v);
      }

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class V, class W>
      void apply(V& z, W& r, typename Dune::template FieldTraits<typename V::ElementType >::real_type reduction)
      {
        using Backend::Native;
        using Backend::native;

        using M = typename GO::Traits::Jacobian;
        M mat(go);
        go.jacobian(z,mat);

        using POP = Dune::PDELab::NonOverlappingDGOnTheFlyOperator<GO,V,W>;
        POP pop(go);
        using PSP = Dune::PDELab::NonoverlappingScalarProduct<GFS,V>;
        PSP psp(gfs,phelper);
        typedef Preconditioner<
          Native<M>,
          Native<V>,
          Native<W>,
          1
          > SeqPrec;
        SeqPrec seqprec(native(mat),steps,1.0);
        typedef DGNonOverlappingWrappedPreconditioner<GFS,SeqPrec> WPREC;
        WPREC wprec(gfs,seqprec);
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
      const GO& go;
      const GFS& gfs;
      PHELPER phelper;
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
    template<class GO>
    class ISTLBackend_DGNOVLP_CG_SSORk
      : public ISTLBackend_DGNOVLP_Base<GO, Dune::SeqSSOR, Dune::CGSolver>
    {
    public:
      /*! \brief make a linear solver object

        \param[in] gfs a grid function space
        \param[in] cc a constraints container object
        \param[in] maxiter maximum number of iterations to do
        \param[in] steps number of SSOR steps to apply as inner iteration
        \param[in] verbose print messages if true
      */
      ISTLBackend_DGNOVLP_CG_SSORk (const GO& go, unsigned maxiter=5000,
                                      int steps=5, int verbose=1)
        : ISTLBackend_DGNOVLP_Base<GO,Dune::SeqSSOR, Dune::CGSolver>(go, maxiter, steps, verbose)
      {}
    };



    template<class GO>
    class ISTLBackend_DGNOVLP_CG_NOPREC
      : public LinearResultStorage
    {
      using GFS = typename GO::Traits::TrialGridFunctionSpace;
      using PHELPER = istl::ParallelHelper<GFS>;

    public:
      /*! \brief make a linear solver object

        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_DGNOVLP_CG_NOPREC (const GO& go_,
                                              unsigned maxiter_=5000,
                                              int verbose_=1)
        : go(go_), gfs(go_.trialGridFunctionSpace()), phelper(gfs,verbose_), maxiter(maxiter_), verbose(verbose_)
      {}

      /*! \brief Compute global norm of a vector.

        \param[in] v the given vector
      */
      template<class Vector>
      typename Vector::ElementType norm (const Vector& v) const
      {
        typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,Vector> PSP;
        PSP psp(gfs,phelper);
        return psp.norm(v);
      }

      /*! \brief solve the given linear system

        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class V, class W>
      void apply(V& z, W& r, typename Dune::template FieldTraits<typename W::ElementType >::real_type reduction)
      {
        using Backend::Native;
        using Backend::native;

        using POP = Dune::PDELab::NonOverlappingDGOnTheFlyOperator<GO,V,W>;
        POP pop(go);
        using PSP = Dune::PDELab::NonoverlappingScalarProduct<GFS,V>;
        PSP psp(gfs,phelper);
        using PRICH = Dune::PDELab::NonoverlappingRichardson<GFS,V,W>;
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

    private:
      const GO& go;
      const GFS& gfs;
      PHELPER phelper;
      unsigned maxiter;
      int verbose;
    };

  } // nampespace PDELab
} // namespace Dune

#endif
