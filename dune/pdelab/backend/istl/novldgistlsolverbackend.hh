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


    template<class GO, class GFS>
    class ISTLBackend_DGNOVLP_CG_NOPREC
      : public SequentialNorm, public LinearResultStorage
    {
      typedef istl::ParallelHelper<GFS> PHELPER;

    public:
      /*! \brief make a linear solver object

        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_DGNOVLP_CG_NOPREC (const GO& go_,
                                              const GFS& gfs_,
                                              unsigned maxiter_=5000,
                                              int verbose_=1)
        : go(go_), gfs(gfs_), phelper(gfs,verbose_), maxiter(maxiter_), verbose(verbose_)
      {}


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
