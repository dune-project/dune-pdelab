#ifndef TWO_LEVEL_SCHWARZ_HH
#define TWO_LEVEL_SCHWARZ_HH

#if HAVE_ARPACKPP && HAVE_SUITESPARSE_UMFPACK

#include <dune/pdelab/boilerplate/pdelab.hh>

#include <dune/common/timer.hh>

#include "coarsespace.hh"

namespace Dune {
  namespace PDELab {

    /*!
    * \brief Two level overlapping Schwarz preconditioner with arbitrary coarse space.
    */
    template<class GFS, class M, class X, class Y>
    class TwoLevelOverlappingAdditiveSchwarz
      : public Dune::Preconditioner<X,Y>
    {
    public:
      typedef Dune::PDELab::Backend::Native<M> ISTLM;

      typedef Dune::BlockVector<Dune::FieldVector<double,1> > COARSE_V;
      typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > COARSE_M;

      // define the category
      virtual Dune::SolverCategory::Category category() const
      {
        return Dune::SolverCategory::overlapping;
      }

      /*! \brief Constructor.

        Constructor gets all parameters to operate the prec.
        \param A The matrix to operate on.
        \param n The number of iterations to perform.
        \param w The relaxation factor.
      */
      TwoLevelOverlappingAdditiveSchwarz (const GFS& gfs_, const M& AF, std::shared_ptr<CoarseSpace<M,X> > coarse_space, bool coarse_space_active = true)
        : gfs(gfs_),
          solverf(Dune::PDELab::Backend::native(AF),false),
          my_rank(gfs.gridView().comm().rank()),
          coarse_space_(coarse_space),
          coarse_solver_ (*coarse_space_->get_coarse_system()),
          coarse_space_active(coarse_space_active)
      { }

      /*!
        \brief Prepare the preconditioner.

        \copydoc Preconditioner::pre(X&,Y&)
      */
      virtual void pre (X& x, Y& b)
      { }

      /*!
        \brief Apply the precondioner.

        \copydoc Preconditioner::apply(X&,const Y&)
      */
      double coarse_time = 0.0;
      int apply_calls = 0;
      bool coarse_space_active = true;

      virtual void apply (X& v, const Y& d)
      {
        // first the subdomain solves
        Y b(d); // need copy, since solver overwrites right hand side
        Dune::InverseOperatorResult result;
        solverf.apply(v,b,result);

        if (!coarse_space_active) {

          Dune::PDELab::AddDataHandle<GFS,X> adddh(gfs,v);
          // Just add local results and return in 1-level Schwarz case
          gfs.gridView().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);

        } else {

          MPI_Barrier(gfs.gridView().comm());
          Dune::Timer timer_setup;

          // coarse defect
          auto coarse_defect = coarse_space_->restrict_defect (d);

          // Solve coarse system
          Dune::InverseOperatorResult result;
          COARSE_V v0(coarse_space_->basis_size(),coarse_space_->basis_size());
          coarse_solver_.apply(v0, *coarse_defect, result);

          // Prolongate coarse solution on local domain
          auto coarse_correction = coarse_space_->prolongate_defect (v0);
          v += *coarse_correction;

          coarse_time += timer_setup.elapsed();
          apply_calls++;

          Dune::PDELab::AddDataHandle<GFS,X> result_addh(gfs,v);
          gfs.gridView().communicate(result_addh,Dune::All_All_Interface,Dune::ForwardCommunication);
        }
      }

      /*!
        \brief Clean up.

        \copydoc Preconditioner::post(X&)
      */
      virtual void post (X& x) {
        if (my_rank == 0) std::cout << "Coarse time CT=" << coarse_time << std::endl;
        if (my_rank == 0) std::cout << "Coarse time per apply CTA=" << coarse_time / (double)apply_calls << std::endl;
      }

    private:

      const GFS& gfs;
      Dune::UMFPack<ISTLM> solverf;
      std::shared_ptr<CoarseSpace<M,X> > coarse_space_;
      Dune::UMFPack<COARSE_M> coarse_solver_;

      int my_rank;
    };

  }
}
#endif

#endif
