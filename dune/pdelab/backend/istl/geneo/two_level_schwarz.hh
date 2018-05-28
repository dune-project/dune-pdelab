#ifndef TWO_LEVEL_SCHWARZ_HH
#define TWO_LEVEL_SCHWARZ_HH

#if HAVE_SUITESPARSE_UMFPACK

#include <dune/common/timer.hh>

#include "coarsespace.hh"

namespace Dune {
  namespace PDELab {
    namespace ISTL {

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
        TwoLevelOverlappingAdditiveSchwarz (const GFS& gfs, const M& AF, std::shared_ptr<CoarseSpace<X> > coarse_space, bool coarse_space_active = true, int verbosity = 0)
          : verbosity_(verbosity),
            coarse_space_active_(coarse_space_active),
            gfs_(gfs),
            solverf_(Dune::PDELab::Backend::native(AF),false),
            coarse_space_(coarse_space),
            coarse_solver_ (*coarse_space_->get_coarse_system()),
            coarse_defect_(coarse_space_->basis_size(), coarse_space_->basis_size()),
            prolongated_(gfs_, 0.0)
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
        virtual void apply (X& v, const Y& d)
        {
          // first the subdomain solves
          Y b(d); // need copy, since solver overwrites right hand side
          Dune::InverseOperatorResult result;
          solverf_.apply(v,b,result);

          if (!coarse_space_active_) {

            Dune::PDELab::AddDataHandle<GFS,X> adddh(gfs_,v);
            // Just add local results and return in 1-level Schwarz case
            gfs_.gridView().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);

          } else {

            gfs_.gridView().comm().barrier();
            Dune::Timer timer_coarse_solve;

            coarse_space_->restrict (d, coarse_defect_);

            // Solve coarse system
            Dune::InverseOperatorResult result;
            COARSE_V v0(coarse_space_->basis_size(),coarse_space_->basis_size());
            coarse_solver_.apply(v0, coarse_defect_, result);

            // Prolongate coarse solution on local domain
            coarse_space_->prolongate(v0, prolongated_);
            v += prolongated_;

            coarse_time_ += timer_coarse_solve.elapsed();
            apply_calls_++;

            Dune::PDELab::AddDataHandle<GFS,X> result_addh(gfs_,v);
            gfs_.gridView().communicate(result_addh,Dune::All_All_Interface,Dune::ForwardCommunication);
          }
        }

        /*!
          \brief Clean up.

          \copydoc Preconditioner::post(X&)
        */
        virtual void post (X& x) {
          if (verbosity_ > 0) std::cout << "Coarse time CT=" << coarse_time_ << std::endl;
          if (verbosity_ > 0) std::cout << "Coarse time per apply CTA=" << coarse_time_ / apply_calls_ << std::endl;
        }

      private:
        int verbosity_;
        bool coarse_space_active_;

        double coarse_time_ = 0.0;
        int apply_calls_ = 0;

        const GFS& gfs_;
        Dune::UMFPack<ISTLM> solverf_;
        std::shared_ptr<CoarseSpace<X> > coarse_space_;
        Dune::UMFPack<COARSE_M> coarse_solver_;

        typename CoarseSpace<X>::COARSE_V coarse_defect_;
        X prolongated_;
      };
    }
  }
}
#endif

#endif
