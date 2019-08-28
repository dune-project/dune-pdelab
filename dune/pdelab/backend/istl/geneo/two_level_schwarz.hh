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
      template<class GFS, class M, class M_EXTERIOR, class X, class Y>
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
        TwoLevelOverlappingAdditiveSchwarz (const GFS& gfs, const M& AF, const M_EXTERIOR& AF_exterior, std::shared_ptr<CoarseSpace<X> > coarse_space, const X& part_unity, bool restrictedMode = false, bool coarse_space_active = true, bool hybrid = false, int verbosity = 0)
          : verbosity_(verbosity),
            hybrid_(hybrid),
            coarse_space_active_(coarse_space_active),
            gfs_(gfs),
            AF_(AF),
            AF_exterior_(AF_exterior),
            solverf_(Dune::PDELab::Backend::native(AF),false),
            coarse_space_(coarse_space),
            part_unity_(part_unity),
            restrictedMode_(restrictedMode),
            coarse_solver_ (*coarse_space_->get_coarse_system()),
            coarse_defect_(coarse_space_->basis_size(), coarse_space_->basis_size()),
            prolongated_(gfs_, 0.0),
            b_(gfs_, 0.0)
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

          if (!coarse_space_active_) {
            // first the subdomain solves
            b_ = d; // need copy, since solver overwrites right hand side
            Dune::InverseOperatorResult result;
            solverf_.apply(v,b_,result);

            Dune::PDELab::AddDataHandle<GFS,X> adddh(gfs_,v);

            if (restrictedMode_) {
              std::transform(
                v.begin(),v.end(),
                             part_unity_.begin(),
                v.begin(),
                std::multiplies<>()
              );
            }

            // Just add local results and return in 1-level Schwarz case
            gfs_.gridView().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);

          } else if (!hybrid_) {
            // first the subdomain solves
            b_ = d; // need copy, since solver overwrites right hand side
            Dune::InverseOperatorResult result;
            solverf_.apply(v,b_,result);

            if (restrictedMode_) {
              using Dune::PDELab::Backend::native;
              for (int j = 0; j < v.N(); j++) {
                for(int j_block = 0; j_block < ISTLM::block_type::rows; j_block++){
                  if (native(part_unity_)[j][j_block] == 0.0)
                    native(v)[j][j_block] *= .5;
                  else
                    native(v)[j][j_block] *= native(part_unity_)[j][j_block];
                }
              }
            }

            //gfs_.gridView().comm().barrier();
            Dune::Timer timer_coarse_solve;

            coarse_space_->restrict (d, coarse_defect_);

            // Solve coarse system
            COARSE_V v0(coarse_space_->basis_size(),coarse_space_->basis_size());
            coarse_solver_.apply(v0, coarse_defect_, result);

            // Prolongate coarse solution on local domain
            coarse_space_->prolongate(v0, prolongated_);
            v += prolongated_;

            coarse_time_ += timer_coarse_solve.elapsed();

            Dune::PDELab::AddDataHandle<GFS,X> result_addh(gfs_,v);
            gfs_.gridView().communicate(result_addh,Dune::All_All_Interface,Dune::ForwardCommunication);
          } else {

            using PDELab::Backend::native;

            // first the subdomain solves
            b_ = d; // need copy, since solver overwrites right hand side
            Dune::InverseOperatorResult result;
            solverf_.apply(v,b_,result);

            if (restrictedMode_) {
              std::transform(
                v.begin(),v.end(),
                             part_unity_.begin(),
                             v.begin(),
                             std::multiplies<>()
              );
            }

            Dune::PDELab::AddDataHandle<GFS,X> corr_addh(gfs_,v);
            gfs_.gridView().communicate(corr_addh,Dune::All_All_Interface,Dune::ForwardCommunication);



            b_ = d;
            native(AF_).mmv(native(v), native(b_));




            //gfs_.gridView().comm().barrier();
            Dune::Timer timer_coarse_solve;

            coarse_space_->restrict (b_, coarse_defect_);

            // Solve coarse system
            COARSE_V v0(coarse_space_->basis_size(),coarse_space_->basis_size());
            coarse_solver_.apply(v0, coarse_defect_, result);

            // Prolongate coarse solution on local domain
            coarse_space_->prolongate(v0, prolongated_);

            Dune::PDELab::AddDataHandle<GFS,X> result_addh(gfs_,prolongated_);
            gfs_.gridView().communicate(result_addh,Dune::All_All_Interface,Dune::ForwardCommunication);

            v += prolongated_;

           coarse_time_ += timer_coarse_solve.elapsed();

           // Second fine solve in order to symmetrize
           /*native(AF_).mmv(native(v), native(c));

           Y v2(d);
           solverf_.apply(v2,c,result);

           Dune::PDELab::AddDataHandle<GFS,X> result_addh2(gfs_,v2);
           gfs_.gridView().communicate(result_addh2,Dune::All_All_Interface,Dune::ForwardCommunication);

           v += v2;*/

          }
          apply_calls_++;
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
        bool hybrid_;

        double coarse_time_ = 0.0;
        int apply_calls_ = 0;

        bool restrictedMode_;
        const X& part_unity_;

        const GFS& gfs_;
        Dune::UMFPack<ISTLM> solverf_;
        std::shared_ptr<CoarseSpace<X> > coarse_space_;
        Dune::UMFPack<COARSE_M> coarse_solver_;

        const M_EXTERIOR& AF_exterior_;
        const M& AF_;

        typename CoarseSpace<X>::COARSE_V coarse_defect_;
        X prolongated_;
        Y b_;
      };
    }
  }
}
#endif

#endif
