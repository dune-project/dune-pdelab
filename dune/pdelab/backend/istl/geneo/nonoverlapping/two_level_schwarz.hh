#ifndef NONOVERLAPPING_TWO_LEVEL_SCHWARZ_HH
#define NONOVERLAPPING_TWO_LEVEL_SCHWARZ_HH

#if HAVE_SUITESPARSE_UMFPACK

#include <dune/common/timer.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/umfpack.hh>
#include <dune/pdelab/gridfunctionspace/dunefunctionsgridfunctionspace.hh>

#include <dune/pdelab/backend/istl/geneo/coarsespace.hh>
#include "overlaptools.hh"

#include <dune/common/parallel/communicator.hh>

namespace Dune {



  namespace PDELab {
    namespace ISTL {

      /*!
      * \brief Two level overlapping Schwarz preconditioner with arbitrary coarse space.
      */
      //template<class GFS, class M, class X, class Y>
      template<typename GridView, typename Matrix, typename Vector>
      class NonoverlappingTwoLevelOverlappingAdditiveSchwarz
       : public Dune::Preconditioner<Vector,Vector>
      {

        template<typename V>
        struct AddGatherScatter
        {
          static typename V::value_type gather(const V& a, int i)
          {
            return a[i]; // I am sending my value
          }
          static void scatter (V& a, typename V::value_type v, int i)
          {
            a[i]+=v; // add what I receive to my value
          }
        };

      public:
        //typedef Dune::PDELab::Backend::Native<M> ISTLM;

        //typedef Dune::BlockVector<Dune::FieldVector<double,1> > COARSE_V;
        //typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > COARSE_M;

        // define the category
        virtual Dune::SolverCategory::Category category() const
        {
          return Dune::SolverCategory::nonoverlapping;
        }

        NonoverlappingTwoLevelOverlappingAdditiveSchwarz (const NonoverlappingOverlapAdapter<GridView, Vector, Matrix>& adapter, std::shared_ptr<Matrix> A, const Vector& part_unity, std::shared_ptr<CoarseSpace<Vector> > coarse_space, bool coarse_space_active = true, int verbosity = 0)
        : verbosity_(verbosity),
          coarse_space_active_(coarse_space_active),
          adapter_(adapter),
          A_(A),
          coarse_space_(coarse_space),
          coarse_solver_ (*coarse_space_->get_coarse_system()),
          coarse_defect_(coarse_space_->basis_size(), coarse_space_->basis_size()),
          prolongated_(adapter_.getExtendedSize())
        {
          const int block_size = Vector::block_type::dimension;

          // Apply Dirichlet conditions to matrix on processor boundaries, inferred from partition of unity
          for (auto rIt=A_->begin(); rIt!=A_->end(); ++rIt){
            for(int block_i = 0; block_i < block_size; block_i++){
              if (part_unity[rIt.index()][block_i] == .0) {
                for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
                {
                  for(int block_j = 0; block_j < block_size; block_j++){
                    (*cIt)[block_i][block_j] = (rIt.index() == cIt.index() && block_i == block_j) ? 1.0 : 0.0;
                  }
                }
              }
            }
          }

          solverf_ = std::make_shared<Dune::UMFPack<Matrix>>(*A_,false);
        }

        /*!
          \brief Prepare the preconditioner.

          \copydoc Preconditioner::pre(X&,Y&)
        */
        virtual void pre (Vector& x, Vector& b)
        { }

        /*!
          \brief Apply the precondioner.

          \copydoc Preconditioner::apply(X&,const Y&)
        */
        virtual void apply (Vector& v, const Vector& d)
        {

          using Attribute = EPISAttribute;
          Dune::AllSet<Attribute> allAttribute;
          auto allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
          allinterface->build(*adapter_.getRemoteIndices(),allAttribute,allAttribute); // all to all communication

          // build up buffered communicator allowing communication over a dof vector
          auto communicator = std::shared_ptr<Dune::BufferedCommunicator>(new Dune::BufferedCommunicator());
          communicator->build<Vector>(*allinterface);

          if (verbosity_ > 2) Dune::printvector(std::cout, d, "defect (local)", "", 1, 10, 17);

          // first the subdomain solves
          Vector b(adapter_.getExtendedSize()), correction(adapter_.getExtendedSize()); // need copy, since solver overwrites right hand side
          adapter_.extendVector(d, b);
          if (verbosity_ > 2) Dune::printvector(std::cout, b, "defect (extended)", "", 1, 10, 17);
          communicator->forward<AddGatherScatter<Vector>>(b,b); // make function known in other subdomains
          if (verbosity_ > 2) Dune::printvector(std::cout, b, "defect (distributed)", "", 1, 10, 17);


          const int block_size = Vector::block_type::dimension;
          Vector b_cpy(b);// FIXME: Avoid this
          for (auto rIt=A_->begin(); rIt!=A_->end(); ++rIt) {
              for(int block_i = 0; block_i < block_size; block_i++){
                  bool isDirichlet = true;
                  for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
                  {
                      for(int block_j = 0; block_j < block_size; block_j++){
                          if ((rIt.index() != cIt.index() || block_i!=block_j) && (*cIt)[block_i][block_j] != 0.0) {
                              isDirichlet = false;
                              break;
                          }
                      }
                      if(!isDirichlet) break;
                  }
                  if (isDirichlet) {
                      b_cpy[rIt.index()] = .0;
                      b[rIt.index()] = .0;
                  }
              }
          }
          if (verbosity_ > 2) Dune::printvector(std::cout, b_cpy, "defect (Dirichlet applied) ", "", 1, 10, 17);

          Dune::InverseOperatorResult result;
          solverf_->apply(correction,b_cpy,result);

          if (verbosity_ > 2) Dune::printvector(std::cout, correction, "correction (1lvl) ", "", 1, 10, 17);

          if (!coarse_space_active_) {


            communicator->forward<AddGatherScatter<Vector>>(correction,correction); // make function known in other subdomains
            std::cout << "size correction: " << correction.N() << " size v: " << v.N() << std::endl;

            adapter_.restrictVector(correction, v);

            //Dune::PDELab::AddDataHandle<GFS,X> adddh(gfs_,v);
            // Just add local results and return in 1-level Schwarz case
            //gfs_.gridView().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);

          } else {

            //gfs_.gridView().comm().barrier();
            Dune::Timer timer_coarse_solve;

            coarse_space_->restrict (b, coarse_defect_);

            if (verbosity_ > 2) Dune::printvector(std::cout, coarse_defect_, "coarse_defect_ ", "", 1, 10, 17);

            typename CoarseSpace<Vector>::COARSE_V v0(coarse_space_->basis_size());
            coarse_solver_.apply(v0, coarse_defect_, result);
            if (verbosity_ > 2) Dune::printvector(std::cout, v0, "v0 ", "");
            coarse_space_->prolongate(v0, prolongated_);

            if (verbosity_ > 2) Dune::printvector(std::cout, prolongated_, "prolongated_ ", "", 1, 10, 17);

            correction += prolongated_;
            if (verbosity_ > 2) Dune::printvector(std::cout, correction, "correction ", "", 1, 10, 17);

            communicator->forward<AddGatherScatter<Vector>>(correction,correction); // make function known in other subdomains
            if (verbosity_ > 2) Dune::printvector(std::cout, correction, "correction (sum) ", "", 1, 10, 17);

            adapter_.restrictVector(correction, v);
            if (verbosity_ > 2) Dune::printvector(std::cout, v, "correction (restricted) ", "", 1, 10, 17);

            //exit(0);

            //communicator->forward<AddGatherScatter<Vector>>(correction,correction); // make function known in other subdomains
            //adapter_.restrictVector(correction, v);


            // Solve coarse system
            /*Dune::InverseOperatorResult result;
            COARSE_V v0(coarse_space_->basis_size(),coarse_space_->basis_size());
            coarse_solver_.apply(v0, coarse_defect_, result);

            // Prolongate coarse solution on local domain
            coarse_space_->prolongate(v0, prolongated_);
            v += prolongated_;*/

            coarse_time_ += timer_coarse_solve.elapsed();

            //Dune::PDELab::AddDataHandle<GFS,X> result_addh(gfs_,v);
            //gfs_.gridView().communicate(result_addh,Dune::All_All_Interface,Dune::ForwardCommunication);
          }
          apply_calls_++;
        }

        /*!
          \brief Clean up.

          \copydoc Preconditioner::post(X&)
        */
        virtual void post (Vector& x) {
          if (verbosity_ > 0) std::cout << "Coarse time CT=" << coarse_time_ << std::endl;
          if (verbosity_ > 0) std::cout << "Coarse time per apply CTA=" << coarse_time_ / apply_calls_ << std::endl;
        }

      private:
        int verbosity_;
        bool coarse_space_active_;

        NonoverlappingOverlapAdapter<GridView, Vector, Matrix> adapter_;

        std::shared_ptr<Matrix> A_ = nullptr;
        std::shared_ptr<Dune::UMFPack<Matrix>> solverf_ = nullptr;

        double coarse_time_ = 0.0;
        int apply_calls_ = 0;

        std::shared_ptr<CoarseSpace<Vector> > coarse_space_;
        Dune::UMFPack<typename CoarseSpace<Vector>::COARSE_M> coarse_solver_;

        typename CoarseSpace<Vector>::COARSE_V coarse_defect_;
        Vector prolongated_;

      };
    }
  }
}
#endif

#endif
