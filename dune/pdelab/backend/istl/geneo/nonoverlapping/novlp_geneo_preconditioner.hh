#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_NOVLP_GENEO_PRECONDITIONER_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_NOVLP_GENEO_PRECONDITIONER_HH

#include <dune/pdelab/backend/istl/geneo/nonoverlapping/geneobasis.hh>

namespace Dune {
  namespace PDELab {

    /**
     * @brief Black-box GenEO preconditioner working on nonoverlapping matrices
     */
    template <typename GO, typename Matrix, typename Vector>
    class NonoverlappingGenEOPreconditioner : public Dune::Preconditioner<Vector,Vector> {

      using GFS = typename GO::Traits::TrialGridFunctionSpace;
      using GV = typename GFS::Traits::GridView;


    public:

      // define the category
      virtual Dune::SolverCategory::Category category() const override
      {
        return prec->category();
      }


      /*!
      *          \brief Prepare the preconditioner.
      *
      *          \copydoc Preconditioner::pre(X&,Y&)
      */
      virtual void pre (Vector& x, Vector& b) override
      {
        prec->pre(x, b);
      }

      /*!
      *          \brief Apply the precondioner.
      *
      *          \copydoc Preconditioner::apply(X&,const Y&)
      */
      virtual void apply (Vector& v, const Vector& d) override
      {
        prec->apply(v, d);
      }
      /*!
      *          \brief Clean up.
      *
      *          \copydoc Preconditioner::post(X&)
      */
      virtual void post (Vector& x) override
      {
        prec->post(x);
      }

      NonoverlappingGenEOPreconditioner(const GO& go, const typename GO::Jacobian& A, int algebraic_overlap, int avg_nonzeros, const double eigenvalue_threshold, int& nev,
                          int nev_arpack = -1, const double shift = 0.001, int verbose = 0)
      {
        using Dune::PDELab::Backend::native;

        const GV& gv = go.trialGridFunctionSpace().gridView();

        Dune::NonoverlappingOverlapAdapter<GV, Vector, Matrix> adapter(gv, native(A), avg_nonzeros, algebraic_overlap);

        auto geneo_matrices = setupGenEOMatrices(go, adapter, A);
        std::shared_ptr<Matrix> A_extended = std::get<0>(geneo_matrices);
        std::shared_ptr<Matrix> A_overlap_extended = std::get<1>(geneo_matrices);
        std::shared_ptr<Vector> part_unity = std::get<2>(geneo_matrices);

        auto subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GO, Matrix, Vector>>(adapter, A_extended, A_overlap_extended, part_unity, eigenvalue_threshold, nev, nev_arpack, shift);

        auto coarse_space = std::make_shared<Dune::PDELab::NonoverlappingSubdomainProjectedCoarseSpace<GV, Matrix, Vector>>(adapter, gv, *A_extended, subdomainbasis, verbose);

        prec = std::make_shared<Dune::PDELab::ISTL::NonoverlappingTwoLevelOverlappingAdditiveSchwarz<GV, Matrix, Vector>>(adapter, A_extended, *part_unity, coarse_space, true, verbose);

      }

    private:

      std::shared_ptr<Dune::PDELab::ISTL::NonoverlappingTwoLevelOverlappingAdditiveSchwarz<GV, Matrix, Vector>> prec = nullptr;

    };

  }
}

#endif
