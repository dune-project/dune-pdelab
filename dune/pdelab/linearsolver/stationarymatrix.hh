// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_LINEARPROBLEM_STATIONARYMATRIX_HH
#define DUNE_PDELAB_LINEARPROBLEM_STATIONARYMATRIX_HH

#include <iostream>
#include <ostream>
#include <memory>

#include <dune/common/timer.hh>
#include <dune/common/deprecated.hh>

#warning <dune/pdelab/linearsolver/stationarymatrix.hh> and StationaryMatrixLinearSolver are deprecated and will be removed after PDELab 2.4. Please use LinearProblemSolver instead.

#include <dune/pdelab/backend/interface.hh>

namespace Dune {
  namespace PDELab {

    //! A class for solving linear problems with stationary matrices.
    /**
     * In apply() it first check whether the marix has already been assembled.
     * If it hasn't,it assembles the matrix and stores it for future
     * applications.  Then it computes the right hand side and solves the
     * problem.
     *
     * \tparam GOS   GridOperatorSpace to use.
     * \tparam SB    Solver backend.
     * \tparam Coeff Type of the matrix/vector entries
     */
    template<class GOS, class SB, class Coeff>
    class
    DUNE_DEPRECATED_MSG("StationaryMatrixLinearSolver is deprecated and will be removed after PDELab 2.4. Please use LinearProblemSolver instead.")
    StationaryMatrixLinearSolver
    {
      typedef typename GOS::template MatrixContainer<Coeff>::Type Matrix;
      using VectorU = Dune::PDELab::Backend::Vector
        <typename GOS::Traits::TrialGridFunctionSpace, Coeff>;
      using VectorV = Dune::PDELab::Backend::Vector
        <typename GOS::Traits::TestGridFunctionSpace, Coeff>;

      const GOS& gos;
      SB& sb;
      std::shared_ptr<Matrix> m;
      Coeff reduction;
      Coeff mindefect;

    public:
      StationaryMatrixLinearSolver(const GOS& gos_, SB& sb_, Coeff reduction_,
                                   Coeff mindefect_ = 1e-99) :
        gos(gos_), sb(sb_), reduction(reduction_), mindefect(mindefect_)
      { }

      void apply (VectorU& x) {
        Dune::Timer watch;
        double timing;

        if(!m) {
          // setup new matrix from sparsity pattern
          watch.reset();

          m.reset(new Matrix(gos));

          timing = watch.elapsed();
          if (gos.trialGridFunctionSpace().gridView().comm().rank()==0)
            std::cout << "=== matrix setup " << timing << " s" << std::endl;

          // assemble matrix
          watch.reset();

          *m = 0.0;
          gos.jacobian(x,*m);

          timing = watch.elapsed();
          if (gos.trialGridFunctionSpace().gridView().comm().rank()==0)
            std::cout << "=== matrix assembly " << timing << " s" << std::endl;
        }
        else {
          if (gos.trialGridFunctionSpace().gridView().comm().rank()==0)
            std::cout << "=== matrix setup skipped" << std::endl
                      << "=== matrix assembly skipped" << std::endl;
        }

        // assemble residual
        watch.reset();

        VectorV r(gos.testGridFunctionSpace(),0.0);
        gos.residual(x,r);  // residual is additive

        timing = watch.elapsed();
        if (gos.trialGridFunctionSpace().gridView().comm().rank()==0)
          std::cout << "=== residual assembly " << timing << " s" << std::endl;

        Coeff defect = sb.norm(r);

        // compute correction
        watch.reset();
        VectorU z(gos.trialGridFunctionSpace(),0.0);
        Coeff red = std::min(reduction,defect/mindefect);
        sb.apply(*m,z,r,red); // solver makes right hand side consistent
        timing = watch.elapsed();

        if (gos.trialGridFunctionSpace().gridView().comm().rank()==0)
          std::cout << "=== solving (reduction: " << red << ") "
                    << timing << " s" << std::endl;

        // and update
        x -= z;
      }

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LINEARPROBLEM_STATIONARYMATRIX_HH
