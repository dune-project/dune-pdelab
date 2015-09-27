// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_SOLVER_HH
#define DUNE_PDELAB_BACKEND_SOLVER_HH

#include <dune/common/fvector.hh>
#include <dune/pdelab/backend/interface.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Backend
    //! \ingroup PDELab
    //! \{

    struct SequentialNorm
    {/*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      template<class V>
      typename Dune::template FieldTraits<typename V::ElementType >::real_type norm(const V& v) const
      {
        return Backend::native(v).two_norm();
      }
    };

    // Status information of a linear solver
    template<class RFType>
    struct LinearSolverResult
    {
      bool converged;            // Solver converged
      unsigned int iterations;   // number of iterations
      double elapsed;            // total user time in seconds
      RFType reduction;          // defect reduction
      RFType conv_rate;          // convergence rate (average reduction per step)

      LinearSolverResult() :
        converged(false), iterations(0), elapsed(0.0), reduction(0.0), conv_rate(0.0) {}
    };

    class LinearResultStorage
    {
    public:
      /*! \brief Return access to result data */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res;
      }

    protected:
      Dune::PDELab::LinearSolverResult<double> res;
    };

    //! \} group Backend

  } // end namespace PDELab
} // end namespace Dune


#endif // DUNE_PDELAB_BACKEND_SOLVER_HH
