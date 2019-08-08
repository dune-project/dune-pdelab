// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_EIGEN_HH
#define DUNE_PDELAB_BACKEND_EIGEN_HH

#if HAVE_EIGEN

#include <dune/pdelab/backend/eigen/descriptors.hh>
#include <dune/pdelab/backend/eigen/vector.hh>
#include <dune/pdelab/backend/eigen/matrix.hh>
#include <dune/pdelab/backend/eigen/solvers.hh>

/** \brief For backward compatibility -- Do not use this! */
namespace Dune {
  namespace PDELab {
    namespace EIGEN {
      using namespace Eigen;
    }
  }
}
#endif //HAVE_EIGEN
#endif // DUNE_PDELAB_BACKEND_EIGEN_HH
