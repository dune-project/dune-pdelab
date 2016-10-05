#ifndef DUNE_PDELAB_BACKEND_ISTL_HH
#define DUNE_PDELAB_BACKEND_ISTL_HH

#include <dune/pdelab/backend/istl/vector.hh>
#include <dune/pdelab/backend/istl/bcrsmatrix.hh>
#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include <dune/pdelab/backend/istl/istlsolverbackend.hh>
#include <dune/pdelab/backend/istl/seq_amg_dg_backend.hh>
#include <dune/pdelab/backend/istl/ovlp_amg_dg_backend.hh>
#include <dune/pdelab/backend/istl/cg_to_dg_prolongation.hh>
#include <dune/pdelab/backend/istl/utility.hh>

/** \brief For backward compatibility -- Do not use this! */
namespace Dune {
  namespace PDELab {
    namespace istl {
      using namespace Dune::PDELab::ISTL;
    }
  }
}

#endif // DUNE_PDELAB_BACKEND_ISTL_HH
