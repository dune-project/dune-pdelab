#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_GENEO_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_GENEO_HH


#include <dune/pdelab/backend/istl/geneo/two_level_schwarz.hh>
#include <dune/pdelab/backend/istl/geneo/nonoverlapping/two_level_schwarz.hh>

#include <dune/pdelab/backend/istl/geneo/subdomainprojectedcoarsespace.hh>
#include <dune/pdelab/backend/istl/geneo/nonoverlapping/subdomainprojectedcoarsespace.hh>

#include <dune/pdelab/backend/istl/geneo/partitionofunity.hh>
#include <dune/pdelab/backend/istl/geneo/localoperator_ovlp_region.hh>
#include <dune/pdelab/backend/istl/geneo/nonoverlapping/geneo_matrix_setup.hh>
#include <dune/pdelab/backend/istl/geneo/nonoverlapping/novlp_operators.hh>
#include <dune/pdelab/backend/istl/geneo/nonoverlapping/novlp_geneo_preconditioner.hh>
#include <dune/pdelab/backend/istl/geneo/nonoverlapping/novlp_geneo_preconditioner_fromfiles.hh>
#include <dune/pdelab/backend/istl/geneo/nonoverlapping/communicator_with_rank.hh>
#include <dune/pdelab/backend/istl/geneo/nonoverlapping/variablesizecommunicator_with_rank.hh>

#include <dune/pdelab/backend/istl/geneo/geneobasis.hh>
#include <dune/pdelab/backend/istl/geneo/nonoverlapping/geneobasis.hh>
#include <dune/pdelab/backend/istl/geneo/liptonbabuskabasis.hh>


#endif // DUNE_PDELAB_BACKEND_ISTL_GENEO_GENEO_HH
