#ifndef DUNE_PDELAB_DIFFUSIONPARAM_HH
#define DUNE_PDELAB_DIFFUSIONPARAM_HH

#include <iostream>
#include <dune/common/deprecated.hh>

namespace Dune {
  namespace PDELab {
    struct DiffusionBoundaryCondition {
      enum Type {
        Neumann = 0,
        Dirichlet = 1
      };
      static bool isDirichlet (int i) DUNE_DEPRECATED
      {
        static const std::ostream & warning =
          std::cerr << "Warning! Please update your code "
                    << "to use the DiffusionBoundaryCondition enum" << std::endl;
        return (i>0);
      }
      static bool isNeumann (int i) DUNE_DEPRECATED
      {
        static const std::ostream & warning =
          std::cerr << "Warning! Please update your code "
                    << "to use the DiffusionBoundaryCondition enum" << std::endl;
        return (i<1);
      }
      static bool isDirichlet (Type i)
      {
        return (i == Dirichlet);
      }
      static bool isNeumann (Type i)
      {
        return (i == Neumann);
      }
    };
  }
}

#endif // DUNE_PDELAB_DIFFUSIONPARAM_HH
