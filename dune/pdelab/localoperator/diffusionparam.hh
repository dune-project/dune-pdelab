#ifndef DUNE_PDELAB_DIFFUSIONPARAM_HH
#define DUNE_PDELAB_DIFFUSIONPARAM_HH

#include <iostream>
#include <dune/common/deprecated.hh>

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    /** \brief Group types and methods to specify the boundary condition of a diffusion problem
     */
    struct DiffusionBoundaryCondition {
      //! Enum for the Boundary condition type of a diffusion problem
      enum Type {
        //! Neumann boundary condition (prescribed flux)
        Neumann = 0,
        //! Dirichlet boundary condition (prescribed value)
        Dirichlet = 1
      };

      //! Test for Dirichlet boundary condition
      static bool isDirichlet (Type i)
      {
        return (i == Dirichlet);
      }
      //! Test for Neumann boundary condition
      static bool isNeumann (Type i)
      {
        return (i == Neumann);
      }
    };
    //! \}
  }
}

#endif // DUNE_PDELAB_DIFFUSIONPARAM_HH
