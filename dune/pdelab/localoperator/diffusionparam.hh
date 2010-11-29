#ifndef DUNE_PDELAB_DIFFUSIONPARAM_HH
#define DUNE_PDELAB_DIFFUSIONPARAM_HH

#include <iostream>
#include <dune/common/deprecated.hh>
#include <dune/common/ftraits.hh>

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
      static bool isDirichlet (int i) DUNE_DEPRECATED
      {
        static const std::ostream & warning =
          std::cerr << "Warning! Please update your code "
                    << "to use the DiffusionBoundaryCondition enum" << std::endl;
        return (i>0);
      }
      //! Test for Neumann boundary condition
      static bool isNeumann (int i) DUNE_DEPRECATED
      {
        static const std::ostream & warning =
          std::cerr << "Warning! Please update your code "
                    << "to use the DiffusionBoundaryCondition enum" << std::endl;
        return (i<1);
      }
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
  template<>
  struct FieldTraits<PDELab::DiffusionBoundaryCondition::Type>
  {
    typedef PDELab::DiffusionBoundaryCondition::Type field_type;
    typedef PDELab::DiffusionBoundaryCondition::Type real_type;
};

}

#endif // DUNE_PDELAB_DIFFUSIONPARAM_HH
