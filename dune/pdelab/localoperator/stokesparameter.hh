#ifndef DUNE_PDELAB_LOCALOPERATOR_STOKESPARAMETER_HH
#define DUNE_PDELAB_LOCALOPERATOR_STOKESPARAMETER_HH

#include <dune/pdelab/constraints/constraintsparameters.hh>

namespace Dune {
    namespace PDELab {

        /**
           These are the boundary condition types as to be returned by
           the employed boundary type function. 

           Possible types:

           <ul>

           <li>\a DoNothing : Do not evaluate boundary integrals.

           <li>\a VelocityDirichlet : Dirichlet conditions for velocity.

           <li>\a PressureDirichlet : Natural Neumann conditions for the
           impulse flux. These are equivalent to a fixed pressure
           condition \b if \f$ \forall i : n \cdot \nabla v_i = 0 \f$.

           </ul>
         */
        struct StokesBoundaryCondition {
            enum Type {
                DoNothing = 0,
                VelocityDirichlet = 1,
                PressureDirichlet = 2,
                SlipVelocity = 3
            };
        };
        
        template<typename BF>
        class StokesVelocityDirichletConstraints
            : public Dune::PDELab::DirichletConstraintsParameters
        {
        private:
            const BF & boundary_function_;
  
        public:

            StokesVelocityDirichletConstraints (const BF bf)
                : boundary_function_(bf)
            {
            }

            template<typename I>
            bool isDirichlet(
                const I & intersection,
                const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                ) const
            {
                StokesBoundaryCondition::Type bctype;
                boundary_function_.evaluate(intersection,coord,bctype);

                return (bctype ==
                    StokesBoundaryCondition::VelocityDirichlet);
            }
            
        };
        
        template<typename BF>
        class StokesPressureDirichletConstraints
            : public Dune::PDELab::DirichletConstraintsParameters
        {
        private:
            const BF & boundary_function_;
  
        public:

            StokesPressureDirichletConstraints (const BF bf)
                : boundary_function_(bf)
            {
            }

            template<typename I>
            bool isDirichlet(
                const I & intersection,
                const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                ) const
            {
                StokesBoundaryCondition::Type bctype;
                boundary_function_.evaluate(intersection,coord,bctype);

                return (bctype ==
                    StokesBoundaryCondition::PressureDirichlet);
            }
            
        };
        
        
    }
}

#endif
