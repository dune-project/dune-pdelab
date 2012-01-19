#ifndef DUNE_PDELAB_LOCALOPERATOR_STOKESPARAMETER_HH
#define DUNE_PDELAB_LOCALOPERATOR_STOKESPARAMETER_HH

#include <dune/pdelab/constraints/constraintsparameters.hh>
#include <dune/common/parametertree.hh>

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

      //! Interface for the parameter class required by the classes
      //! TaylorHoodNavierStokes and TaylorHoodNavierStokesJacobian.
      template <class BF, class NF, class RF>
      class TaylorHoodNavierStokesDefaultParameters
      {
      public:

        typedef StokesBoundaryCondition BCType;
        typedef BF BCTypeFunction;
        typedef NF NeumannStressFunction;

        TaylorHoodNavierStokesDefaultParameters(Dune::ParameterTree config, const BF & _bf, const NF & _nf):
          rho_(config.get<double>("rho")), 
          mu_(config.get<double>("mu")), 
          bf_(_bf), nf_(_nf)
        {}

        RF rho() const{ return rho_; }
        RF mu()  const{ return mu_; }

        template<typename IG, typename LC>
        typename BCType::Type
        bcType(const IG & ig, const LC & x) const
        {
          typename BCType::Type y;
          bf_.evaluate(ig,x,y);
          return y;
        }
        
        template<typename IG, typename LC, typename V>
        V
        stress(const IG & ig, const LC & x, V normal) const
        {
          typename NF::Traits::RangeType y(0);
          nf_.evaluate(*(ig.inside()),ig.geometryInInside().global(x),y);
          normal *= y;
          return normal;
        }

      private:
        const RF rho_;
        const RF mu_;
        const BF & bf_;
        const NF & nf_;
      };


        template<typename PRM>
        class StokesVelocityDirichletConstraints
            : public Dune::PDELab::DirichletConstraintsParameters
        {
        private:
            const PRM & prm_;

        public:

            StokesVelocityDirichletConstraints (const PRM & _prm)
              : prm_(_prm) { }

            template<typename I>
            bool isDirichlet(
                const I & intersection,
                const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                ) const
            {
                StokesBoundaryCondition::Type bctype =prm_.bcType(intersection,coord);
                return (bctype == StokesBoundaryCondition::VelocityDirichlet);
            }
        };

        template<typename PRM>
        class StokesPressureDirichletConstraints
            : public Dune::PDELab::DirichletConstraintsParameters
        {
        private:
            const PRM & prm_;

        public:

            StokesPressureDirichletConstraints (const PRM & _prm)
              : prm_(_prm) { }

            template<typename I>
            bool isDirichlet(
                const I & intersection,
                const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                ) const
            {
                StokesBoundaryCondition::Type bctype =prm_.bcType(intersection,coord);
                return (bctype == StokesBoundaryCondition::PressureDirichlet);
            }
        };
        
    }
}

#endif
