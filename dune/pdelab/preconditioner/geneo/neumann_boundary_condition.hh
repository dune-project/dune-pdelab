#ifndef DUNE_PDELAB_NEUMANN_BOUNDARY_CONDITION_HH
#define DUNE_PDELAB_NEUMANN_BOUNDARY_CONDITION_HH

namespace Dune {
  namespace PDELab {

    /*!
     * \brief Exclusive Neumann boundary condition.
     */
    class PureNeumannBoundaryCondition
      :
      public Dune::PDELab::FluxConstraintsParameters,
      public Dune::PDELab::DirichletConstraintsParameters   /*@\label{bcp:base}@*/
    {
    public:

      PureNeumannBoundaryCondition()
      {}

      template<typename I>
      bool isDirichlet(const I & ig               /*@\label{bcp:name}@*/
                       , const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                       ) const
      {
        return false;
      }

      template<typename I>
      bool isNeumann(const I & ig,   /*@\label{bcp:name}@*/
                     const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                     ) const
      {
        return true;
      }

    };
  }
}

#endif //DUNE_PDELAB_NEUMANN_BOUNDARY_CONDITION_HH
