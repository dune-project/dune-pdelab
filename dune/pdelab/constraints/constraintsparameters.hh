// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_CONSTRAINTSPARAMETERS_HH
#define DUNE_PDELAB_CONSTRAINTSPARAMETERS_HH

#include <dune/pdelab/common/typetree.hh>
#include <dune/common/fvector.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //! Interface for the constraints parameters describing dirichlet constraints
    struct DirichletConstraintsParameters :
      public TypeTree::LeafNode
    {
      template<typename I>
      bool isDirichlet(const I & intersection, const FieldVector<typename I::ctype, I::dimension-1> & coord)
      {
        return true;
      }
    };
    
    //! Interface for the constraints parameters describing flux constraints, e.g. needed for RT0
    struct FluxConstraintsParameters :
      public TypeTree::LeafNode
    {
      template<typename I>
      bool isNeumann(const I & intersection, const FieldVector<typename I::ctype, I::dimension-1> & coord)
      {
        return true;
      }
    };
        
    //! \}

  }
}

#endif // DUNE_PDELAB_CONSTRAINTSPARAMETERS_HH

