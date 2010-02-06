// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_EDGES02DFEM_HH
#define DUNE_PDELAB_EDGES02DFEM_HH

#include<vector>
#include<dune/localfunctions/whitney/edges02d.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

	//! wrap up element from local functions
    //! \ingroup FiniteElementMap

    template<typename GV, typename R>
    class EdgeS02DLocalFiniteElementMap
      : public EdgeS0LocalFiniteElementMap<
          GV,
          EdgeS02DLocalFiniteElement<typename GV::Grid::ctype,R>,
          EdgeS02DLocalFiniteElementMap<GV,R>
        >
    {
      typedef EdgeS02DLocalFiniteElementMap<GV,R> This;
      typedef EdgeS02DLocalFiniteElement<typename GV::Grid::ctype,R> FE;
      typedef EdgeS0LocalFiniteElementMap<GV, FE, This> Base;
     
    public:
	  EdgeS02DLocalFiniteElementMap(const GV& gv) : Base(gv) {}
    };

  }
}

#endif // DUNE_PDELAB_EDGES02DFEM_HH
