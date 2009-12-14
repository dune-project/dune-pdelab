// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_EDGES03DFEM_HH
#define DUNE_PDELAB_EDGES03DFEM_HH

#include<vector>
#include<dune/localfunctions/edges03d.hh>
#include"finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

	//! wrap up element from local functions
    //! \ingroup FiniteElementMap

    template<typename GV, typename R>
    class EdgeS03DLocalFiniteElementMap
      : public EdgeS0LocalFiniteElementMap<
          GV,
          EdgeS03DLocalFiniteElement<typename GV::Grid::ctype,R>,
          EdgeS03DLocalFiniteElementMap<GV,R>
        >
    {
      typedef EdgeS03DLocalFiniteElementMap<GV,R> This;
      typedef EdgeS03DLocalFiniteElement<typename GV::Grid::ctype,R> FE;
      typedef EdgeS0LocalFiniteElementMap<GV, FE, This> Base;
     
    public:
	  EdgeS03DLocalFiniteElementMap(const GV& gv) : Base(gv) {}
    };

  }
}

#endif // DUNE_PDELAB_EDGES03DFEM_HH
