// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <dune/pdelab/finiteelementmap/pkfem.hh>

#ifdef USE_PK_FEM_FACTORY

struct PKFEMFactory
{

  template<typename GV, typename DF, typename RF, Dune::GeometryType::BasicType basic_type>
  struct FEM
  {
    typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,RF,FEM_FACTORY_ORDER,GV::dimension> type;
    typedef Dune::shared_ptr<type> pointer;
  };

  template<typename GV, typename DF, typename RF, Dune::GeometryType::BasicType basic_type>
  static typename FEM<GV,DF,RF,basic_type>::pointer create(const GV& gv)
  {
    return Dune::make_shared<typename FEM<GV,DF,RF,basic_type>::type>(gv);
  }

};

#endif
