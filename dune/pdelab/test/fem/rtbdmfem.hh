// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_TEST_FEM_RTBDMFEM_HH
#define DUNE_PDELAB_TEST_FEM_RTBDMFEM_HH

#include <dune/pdelab/finiteelementmap/brezzidouglasmarinifem.hh>
#include <dune/pdelab/finiteelementmap/raviartthomasfem.hh>

#ifdef USE_RT_BDM_FEM_FACTORY

struct RTBDMFEMFactory
{

  template<typename GV, typename DF, typename RF, Dune::GeometryType::BasicType basic_type>
  struct FEM
  {
    typedef FEM_FACTORY_FEM_CLASS<GV,DF,RF,FEM_FACTORY_ORDER,basic_type> type;
    typedef std::shared_ptr<type> pointer;
  };

  template<typename GV, typename DF, typename RF, Dune::GeometryType::BasicType basic_type>
  static typename FEM<GV,DF,RF,basic_type>::pointer create(const GV& gv)
  {
    return Dune::make_shared<typename FEM<GV,DF,RF,basic_type>::type>(gv);
  }

};

#endif

#endif // DUNE_PDELAB_TEST_FEM_RTBDMFEM_HH
