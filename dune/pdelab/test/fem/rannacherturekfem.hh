// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_TEST_FEM_RANNACHERTUREKFEM_HH
#define DUNE_PDELAB_TEST_FEM_RANNACHERTUREKFEM_HH

#include <dune/pdelab/finiteelementmap/rannacherturekfem.hh>

#ifdef USE_RANNACHER_TUREK_FEM_FACTORY

struct RannacherTurekFEMFactory
{

  template<typename GV, typename DF, typename RF, Dune::GeometryType::BasicType basic_type>
  struct FEM
  {
    typedef Dune::PDELab::RannacherTurekLocalFiniteElementMap<DF,RF,GV::dimension> type;
    typedef std::shared_ptr<type> pointer;
  };

  template<typename GV, typename DF, typename RF, Dune::GeometryType::BasicType basic_type>
  static typename FEM<GV,DF,RF,basic_type>::pointer create(const GV& gv)
  {
    return Dune::make_shared<typename FEM<GV,DF,RF,basic_type>::type>();
  }

};

#endif

#endif // DUNE_PDELAB_TEST_FEM_RANNACHERTUREKFEM_HH
