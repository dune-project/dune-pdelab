// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_TEST_FEM_OPBFEM_HH
#define DUNE_PDELAB_TEST_FEM_OPBFEM_HH

#include <dune/pdelab/finiteelementmap/opbfem.hh>

#ifdef USE_OPB_FEM_FACTORY

struct OPBFEMFactory
{

  template<typename GV, typename DF, typename RF, Dune::GeometryType::BasicType basic_type>
  struct FEM
  {

#if HAVE_GMP
    typedef Dune::GMPField<512> CFT;
#else
#warning Testing OPBLocalFiniteElementMap without GMP!
    typedef double CFT;
#endif

    typedef Dune::PDELab::OPBLocalFiniteElementMap<DF,RF,FEM_FACTORY_ORDER,GV::dimension,basic_type,CFT> type;
    typedef std::shared_ptr<type> pointer;
  };

  template<typename GV, typename DF, typename RF, Dune::GeometryType::BasicType basic_type>
  static typename FEM<GV,DF,RF,basic_type>::pointer create(const GV& gv)
  {
#if !HAVE_GMP
    std::cerr << "Warning: Testing OPBLocalFiniteElementMap without GMP!" << std::endl;
#endif
    return Dune::make_shared<typename FEM<GV,DF,RF,basic_type>::type>();
  }

};

#endif

#endif // DUNE_PDELAB_TEST_FEM_OPBFEM_HH
