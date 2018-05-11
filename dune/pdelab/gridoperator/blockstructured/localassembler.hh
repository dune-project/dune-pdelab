//
// Created by marckoch on 09/05/18.
//

#ifndef DUNE_BLOCKSTRUCTURED_LOCALASSEMBLER_HH
#define DUNE_BLOCKSTRUCTURED_LOCALASSEMBLER_HH

#include <dune/pdelab/gridoperator/default/localassembler.hh>
#include <dune/pdelab/gridfunctionspace/blockstructured/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/blockstructured/lfsindexcache.hh>

namespace Dune{
  namespace Blockstructured{

    template<typename GO, typename LOP, bool nonoverlapping_mode = false>
    class BlockstructuredLocalAssembler :
        public Dune::PDELab::DefaultLocalAssembler<GO, LOP, nonoverlapping_mode>{
    public:

      using Base =  Dune::PDELab::DefaultLocalAssembler<GO, LOP, nonoverlapping_mode>;

      using CU = typename Base::CU;
      using CV = typename Base::CV;

      using LFSU = BlockstructuredLocalFunctionSpace<typename Base::GFSU>;
      using LFSV = BlockstructuredLocalFunctionSpace<typename Base::GFSV>;

      using LFSUCache = LFSIndexCache<LFSU, CU>;
      using LFSVCache = LFSIndexCache<LFSV, CV>;

      using LocalResidualAssemblerEngine = Dune::PDELab::DefaultLocalResidualAssemblerEngine<BlockstructuredLocalAssembler>;
      using LocalJacobianApplyAssemblerEngine = Dune::PDELab::DefaultLocalJacobianApplyAssemblerEngine<BlockstructuredLocalAssembler>;

      //! Constructor with empty constraints
      BlockstructuredLocalAssembler (LOP & lop, shared_ptr<typename GO::BorderDOFExchanger> border_dof_exchanger)
        : Base(lop, border_dof_exchanger), residual_engine(*this), jacobian_apply_engine(*this)
      {}

      //! Constructor for non trivial constraints
      BlockstructuredLocalAssembler (LOP & lop, const typename Base::CU& cu, const typename Base::CV& cv,
                                     shared_ptr<typename GO::BorderDOFExchanger> border_dof_exchanger)
        : Base(lop, cu, cv, border_dof_exchanger), residual_engine(*this), jacobian_apply_engine(*this)
      {}

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalResidualAssemblerEngine & localResidualAssemblerEngine
          (typename Base::Traits::Residual & r, const typename Base::Traits::Solution & x)
      {
        residual_engine.setResidual(r);
        residual_engine.setSolution(x);
        return residual_engine;
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalJacobianApplyAssemblerEngine & localJacobianApplyAssemblerEngine
          (typename Base::Traits::Residual & r, const typename Base::Traits::Solution & x)
      {
        jacobian_apply_engine.setResidual(r);
        jacobian_apply_engine.setSolution(x);
        return jacobian_apply_engine;
      }

    private:
      LocalResidualAssemblerEngine residual_engine;
      LocalJacobianApplyAssemblerEngine jacobian_apply_engine;
    };
  }
}


#endif //DUNE_PDELAB_LOCALASSEMBLER_HH
