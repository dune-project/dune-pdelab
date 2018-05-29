// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDOPERATOR_BLOCKSTRUCTURED_HH
#define DUNE_PDELAB_GRIDOPERATOR_BLOCKSTRUCTURED_HH

#include <dune/pdelab/gridoperator/blockstructured/localassembler.hh>
#include <dune/pdelab/gridoperator/blockstructured/assembler.hh>
#include "gridoperator.hh"

namespace Dune{
  namespace Blockstructured{
        /**
       \brief Standard grid operator implementation

       \tparam GFSU GridFunctionSpace for ansatz functions
       \tparam GFSV GridFunctionSpace for test functions
       \tparam MB The matrix backend to be used for representation of the jacobian
       \tparam DF The domain field type of the operator
       \tparam RF The range field type of the operator
       \tparam JF The jacobian field type
       \tparam CU   Constraints maps for the individual dofs (trial space)
       \tparam CV   Constraints maps for the individual dofs (test space)

    */
    template<typename GFSU, typename GFSV, typename LOP,
             typename MB, typename DF, typename RF, typename JF,
             typename CU=Dune::PDELab::EmptyTransformation,
             typename CV=Dune::PDELab::EmptyTransformation
             >
    class BlockstructuredGridOperator
        : public Dune::PDELab::GridOperator<GFSU,GFSV,LOP,MB,DF,RF,JF,CU,CV>
    {
    public:

      using Base = Dune::PDELab::GridOperator<GFSU,GFSV,LOP,MB,DF,RF,JF,CU,CV>;

      //! The global assembler type
      using Assembler = BlockstructuredAssembler<GFSU,GFSV,CU,CV>;

      //! The global assembler type
      using LocalAssembler = BlockstructuredLocalAssembler<BlockstructuredGridOperator, LOP,
          GFSU::Traits::EntitySet::Partitions::partitionIterator() == InteriorBorder_Partition>;


      //! Constructor for non trivial constraints
      BlockstructuredGridOperator(const GFSU & gfsu_, const CU & cu_, const GFSV & gfsv_, const CV & cv_, LOP & lop_, const MB& mb_ = MB())
        : Base(gfsu_, cu_, gfsv_, cv_, lop_, mb_), global_assembler(gfsu_,gfsv_,cu_,cv_)
          , dof_exchanger(std::make_shared<typename Base::BorderDOFExchanger>(*this))
          , local_assembler(lop_, cu_, cv_,dof_exchanger)
      {}

      //! Constructor for empty constraints
      BlockstructuredGridOperator(const GFSU & gfsu_, const GFSV & gfsv_, LOP & lop_, const MB& mb_ = MB())
        : Base(gfsu_, gfsv_, lop_, mb_), global_assembler(gfsu_,gfsv_)
          , dof_exchanger(std::make_shared<typename Base::BorderDOFExchanger>(*this))
          , local_assembler(lop_, dof_exchanger)
      {}


      //! Assemble residual
      void residual(const typename Base::Domain & x, typename Base::Range & r) const
      {
        auto & residual_engine = local_assembler.localResidualAssemblerEngine(r,x);
        global_assembler.assemble(residual_engine);
      }

      //! Apply jacobian matrix without explicitly assembling it
      void jacobian_apply(const typename Base::Domain & z, typename Base::Range & r) const
      {
        auto & jacobian_apply_engine = local_assembler.localJacobianApplyAssemblerEngine(r,z);
        global_assembler.assemble(jacobian_apply_engine);
      }

    private:
      Assembler global_assembler;
      shared_ptr<typename Base::BorderDOFExchanger> dof_exchanger;

      mutable LocalAssembler local_assembler;

    };
  }
}

#endif //DUNE_PDELAB_GRIDOPERATOR_BLOCKSTRUCTURED_HH
