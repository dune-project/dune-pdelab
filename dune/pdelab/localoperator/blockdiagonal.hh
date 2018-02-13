// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_BLOCKDIAGONAL_HH
#define DUNE_PDELAB_LOCALOPERATOR_BLOCKDIAGONAL_HH

#include <memory>
#include <utility>

#include <dune/typetree/traversal.hh>
#include <dune/typetree/childextraction.hh>

#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/pattern.hh>

namespace Dune {
  namespace PDELab {

#ifndef DOXYGEN

    namespace impl {

      // TypeTree visitor for applying an operation to a pair of ansatz and test local function space
      template<typename LFSV, typename Operation>
      struct ApplyBlockOperation
        : public TypeTree::TreeVisitor
        , public TypeTree::StaticTraversal
      {

        template<typename LeafLFSU, typename TreePath>
        void leaf(LeafLFSU& leaf_lfsu, TreePath tree_path)
        {
          auto& leaf_lfsv = child(_lfsv,tree_path);
          _operation(tree_path,leaf_lfsu,leaf_lfsv);
        }

        ApplyBlockOperation(const LFSV& lfsv, Operation operation)
          : _lfsv(lfsv)
          , _operation(operation)
        {}

        const LFSV& _lfsv;
        Operation _operation;

      };

      template<typename LFSV, typename Operation>
      auto applyBlockOperation(const LFSV& lfsv, Operation operation)
      {
        return ApplyBlockOperation<LFSV,Operation>(lfsv,operation);
      }

    } // namespace impl

#endif // DOXYGEN

    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    //! Block diagonal extension of scalar local operator.
    /**
     * This adapter class takes an existing local operator that only has volume methods and extends
     * it in a block diagonal fashion to trees of function spaces.
     *
     * The scalar operator is stored internally as a shared_ptr, so you can either construct the
     * adapter from a shared_ptr to a scalar space, or you can pass it constructor arguments for
     * the scalar space and let the adapter construct the scalar space from them.
     *
     * The wrapper also implements the instationary local operator interface by forwarding all calls
     * to the scalar operator.
     *
     * If the operator is used to create a matrix, the sparsity pattern will contain all
     * off-diagonal entries for each grid cell.
     */
    template<typename ScalarLOP>
    class BlockDiagonalLocalOperatorFullCoupling
      : public FullVolumePattern
      , public LocalOperatorDefaultFlags
    {

    public:

      using RealType = typename ScalarLOP::RealType;

      static constexpr bool doPatternVolume = true;
      static constexpr bool doAlphaVolume = true;

      //! Constructs the adapter by wrapping an existing shared_ptr to the scalar operator.
      BlockDiagonalLocalOperatorFullCoupling(const std::shared_ptr<ScalarLOP>& scalar_lop)
        : _scalar_lop(scalar_lop)
      {}

      //! Constructs the adapter by wrapping an existing shared_ptr to the scalar operator.
      BlockDiagonalLocalOperatorFullCoupling(std::shared_ptr<ScalarLOP>& scalar_lop)
        : _scalar_lop(scalar_lop)
      {}

      //! Constructs the adapter and creates a scalar operator with the given arguments.
      template<typename... ScalarOperatorArgs>
      BlockDiagonalLocalOperatorFullCoupling(ScalarOperatorArgs&&... scalarOperatorArgs)
        : _scalar_lop(std::make_shared<ScalarLOP>(std::forward<ScalarOperatorArgs>(scalarOperatorArgs)...))
      {}

      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        auto visitor = impl::applyBlockOperation(
          lfsv,
          [&](auto tree_path, auto& lfsu, auto& lfsv)
          {
            _scalar_lop->alpha_volume(eg,lfsu,x,lfsv,r);
          });
        TypeTree::applyToTree(lfsu,visitor);
      }

      template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
      void jacobian_apply_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y) const
      {
        auto visitor = impl::applyBlockOperation(
          lfsv,
          [&](auto tree_path, auto& lfsu, auto& lfsv)
          {
            _scalar_lop->jacobian_apply_volume(eg,lfsu,x,lfsv,y);
          });
        TypeTree::applyToTree(lfsu,visitor);
      }

      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, M& mat) const
      {
        auto visitor = impl::applyBlockOperation(
          lfsv,
          [&](auto tree_path, auto& lfsu, auto& lfsv)
          {
            _scalar_lop->jacobian_volume(eg,lfsu,x,lfsv,mat);
          });
        TypeTree::applyToTree(lfsu,visitor);
      }

      void setTime(RealType t)
      {
        _scalar_lop->setTime(t);
      }

      RealType getTime() const
      {
        return _scalar_lop->getTime();
      }

      void preStep(RealType time, RealType dt, int stages)
      {
        _scalar_lop->preStep(time,dt,stages);
      }

      void postStep()
      {
        _scalar_lop->postStep();
      }

      void preStage(RealType time, int r)
      {
        _scalar_lop->preStage(time,r);
      }

      int getStage() const
      {
        return _scalar_lop->getStage();
      }

      void postStage()
      {
        _scalar_lop->postStage();
      }

      RealType suggestTimestep(RealType dt) const
      {
        return _scalar_lop->suggestTimeStep(dt);
      }

    private:

      std::shared_ptr<ScalarLOP> _scalar_lop;

    };

    //! \} group LocalOperator

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_BLOCKDIAGONAL_HH
