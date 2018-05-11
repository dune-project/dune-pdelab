//
// Created by marckoch on 09/05/18.
//

#ifndef DUNE_PDELAB_LOCALFUNCTIONSPACE_HH
#define DUNE_PDELAB_LOCALFUNCTIONSPACE_HH

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

namespace Dune{
  namespace Blockstructured{

    //! single component local function space
    template<typename GFS, typename DOFIndex>
    class BlockstructuredLeafLocalFunctionSpaceNode
        : public Dune::PDELab::LeafLocalFunctionSpaceNode<GFS, DOFIndex>
    {

      using Base = Dune::PDELab::LeafLocalFunctionSpaceNode<GFS, DOFIndex>;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ClearSizeVisitor;

      template<typename,bool>
      friend struct ComputeSizeVisitor;

      template<typename,bool>
      friend struct FillIndicesVisitor;

    public:
      using Traits = typename Base::Traits;

      typedef FiniteElementInterfaceSwitch<
          typename Traits::FiniteElementType
      > FESwitch;

    public:

      //! \brief initialize with grid function space
      template<typename Transformation>
      BlockstructuredLeafLocalFunctionSpaceNode (std::shared_ptr<const GFS> gfs, const Transformation& t)
          : Base(gfs, t)
      {
      }

      template<typename Transformation>
      BlockstructuredLeafLocalFunctionSpaceNode (const GFS& gfs, const Transformation& t)
          : Base(stackobject_to_shared_ptr(gfs), t)
      {
      }

      //! Calculates the multiindices associated with the given entity.
      template<typename Entity, typename DOFIndexIterator, bool fast>
      void dofIndices(const Entity& e, DOFIndexIterator, DOFIndexIterator, std::integral_constant<bool,fast>)
      {            // get layout of entity
        const typename FESwitch::Coefficients &coeffs =
            FESwitch::coefficients(this->finiteElement());

        using EntitySet = typename GFS::Traits::EntitySet;
        auto es = this->gridFunctionSpace().entitySet();

        auto refEl = Dune::ReferenceElements<double,EntitySet::dimension>::general(this->pfe->type());

        for (int c = 0; c < refEl.dimension + 1; ++c) {
          _dof_indices_sc[c].resize(refEl.size(c));
          for (int s = 0; s < refEl.size(c); ++s) {
            // get geometry type of subentity
            auto gt = refEl.type(s, c);

            // evaluate consecutive index of subentity
            auto index = es.indexSet().subIndex(e, s, c);
            GFS::Ordering::Traits::DOFIndexAccessor::store(_dof_indices_sc[c][s], gt, index, 0);
          }
        }
      }

      template<bool fast = false>
      void bind (const typename Traits::Element& e, std::integral_constant<bool,fast> fast_ = std::integral_constant<bool,fast>{})
      {
        Dune::PDELab::GridViewLocalFunctionSpaceBaseNode<GFS,DOFIndex>::bind(*this,e,fast_);
      }

      std::array<std::vector<typename Base::Traits::DOFIndex>,3> _dof_indices_sc;
    };

    template<typename GFS>
    struct gfs_to_blockstructured_lfs {

      //! The MultiIndex type that will be used in the resulting LocalFunctionSpace tree.
      //typedef Dune::PDELab::MultiIndex<std::size_t,TypeTree::TreeInfo<GFS>::depth> MultiIndex;
      typedef typename Dune::PDELab::build_dof_index_type<GFS>::type DOFIndex;

    };

    // Register LeafGFS -> LocalFunctionSpace transformation
    template<typename GridFunctionSpace, typename Params>
    TypeTree::GenericLeafNodeTransformation<
        GridFunctionSpace,
        gfs_to_blockstructured_lfs<Params>,
        BlockstructuredLeafLocalFunctionSpaceNode<GridFunctionSpace,typename gfs_to_blockstructured_lfs<Params>::DOFIndex>
    >
    registerNodeTransformation(GridFunctionSpace* gfs, gfs_to_blockstructured_lfs<Params>* t,
                               Dune::PDELab::LeafGridFunctionSpaceTag* tag);



    template <typename GFS, typename TAG=Dune::PDELab::AnySpaceTag>
    class BlockstructuredLocalFunctionSpace :
        public TypeTree::TransformTree<GFS,gfs_to_blockstructured_lfs<GFS> >::Type
    {
      typedef typename TypeTree::TransformTree<GFS,gfs_to_blockstructured_lfs<GFS> >::Type BaseT;
      typedef typename BaseT::Traits::IndexContainer::size_type I;
      typedef typename BaseT::Traits::IndexContainer::size_type LocalIndex;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ClearSizeVisitor;

      template<typename,bool>
      friend struct ComputeSizeVisitor;

      template<typename,bool>
      friend struct FillIndicesVisitor;

    public:

      BlockstructuredLocalFunctionSpace(const GFS & gfs)
          : BaseT(TypeTree::TransformTree<GFS,gfs_to_blockstructured_lfs<GFS> >::transform(gfs))
      {
        this->_dof_indices = &(this->_dof_index_storage);
        this->setup(*this);
      }

      BlockstructuredLocalFunctionSpace(std::shared_ptr<const GFS> pgfs)
      : BaseT(*TypeTree::TransformTree<GFS,gfs_to_blockstructured_lfs<GFS> >::transform_storage(pgfs))
      {
        this->_dof_indices = &(this->_dof_index_storage);
        this->setup(*this);
      }

      BlockstructuredLocalFunctionSpace(const BlockstructuredLocalFunctionSpace & lfs)
          : BaseT(lfs)
      {
        // We need to reset the DOFIndex storage pointers in the new LFS tree,
        // as they are still pointing to the _dof_index_storage of the
        // old tree.
        this->_dof_indices = &(this->_dof_index_storage);
        this->setup(*this);
      }

      LocalIndex localIndex (LocalIndex index) const
      {
        return LocalIndex(BaseT::localIndex(index));
      }

    private:
      // we don't support getChild yet, so let's hide it!
      template<int i>
      void getChild () const;
      template<int i>
      void child () const;
    };
  }
}

#endif //DUNE_PDELAB_LOCALFUNCTIONSPACE_HH
