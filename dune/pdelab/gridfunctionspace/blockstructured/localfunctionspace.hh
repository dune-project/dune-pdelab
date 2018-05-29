//
// Created by marckoch on 09/05/18.
//

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_LOCALFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_LOCALFUNCTIONSPACE_HH

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/blockstructured/visitors.hh>

namespace Dune{
  namespace Blockstructured{

    template<typename GFS>
    struct gfs_to_blockstructured_lfs {

      //! The MultiIndex type that will be used in the resulting LocalFunctionSpace tree.
      //typedef Dune::PDELab::MultiIndex<std::size_t,TypeTree::TreeInfo<GFS>::depth> MultiIndex;
      typedef typename Dune::PDELab::build_dof_index_type<GFS>::type DOFIndex;
    };

    template<typename NodeType>
    void bind(NodeType &node, const typename NodeType::Traits::Element &e) {
      typedef typename NodeType::Traits::Element Element;
      // compute sizes
      ComputeSizeVisitor<Element> csv(e);
      TypeTree::applyToTree(node, csv);


      // initialize iterators and fill indices
      FillIndicesVisitor<Element> fiv(e);
      TypeTree::applyToTree(node, fiv);
    }

    // local function space for a power grid function space
    template<typename GFS, typename DOFIndex, typename ChildLFS, std::size_t k>
    class PowerLocalFunctionSpaceNode :
        public Dune::PDELab::PowerLocalFunctionSpaceNode<GFS, DOFIndex, ChildLFS, k>
    {
      using Base = Dune::PDELab::PowerLocalFunctionSpaceNode<GFS, DOFIndex, ChildLFS, k>;

    public:
      using Traits = typename Base::Traits;

      using Base::Base;

      //! \brief bind local function space to entity
      void bind (const typename Traits::Element& e)
      {
        // call method on base class, this avoid the barton neckman trick
        Dune::Blockstructured::bind(*this,e);
      }
    };


    // transformation template, we need a custom template in order to inject the DOFIndex type into the LocalFunctionSpace
    template<typename SourceNode, typename Transformation>
    struct power_gfs_to_lfs_template
    {
      template<typename TC>
      struct result
      {
        typedef PowerLocalFunctionSpaceNode<SourceNode,typename Transformation::DOFIndex,TC,TypeTree::StaticDegree<SourceNode>::value> type;
      };
    };

    // register PowerGFS -> LocalFunctionSpace transformation
    template<typename PowerGridFunctionSpace, typename Params>
    TypeTree::TemplatizedGenericPowerNodeTransformation<
        PowerGridFunctionSpace,
        gfs_to_blockstructured_lfs<Params>,
        power_gfs_to_lfs_template<PowerGridFunctionSpace,gfs_to_blockstructured_lfs<Params> >::template result
    >
    registerNodeTransformation(PowerGridFunctionSpace* pgfs, gfs_to_blockstructured_lfs<Params>* t,
                               Dune::PDELab::PowerGridFunctionSpaceTag* tag);


    // local function space for a power grid function space
    template<typename GFS, typename DOFIndex, typename... Children>
    class CompositeLocalFunctionSpaceNode
        : public Dune::PDELab::CompositeLocalFunctionSpaceNode<GFS, DOFIndex,Children...>
    {
      using Base = Dune::PDELab::CompositeLocalFunctionSpaceNode<GFS, DOFIndex,Children...>;

    public:
      using Traits = typename Base::Traits;

      using Base::Base;

      //! \brief bind local function space to entity
      void bind (const typename Traits::Element& e)
      {
        // call method on base class, this avoid the barton neckman trick
        Dune::Blockstructured::bind(*this,e);
      }
    };

    // transformation template, we need a custom template in order to inject the MultiIndex type into the LocalFunctionSpace
    template<typename SourceNode, typename Transformation>
    struct composite_gfs_to_lfs_template
    {
      template<typename... TC>
      struct result
      {
        typedef CompositeLocalFunctionSpaceNode<SourceNode,typename Transformation::DOFIndex,TC...> type;
      };
    };

    // register CompositeGFS -> LocalFunctionSpace transformation (variadic version)
    template<typename CompositeGridFunctionSpace, typename Params>
    TypeTree::TemplatizedGenericCompositeNodeTransformation<
        CompositeGridFunctionSpace,
        gfs_to_blockstructured_lfs<Params>,
        composite_gfs_to_lfs_template<CompositeGridFunctionSpace,gfs_to_blockstructured_lfs<Params> >::template result
    >
    registerNodeTransformation(CompositeGridFunctionSpace* cgfs, gfs_to_blockstructured_lfs<Params>* t,
                               Dune::PDELab::CompositeGridFunctionSpaceTag* tag);


    //! single component local function space
    template<typename GFS, typename DOFIndex>
    class LeafLocalFunctionSpaceNode
        : public Dune::PDELab::LeafLocalFunctionSpaceNode<GFS, DOFIndex>
    {
      using Base = Dune::PDELab::LeafLocalFunctionSpaceNode<GFS, DOFIndex>;

    public:
      using Traits = typename Base::Traits;

      typedef FiniteElementInterfaceSwitch<
          typename Traits::FiniteElementType
      > FESwitch;

    public:
      using Base::Base;

      //! Calculates the multiindices associated with the given entity.
      template<typename Entity, typename DOFIndexIterator>
      void dofIndices(const Entity& e, DOFIndexIterator /*it*/, DOFIndexIterator /*endit*/, std::integral_constant<bool,false>)
      {
        using EntitySet = typename GFS::Traits::EntitySet;
        auto es = this->gridFunctionSpace().entitySet();

        auto refEl = Dune::ReferenceElements<double,EntitySet::dimension>::general(this->pfe->type());

        auto& subentityWiseDOFs = *this->_subentityWiseDOFs_ptr;

        subentityWiseDOFs[this->offsetLeafs].clear();
        for (int c = 0; c < refEl.dimension + 1; ++c) {
          for (int s = 0; s < refEl.size(c); ++s) {
            // get geometry type of subentity
            auto gt = refEl.type(s, c);

            // evaluate consecutive index of subentity
            auto index = es.indexSet().subIndex(e, s, c);
            using DOFIndexAccessor = typename GFS::Ordering::Traits::DOFIndexAccessor;
            DOFIndexAccessor::store(subentityWiseDOFs[this->offsetLeafs].index(s, c), gt, index, 0);
          }
        }
      }

      void bind (const typename Traits::Element& e)
      {
        Dune::Blockstructured::bind(*this,e);
      }
    };

    // Register LeafGFS -> LocalFunctionSpace transformation
    template<typename GridFunctionSpace, typename Params>
    TypeTree::GenericLeafNodeTransformation<
        GridFunctionSpace,
        gfs_to_blockstructured_lfs<Params>,
        LeafLocalFunctionSpaceNode<GridFunctionSpace,typename gfs_to_blockstructured_lfs<Params>::DOFIndex>
    >
    registerNodeTransformation(GridFunctionSpace* gfs, gfs_to_blockstructured_lfs<Params>* t,
                               Dune::PDELab::LeafGridFunctionSpaceTag* tag);


    template <typename GFS, typename TAG=Dune::PDELab::AnySpaceTag>
    class LocalFunctionSpace:
        public TypeTree::TransformTree<GFS,gfs_to_blockstructured_lfs<GFS> >::Type
    {
      static_assert(std::is_same<TAG, Dune::PDELab::AnySpaceTag>::value, "Use this LFS only with AnySpaceTag");

      typedef typename TypeTree::TransformTree<GFS,gfs_to_blockstructured_lfs<GFS> >::Type BaseT;
      typedef typename BaseT::Traits::IndexContainer::size_type I;
      typedef typename BaseT::Traits::IndexContainer::size_type LocalIndex;

    public:

      LocalFunctionSpace(const GFS & gfs)
          : BaseT(TypeTree::TransformTree<GFS,gfs_to_blockstructured_lfs<GFS> >::transform(gfs))
          , maxLocalSize(gfs.ordering().maxLocalSize())
      {
        this->_dof_indices = &(this->_dof_index_storage);
        this->_subentityWiseDOFs_ptr = &(this->_subentityWiseDOFs);
        this->setup();
      }

      LocalFunctionSpace(std::shared_ptr<const GFS> pgfs)
      : BaseT(*TypeTree::TransformTree<GFS,gfs_to_blockstructured_lfs<GFS> >::transform_storage(pgfs))
      , maxLocalSize(pgfs->ordering().maxLocalSize())
      {
        this->_dof_indices = &(this->_dof_index_storage);
        this->_subentityWiseDOFs_ptr = &(this->_subentityWiseDOFs);
        this->setup();
      }

      LocalFunctionSpace(const LocalFunctionSpace & lfs)
          : BaseT(lfs)
          , maxLocalSize(lfs.maxSize())
      {
        // We need to reset the DOFIndex storage pointers in the new LFS tree,
        // as they are still pointing to the _dof_index_storage of the
        // old tree.
        this->_dof_indices = &(this->_dof_index_storage);
        this->_subentityWiseDOFs_ptr = &(this->_subentityWiseDOFs);
        this->setup();
      }

      LocalIndex localIndex (LocalIndex index) const
      {
        return LocalIndex(BaseT::localIndex(index));
      }

      //! \brief get maximum possible size (which is maxLocalSize from grid function space)
      typename BaseT::Traits::IndexContainer::size_type maxSize () const
      {
        // _dof_indices is always as large as the max local size of the root GFS
        return maxLocalSize;
      }

      //! \brief get size of an appropriate local vector object
      /**
         this is the number of dofs of the complete local function
         space tree, i.e. the size() of the root node. The local
         vector objects must always have this size and the localIndex
         method maps into the range [0,localVectorSize()[
       */
      typename BaseT::Traits::IndexContainer::size_type localVectorSize () const
      {
        return maxLocalSize;
      }

      void setup()
      {
        this->_subentityWiseDOFs_ptr->resize(Dune::TypeTree::TreeInfo<GFS>::leafCount);
        TypeTree::applyToTree(*this, PropagateGlobalStorageVisitor<>());
        BaseT::setup(*this);
      }

    private:

      typename BaseT::Traits::IndexContainer::size_type maxLocalSize;
    };
  }
}

#endif //DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_LOCALFUNCTIONSPACE_HH
