// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDVIEWORDERING_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDVIEWORDERING_HH

#include <dune/pdelab/common/typetree.hh>
#include <dune/pdelab/gridfunctionspace/orderingutility.hh>
#include <dune/pdelab/gridfunctionspace/localorderingdynamicbase.hh>
#include <dune/pdelab/gridfunctionspace/orderingdynamicbase.hh>
#include <dune/pdelab/gridfunctionspace/directleaflocalordering.hh>

namespace Dune {
  namespace PDELab {



    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //! Dummy ordering for leaf gridfunctionspaces
    template<typename LocalOrdering>
    class LeafGridViewOrdering
      : public TypeTree::VariadicCompositeNode<LocalOrdering>
      , public VirtualOrderingBase<typename LocalOrdering::Traits::DOFIndex,typename LocalOrdering::Traits::ContainerIndex>
      , public OrderingBase<typename LocalOrdering::Traits::DOFIndex,typename LocalOrdering::Traits::ContainerIndex>
    {
    public:
      typedef typename LocalOrdering::Traits Traits;

      static const bool has_dynamic_ordering_children = false;

      static const bool consume_tree_index = false;

    private:

      typedef typename Traits::GridView GV;

      typedef TypeTree::VariadicCompositeNode<LocalOrdering> NodeT;

      typedef OrderingBase<typename LocalOrdering::Traits::DOFIndex,typename LocalOrdering::Traits::ContainerIndex> BaseT;

    public:

      LocalOrdering& localOrdering()
      {
        return this->template child<0>();
      }

      const LocalOrdering& localOrdering() const
      {
        return this->template child<0>();
      }


      LeafGridViewOrdering(const typename NodeT::NodeStorage& localOrdering, bool container_blocked)
        : NodeT(localOrdering)
        , BaseT(*this,container_blocked,this)
        , _gv(this->template child<0>().gridView())
      {}

      virtual void map_index_dynamic(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {
        map_index(di,ci);
      }

      typename Traits::ContainerIndex map_index(const typename Traits::DOFIndex& di) const
      {
        typename Traits::ContainerIndex ci;
        map_index(di.view(),ci);
        return ci;
      }

      void map_index(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {
        const typename Traits::SizeType geometry_type_index = di.entityIndex()[0];
        const typename Traits::SizeType entity_index = di.entityIndex()[1];
        assert (di.treeIndex().size() == 1);
        ci.push_back(di.treeIndex().back());
        if (_container_blocked)
          {
            ci.push_back(localOrdering()._gt_entity_offsets[geometry_type_index] + entity_index);
          }
        else if (localOrdering()._fixed_size)
          {
            ci.back() += _gt_dof_offsets[geometry_type_index] + entity_index * localOrdering()._gt_dof_offsets[geometry_type_index];
          }
        else
          {
            ci.back() += localOrdering()._entity_dof_offsets[localOrdering()._gt_entity_offsets[geometry_type_index] + entity_index];
          }
      }


      template<typename ItIn, typename ItOut>
      void map_indices(const ItIn begin, const ItIn end, ItOut out) const
      {
        typedef typename Traits::SizeType size_type;
        if (_container_blocked)
          {
            for (ItIn in = begin; in != end; ++in, ++out)
              {
                assert(in->treeIndex().size() == 1);
                out->push_back(in->treeIndex().back());
                out->push_back(localOrdering()._gt_entity_offsets[in->entityIndex()[0]] + in->entityIndex()[1]);
              }
          }
        else if (localOrdering()._fixed_size)
          {
            for (ItIn in = begin; in != end; ++in, ++out)
              {
                assert(in->treeIndex().size() == 1);
                out->push_back(in->treeIndex().back());
                const size_type geometry_type_index = in->entityIndex()[0];
                const size_type entity_index = in->entityIndex()[1];
                out->back() += _gt_dof_offsets[geometry_type_index] + entity_index * localOrdering()._gt_dof_offsets[geometry_type_index];
              }
          }
        else
          {
            for (ItIn in = begin; in != end; ++in, ++out)
              {
                assert(in->treeIndex().size() == 1);
                out->push_back(in->treeIndex().back());
                const size_type geometry_type_index = in->entityIndex()[0];
                const size_type entity_index = in->entityIndex()[1];
                out->back() += localOrdering()._entity_dof_offsets[localOrdering()._gt_entity_offsets[geometry_type_index] + entity_index];
              }
          }
      }


      void update()
      {
        _recursive_update();
      }

      //! update internal data structures
      void _recursive_update()
      {
        LocalOrdering& lo = localOrdering();
        lo.update_a_priori_fixed_size();
        _fixed_size = lo._fixed_size;

        const std::size_t dim = GV::dimension;

        typedef typename Traits::SizeType size_type;
        typedef std::vector<GeometryType> GTVector;
        GTVector geom_types;

        for (size_type cc = 0; cc <= dim; ++cc)
          {
            const GTVector& per_codim_geom_types = _gv.indexSet().geomTypes(cc);
            std::copy(per_codim_geom_types.begin(),per_codim_geom_types.end(),std::back_inserter(geom_types));
          }

        if (lo._fixed_size)
          {
            lo.update_fixed_size(geom_types);

            _gt_dof_offsets.assign(GlobalGeometryTypeIndex::size(dim) + 1,0);

            const GTVector::const_iterator end_it = geom_types.end();
            for (GTVector::const_iterator it = geom_types.begin(); it != end_it; ++it)
              {
                const size_type gt_index = GlobalGeometryTypeIndex::index(*it);
                _gt_dof_offsets[gt_index + 1] = lo.size(gt_index,0) * _gv.indexSet().size(*it);
              }
            std::partial_sum(_gt_dof_offsets.begin(),_gt_dof_offsets.end(),_gt_dof_offsets.begin());
            _block_count = _size = _gt_dof_offsets.back();
            _max_local_size = lo.maxLocalSize();
          }
        else
          {
            assert(!"Not implemented yet!");
            // FIXME
          }
      }

      //! dofs are blocked per entity/intersection on the leafs
      bool blocked() const { return true; }

      //! \brief whether all entites of the same geometry type/all
      //!        intersections have the same number of dofs
      /**
       * On the leaf this is realized by iterating over the grid during update
       * an checking.
       *
       * \note Even if fixedSize()==true the number of dofs may still vary
       *       between entities od different geometry type or between entities
       *       and intersections.
       */
      bool fixedSize() const { return _fixed_size; }

      //! \brief maximum number of dofs attached to any given element and all
      //!        of its subentities and intersections
      /**
       * This is generally not an exact maximum and may be bigger than the
       * actual maximum.  There is however one special case: it is guaranteed
       * to be the exact maximum for fixedSize()==true.
       */

#if 0
      //! \brief number of indices attached to a given entity (of arbitrary
      //!        codimension)
      /**
       * \note If the grid does not support a given entity type, it may still
       *       be possible to get this information using entitySize(const
       *       Element &e, std::size_t codim, std::size_t subentity).
       */
      template<class Entity>
      typename Traits::SizeType entitySize(const Entity &e) const { return gfs.entitySize(e); }
      //! number of indices attached to a given subentity of an element
      /**
       * This method determines the number of indices attached to a subentity
       * of the given codim 0 entity.  If the grid (and the ordering) directly
       * supports entities of the given codimension, this is equivalent to
       * calling entitySize((*e.subEntity<codim>(subentity)).
       */
      template<class Element>
      typename Traits::SizeType entitySize(const Element &e, std::size_t codim,
                          std::size_t subentity) const
      { return gfs.entitySize(e, codim, subentity); }
      //! number of indices attached to a given intersection
      template<class Intersection>
      typename Traits::SizeType intersectionSize(const Intersection &i) const
      { return gfs.intersectionSize(i); }

      //! \brief offset of the block of dofs attached to a given entity (of
      //!        arbitrary codimension)
      /**
       * \note If the grid does not support a given entity type, it may still
       *       be possible to get this information using entityOffset(const
       *       Element &e, std::size_t codim, std::size_t subentity).
       *
       * \throw NotImplemented        If this EntityType is not supported by
       *                              the ordering.
       * \throw InvalidStateException If blocked()==false.
       */
      template<class Entity>
      typename Traits::SizeType entityOffset(const Entity &e) const
      { return gfs.entityOffset(e); }
      //! \brief offset of the blocks of dofs attached to a given subentity of
      //!        an element
      /**
       * This method determines the starting offset of the block of dofs
       * attached to a subentity of the given codim 0 entity.  If the grid
       * (and the ordering) directly support entities of the given
       * codimension, this is equivalent to calling
       * entityOffset(*e.subEntity<codim>(subentity)).
       */
      template<class Element>
      typename Traits::SizeType entityOffset(const Element &e, std::size_t codim,
                            std::size_t subentity) const
      { return gfs.entityOffset(e, codim, subentity); }
      //! offset of the block of dofs attached to a given intersection
      template<class Intersection>
      typename Traits::SizeType intersectionOffset(const Intersection &i) const
      { return gfs.intersectionOffset(i); }

#endif // 0

    private:

      using BaseT::_max_local_size;
      using BaseT::_size;
      using BaseT::_block_count;
      using BaseT::_container_blocked;
      using BaseT::_fixed_size;

      typename Traits::GridView _gv;
      std::vector<typename Traits::SizeType> _gt_dof_offsets;

    };

    template<typename GFS, typename Transformation, typename SizeTag>
    struct leaf_gfs_to_ordering_descriptor;

    struct GridFunctionGeneralMapper;

    template<typename GFS, typename Transformation>
    struct leaf_gfs_to_ordering_descriptor<GFS,Transformation,GridFunctionGeneralMapper>
    {

      static const bool recursive = false;

      typedef DirectLeafLocalOrdering<GFS,
                                      typename Transformation::DOFIndex,
                                      typename Transformation::ContainerIndex
                                      > LocalOrdering;

      typedef LeafGridViewOrdering<LocalOrdering> GridViewOrdering;

      typedef GridViewOrdering transformed_type;
      typedef shared_ptr<transformed_type> transformed_storage_type;

      static transformed_type transform(const GFS& gfs, const Transformation& t)
      {
        return transformed_type(make_tuple(make_shared<LocalOrdering>(gfs)),gfs.backend().blocked());
      }

      static transformed_storage_type transform_storage(shared_ptr<const GFS> gfs, const Transformation& t)
      {
        return make_shared<transformed_type>(make_tuple(make_shared<LocalOrdering>(*gfs)),gfs->backend().blocked());
      }

    };



    template<typename GridFunctionSpace, typename Params>
    leaf_gfs_to_ordering_descriptor<
      GridFunctionSpace,
      gfs_to_ordering<Params>,
      typename GridFunctionSpace::SizeTag
      >
    lookupNodeTransformation(GridFunctionSpace* gfs, gfs_to_ordering<Params>* t, LeafGridFunctionSpaceTag tag);


    struct collect_a_priori_fixed_size
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp)
      {
        node._fixed_size = node.gridFunctionSpace().fixedSize();
        any = any || node._fixed_size;
      }

      template<typename Node, typename TreePath>
      void pre(Node& node, TreePath tp) const
      {
        node._fixed_size = true;
      }

      template<typename Node, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(Node& node, const Child& child, TreePath tp, ChildIndex childIndex) const
      {
        node._fixed_size = node._fixed_size && child._fixed_size;
      }

      collect_a_priori_fixed_size()
        : any(false)
      {}

      bool any;

    };


    template<typename GV>
    struct update_fixed_size
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      typedef std::vector<Dune::GeometryType> GTVector;

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp) const
      {
        if (node._fixed_size)
          {
            typedef typename Node::Traits::SizeType size_type;
            const size_type dim = GV::dimension;
            node._codim_used.assign(dim,false);
            node._gt_used.assign(GlobalGeometryTypeIndex::size(dim),false);
            node._gt_dof_offsets.assign(GlobalGeometryTypeIndex::size(dim),0);
            for (GTVector::const_iterator it = geom_types.begin(); it != geom_types.end(); ++it)
              {
                size_type size = node.gridFunctionSpace().size(*it);
                node._gt_dof_offsets[GlobalGeometryTypeIndex::index(*it)] = size;
                node._gt_used[GlobalGeometryTypeIndex::index(*it)] = size > 0;
                node._codim_used[dim - it->dim()] = node._codim_used[dim - it->dim()] || (size > 0);
              }
            node._max_local_size = node.gridFunctionSpace().fixedMaxLocalSize();
          }
      }

      template<typename Node, typename TreePath>
      void pre(Node& node, TreePath tp) const
      {
        if (node._fixed_size)
          {
            typedef typename Node::Traits::SizeType size_type;
            const size_type dim = GV::dimension;
            node._codim_used.assign(dim,false);
            node._gt_used.assign(Dune::GlobalGeometryTypeIndex::size(dim),false);
            node._gt_dof_offsets.assign(Dune::GlobalGeometryTypeIndex::size(dim) * Node::CHILDREN,0);
            node._max_local_size = 0;
          }
      }

      template<typename Node, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(Node& node, const Child& child, TreePath tp, ChildIndex childIndex) const
      {
        if (node._fixed_size)
          {
            std::transform(node._codim_used.begin(),
                           node._codim_used.end(),
                           child._codim_used.begin(),
                           node._codim_used.begin(),
                           std::logical_or<bool>());
            std::transform(node._gt_used.begin(),
                           node._gt_used.end(),
                           child._gt_used.begin(),
                           node._gt_used.begin(),
                           std::logical_or<bool>());

            node._max_local_size += child._max_local_size;

            typedef typename Node::Traits::SizeType size_type;

            const size_type per_gt_size = child._child_count > 0 ? child._child_count : 1;
            const size_type size_offset = child._child_count > 0 ? child._child_count - 1 : 0;

            for (size_type gt = 0; gt < Dune::GlobalGeometryTypeIndex::size(GV::dimension); ++gt)
              node._gt_dof_offsets[gt * Node::CHILDREN + childIndex] = child._gt_dof_offsets[gt * per_gt_size + size_offset];
          }
      }

      template<typename Node, typename TreePath>
      void post(Node& node, TreePath tp) const
      {
        if (node._fixed_size)
          {
            typedef typename std::vector<typename Node::Traits::SizeType>::iterator iterator;

            iterator next_gt_it = node._gt_dof_offsets.begin() + Node::CHILDREN;
            const iterator end_it = node._gt_dof_offsets.end();

            for (iterator it = node._gt_dof_offsets.begin();
                 it != end_it;
                 it += Node::CHILDREN, next_gt_it += Node::CHILDREN)
              std::partial_sum(it,next_gt_it,it);
          }
      }

      update_fixed_size(const GV gv_, const GTVector& geom_types_)
        : gv(gv_)
        , geom_types(geom_types_)
      {}

      GV gv;
      const GTVector& geom_types;

    };


    struct pre_collect_used_geometry_types
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp) const
      {
        if (!node._fixed_size)
          {
            node._codim_used.assign(dim,false);
            node._gt_used.assign(Dune::GlobalGeometryTypeIndex::size(dim),false);
            node._gt_dof_offsets.assign(Dune::GlobalGeometryTypeIndex::size(dim) * node._child_count,0);
            node._gt_entity_offsets.assign(Dune::GlobalGeometryTypeIndex::size(dim) + 1,0);
          }
      }

      template<typename Node, typename TreePath>
      void pre(Node& node, TreePath tp) const
      {
        leaf(node,tp);
      }

      pre_collect_used_geometry_types(std::size_t dimension)
        : dim(dimension)
      {}

      const std::size_t dim;

    };


    template<typename Cell>
    struct collect_used_geometry_types_from_cell
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp) const
      {
        if (!node._fixed_size)
          node.collect_used_geometry_types_from_cell(cell);
      }

      collect_used_geometry_types_from_cell(const Cell& cell_)
        : cell(cell_)
        , ref_el(Dune::GenericReferenceElements<typename Cell::ctype,Cell::dimension>::general(cell_.type()))
      {}

      const Cell& cell;
      const Dune::GenericReferenceElement<typename Cell::ctype,Cell::dimension>& ref_el;

    };


    template<typename GV>
    struct post_collect_used_geometry_types
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      typedef std::vector<Dune::GeometryType> GTVector;


      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp) const
      {
        if (!node._fixed_size)
          {
            typedef typename Node::Traits::SizeType size_type;

            for (GTVector::const_iterator it = geom_types.begin(); it != geom_types.end(); ++it)
              {
                if (node._gt_used[Dune::GlobalGeometryTypeIndex::index(*it)])
                  node._gt_entity_offsets[Dune::GlobalGeometryTypeIndex::index(*it) + 1] = gv.indexSet().size(*it);
              }

            std::partial_sum(node._gt_entity_offsets.begin(),node._gt_entity_offsets.end(),node._gt_entity_offsets.begin());
            node._entity_dof_offsets.assign(node._gt_entity_offsets.back() * node._child_count,0);
          }
      }

      template<typename Node, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(Node& node, const Child& child, TreePath tp, ChildIndex childIndex) const
      {
        if (!node._fixed_size)
          {
            std::transform(node._codim_used.begin(),
                           node._codim_used.end(),
                           child._codim_used.begin(),
                           node._codim_used.begin(),
                           std::logical_or<bool>());
            std::transform(node._gt_used.begin(),
                           node._gt_used.end(),
                           child._gt_used.begin(),
                           node._gt_used.begin(),
                           std::logical_or<bool>());
          }
      }

      template<typename Node, typename TreePath>
      void post(Node& node, TreePath tp) const
      {
        leaf(node,tp);
      }

      post_collect_used_geometry_types(const GV& gv_, const GTVector& geom_types_)
        : gv(gv_)
        , geom_types(geom_types_)
      {}

      GV gv;
      const GTVector& geom_types;

    };


    template<typename GV>
    struct extract_per_entity_sizes_from_cell
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      static const std::size_t dim = GV::dimension;
      typedef typename GV::template Codim<0>::Entity Cell;
      typedef std::size_t size_type;

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp)
      {
        if (!node._fixed_size)
          node.extract_per_entity_sizes_from_cell(*cell,gt_sizes);
      }

      extract_per_entity_sizes_from_cell(const GV& gv_)
        : gv(gv_)
        , cell(nullptr)
        , ref_el(nullptr)
        , gt_sizes(Dune::GlobalGeometryTypeIndex::size(dim),0)
      {}

      void set_cell(const Cell& cell_)
      {
        cell = &cell_;
        ref_el = &(Dune::GenericReferenceElements<typename GV::ctype,dim>::general(cell_.type()));
      }

      GV gv;
      const Cell* cell;
      const Dune::GenericReferenceElement<typename GV::ctype,dim>* ref_el;
      std::vector<size_type> gt_sizes;

    };


    template<typename GV>
    struct post_extract_per_entity_sizes
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      typedef std::vector<GeometryType> GTVector;


      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp) const
      {
        if (!node._fixed_size)
          {
            if (node._fixed_size_possible)
              {
                node._entity_dof_offsets = std::vector<typename Node::Traits::SizeType>();
                node._fixed_size = true;
              }
          }
      }

      template<typename Node, typename TreePath>
      void pre(Node& node, TreePath tp) const
      {
        if (!node._fixed_size)
          node._fixed_size_possible = true;
      }


      template<typename Node, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(Node& node, const Child& child, TreePath tp, ChildIndex childIndex) const
      {
        if (!node._fixed_size)
          node._fixed_size_possible = node._fixed_size_possible && child._fixed_size;
      }


      template<typename Node, typename TreePath>
      void post(Node& node, TreePath tp) const
      {
        if (!node._fixed_size)
          {
            typedef typename Node::Traits::SizeType size_type;
            const size_type dim = GV::dimension;

            if (node._fixed_size_possible)
              {
                for (size_type gt = 0; gt < GlobalGeometryTypeIndex::size(GV::dimension); ++gt)
                  {
                    size_type carry = 0;
                    for (size_type child_index = 0; child_index < Node::CHILDREN; ++child_index)
                      node._gt_dof_offsets[gt * Node::CHILDREN + child_index] = (carry += node.dynamic_child(child_index)._gt_dof_offsets[gt * node.dynamic_child(child_index)._child_count + node.dynamic_child(child_index)._child_count - 1]);
                  }
                node._fixed_size = true;
              }
            else
              {
                typedef typename Node::Traits::SizeType size_type;

                for (GTVector::const_iterator it = geom_types.begin(); it != geom_types.end(); ++it)
                  {
                    if (node._gt_used[GlobalGeometryTypeIndex::index(*it)])
                      node._gt_entity_offsets[GlobalGeometryTypeIndex::index(*it) + 1] = gv.indexSet().size(*it);
                  }

                std::partial_sum(node._gt_entity_offsets.begin(),node._gt_entity_offsets.end(),node._gt_entity_offsets.begin());
                node._entity_dof_offsets.assign(node._gt_entity_offsets.back() * node._child_count,0);

                for (size_type geometry_type_index = 0; geometry_type_index < GlobalGeometryTypeIndex::size(dim); ++geometry_type_index)
                  {
                    if (!node._gt_used[geometry_type_index])
                      continue;
                    for (size_type entity_index = 0; entity_index < node._gt_entity_offsets[geometry_type_index+1]; ++entity_index)
                      {
                        size_type carry = 0;
                        for (size_type child_index = 0; child_index < Node::CHILDREN; ++child_index)
                          node._entity_dof_offsets[node._gt_entity_offsets[geometry_type_index] + entity_index] = (carry += node.dynamic_child(child_index).size(geometry_type_index,entity_index));
                      }
                  }
              }
          }
      }

      post_extract_per_entity_sizes(const GV& gv_, const GTVector& geom_types_)
        : gv(gv_)
        , geom_types(geom_types_)
      {}

      GV gv;
      const GTVector& geom_types;

    };


    template<typename LocalOrdering>
    class GridViewOrdering
      : public TypeTree::VariadicCompositeNode<LocalOrdering>
      , public VirtualOrderingBase<typename LocalOrdering::Traits::DOFIndex, typename LocalOrdering::Traits::ContainerIndex>
      , public OrderingBase<typename LocalOrdering::Traits::DOFIndex, typename LocalOrdering::Traits::ContainerIndex>
    {
    public:
      typedef typename LocalOrdering::Traits Traits;

      static const bool has_dynamic_ordering_children = false;

      static const bool consume_tree_index = false;

    private:

      typedef TypeTree::VariadicCompositeNode<LocalOrdering> NodeT;
      typedef OrderingBase<typename LocalOrdering::Traits::DOFIndex, typename LocalOrdering::Traits::ContainerIndex> BaseT;

      typedef typename Traits::GridView GV;

    public:
      //! Construct ordering object
      /**
       * In general, an ordering object is not properly setup after
       * construction.  This must be done by a seperate call to update().
       * This particular ordering however can be used right away.
       */
      GridViewOrdering(const typename NodeT::NodeStorage& local_ordering, bool container_blocked)
        : NodeT(local_ordering)
        , BaseT(*this,container_blocked,this)
        , _gv(localOrdering().gridView())
      {}

      LocalOrdering& localOrdering()
      {
        return this->template child<0>();
      }

      const LocalOrdering& localOrdering() const
      {
        return this->template child<0>();
      }

      virtual void map_index_dynamic(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {
        map_index(di,ci);
      }

      typename Traits::ContainerIndex map_index(const typename Traits::DOFIndex& di) const
      {
        typename Traits::ContainerIndex ci;
        map_index(di.view(),ci);
        return ci;
      }

      void map_index(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {
        typedef typename Traits::SizeType size_type;
        const size_type geometry_type_index = di.entityIndex()[0];
        const size_type entity_index = di.entityIndex()[1];
        localOrdering().map_local_index(geometry_type_index,entity_index,di.treeIndex(),ci);
        if (_container_blocked)
          {
            ci.push_back(_gt_entity_offsets[geometry_type_index] + entity_index);
          }
        else if (_fixed_size)
          {
            ci.back() += _gt_dof_offsets[geometry_type_index] + entity_index * localOrdering().size(geometry_type_index,entity_index);
          }
        else
          {
            ci.back() += _entity_dof_offsets[_gt_entity_offsets[geometry_type_index] + entity_index];
          }
      }

      template<typename ItIn, typename ItOut>
      void map_indices(const ItIn begin, const ItIn end, ItOut out) const
      {
        typedef typename Traits::SizeType size_type;
        if (_container_blocked)
          {
            if (_fixed_size)
              for (ItIn in = begin; in != end; ++in, ++out)
                {
                  const size_type geometry_type_index = in->entityIndex()[0];
                  const size_type entity_index = in->entityIndex()[1];
                  out->push_back(_gt_dof_offsets[geometry_type_index] + entity_index);
                }
            else
              for (ItIn in = begin; in != end; ++in, ++out)
                {
                  out->push_back(_gt_entity_offsets[in->entityIndex()[0]] + in->entityIndex()[1]);
                }
          }
        else if (_fixed_size)
          {
            for (ItIn in = begin; in != end; ++in, ++out)
              {
                const size_type geometry_type_index = in->entityIndex()[0];
                const size_type entity_index = in->entityIndex()[1];
                out->back() += _gt_dof_offsets[geometry_type_index] + entity_index * localOrdering().size(geometry_type_index,entity_index);
              }
          }
        else
          {
            for (ItIn in = begin; in != end; ++in, ++out)
              {
                const size_type geometry_type_index = in->entityIndex()[0];
                const size_type entity_index = in->entityIndex()[1];
                out->back() += _entity_dof_offsets[_gt_entity_offsets[geometry_type_index] + entity_index];
              }
          }
      }

      void update()
      {

        typedef std::vector<GeometryType> GTVector;
        typedef typename Traits::SizeType size_type;
        const size_type dim = GV::dimension;
        GTVector geom_types;

        for (size_type cc = 0; cc <= dim; ++cc)
          {
            const GTVector& per_codim_geom_types = _gv.indexSet().geomTypes(cc);
            std::copy(per_codim_geom_types.begin(),per_codim_geom_types.end(),std::back_inserter(geom_types));
          }

        // Do we already know that we have fixed per-GeometryType sizes?
        collect_a_priori_fixed_size fixed_size_collector;
        TypeTree::applyToTree(localOrdering(),fixed_size_collector);
        _fixed_size = localOrdering().fixedSize();

        typedef std::vector<GeometryType> GTVector;
        const size_type gt_index_count = GlobalGeometryTypeIndex::size(GV::dimension);

        if (fixed_size_collector.any)
          {
            // collect used GeometryTypes
            TypeTree::applyToTree(localOrdering(),update_fixed_size<GV>(_gv,geom_types));

            _gt_dof_offsets.resize(gt_index_count + 1);
            _block_count = 0;

            const GTVector::const_iterator end_it = geom_types.end();
            for (GTVector::const_iterator it = geom_types.begin(); it != end_it; ++it)
              {
                const size_type gt_index = GlobalGeometryTypeIndex::index(*it);
                const size_type gt_size = localOrdering().size(gt_index,0);
                const size_type gt_entity_count = _gv.indexSet().size(*it);
                _gt_dof_offsets[gt_index + 1] = gt_size * gt_entity_count;
                _block_count += (gt_size > 0) * gt_entity_count;
              }
            std::partial_sum(_gt_dof_offsets.begin(),_gt_dof_offsets.end(),_gt_dof_offsets.begin());
            _size = _gt_dof_offsets.back();
            if (!_container_blocked)
              _block_count = _size;
          }
        else
          {
            TypeTree::applyToTree(localOrdering(),pre_collect_used_geometry_types(GV::dimension));

            typedef typename GV::template Codim<0>::Iterator CellIterator;
            typedef typename GV::template Codim<0>::Entity Cell;

            const CellIterator end_it = _gv.template end<0>();
            for (CellIterator it = _gv.template begin<0>(); it != end_it; ++it)
              {
                TypeTree::applyToTree(localOrdering(),collect_used_geometry_types_from_cell<Cell>(*it));
              }
            TypeTree::applyToTree(localOrdering(),post_collect_used_geometry_types<GV>(_gv,geom_types));
            // allocate

            //TypeTree::applyToTree(localOrdering(),pre_extract_per_entity_sizes<GV>(_gv));
            extract_per_entity_sizes_from_cell<GV> visitor(_gv);
            for (CellIterator it = _gv.template begin<0>(); it != end_it; ++it)
              {
                visitor.set_cell(*it);
                TypeTree::applyToTree(localOrdering(),visitor);
              }
            TypeTree::applyToTree(localOrdering(),post_extract_per_entity_sizes<GV>(_gv,geom_types));

            if (localOrdering().fixedSize())
              {
                _fixed_size = true;
                _gt_dof_offsets.resize(gt_index_count + 1);

                _block_count = 0;

                for (GTVector::const_iterator it = geom_types.begin(); it != geom_types.end(); ++it)
                  {
                    const size_type gt_index = GlobalGeometryTypeIndex::index(*it);
                    const size_type gt_size = localOrdering().size(gt_index,0);
                    const size_type gt_entity_count = _gv.indexSet().size(*it);
                    _gt_dof_offsets[gt_index + 1] = gt_size * gt_entity_count;
                    _block_count += (gt_size > 0) * gt_entity_count;
                  }

                std::partial_sum(_gt_dof_offsets.begin(),_gt_dof_offsets.end(),_gt_dof_offsets.begin());
                _size = _gt_dof_offsets.back();

                if (!_container_blocked)
                  _block_count = _size;

              }
            else
              {
                _gt_entity_offsets.assign(gt_index_count + 1,0);

                for (GTVector::const_iterator it = geom_types.begin(); it != geom_types.end(); ++it)
                  {
                    if (!localOrdering().contains(*it))
                      continue;
                    const size_type gt_index = GlobalGeometryTypeIndex::index(*it);
                    _gt_entity_offsets[gt_index + 1] = _gv.indexSet().size(*it);
                  }

                std::partial_sum(_gt_entity_offsets.begin(),_gt_entity_offsets.end(),_gt_entity_offsets.begin());
                _entity_dof_offsets.assign(_gt_entity_offsets.back()+1,0);
                _block_count = 0;

                size_type carry = 0;
                for (GTVector::const_iterator it = geom_types.begin(); it != geom_types.end(); ++it)
                  {
                    if (!localOrdering().contains(*it))
                      continue;
                    const size_type gt_index = Dune::GlobalGeometryTypeIndex::index(*it);
                    size_type entity_pos = _gt_entity_offsets[gt_index] + 1;
                    const size_type entity_count = _gt_entity_offsets[gt_index + 1] - entity_pos + 1;
                    for (size_type entity_index = 0; entity_index < entity_count; ++entity_index, ++entity_pos)
                      {
                        const size_type size = localOrdering().size(gt_index,entity_index);
                        _entity_dof_offsets[entity_pos] = (carry += size);
                        _block_count += (size > 0);
                      }
                  }
                _size = _entity_dof_offsets.back();

                if (!_container_blocked)
                  _block_count = _size;
              }
          }
        _max_local_size = localOrdering().maxLocalSize();
      }

    private:

      using BaseT::_container_blocked;
      using BaseT::_fixed_size;
      using BaseT::_max_local_size;
      using BaseT::_child_offsets;
      using BaseT::_size;
      using BaseT::_block_count;

      typename Traits::GridView _gv;
      std::vector<typename Traits::SizeType> _gt_dof_offsets;
      std::vector<typename Traits::SizeType> _gt_entity_offsets;
      std::vector<typename Traits::SizeType> _entity_dof_offsets;

    };


   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_LEAFORDERING_HH
