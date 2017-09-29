// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_GRIDVIEWORDERING_HH
#define DUNE_PDELAB_ORDERING_GRIDVIEWORDERING_HH

#include <dune/typetree/typetree.hh>

#include <dune/pdelab/ordering/utility.hh>
#include <dune/pdelab/ordering/localorderingbase.hh>
#include <dune/pdelab/ordering/orderingbase.hh>
#include <dune/pdelab/ordering/leaflocalordering.hh>
#include <dune/pdelab/ordering/lexicographicordering.hh>
#include <dune/pdelab/ordering/leafgridviewordering.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Ordering
    //! \{

    template<typename Codims>
    struct collect_used_codims
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp)
      {
        node.collect_used_codims(codims);
      }

      collect_used_codims(Codims& codims_)
        : codims(codims_)
      {}

      Codims& codims;

    };


    struct collect_a_priori_fixed_size
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp)
      {
        node.update_a_priori_fixed_size();
        any = any || node._fixed_size;
        all = all && node._fixed_size;
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
        , all(true)
      {}

      bool any;
      bool all;

    };


    template<typename ES>
    struct update_fixed_size
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp) const
      {
        if (node._fixed_size)
          {
            typedef typename Node::Traits::SizeType size_type;
            const size_type dim = ES::dimension;
            node._codim_used.reset();
            node._gt_used.assign(GlobalGeometryTypeIndex::size(dim),false);
            node._gt_dof_offsets.assign(GlobalGeometryTypeIndex::size(dim),0);
            for (const auto& gt : es.indexSet().types())
              {
                size_type size = node.finiteElementMap().size(gt);
                node._gt_dof_offsets[GlobalGeometryTypeIndex::index(gt)] = size;
                node._gt_used[GlobalGeometryTypeIndex::index(gt)] = size > 0;
                node._codim_used[dim - gt.dim()] = node._codim_used[dim - gt.dim()] || (size > 0);
              }
            node._max_local_size = node.finiteElementMap().maxLocalSize();
          }
      }

      template<typename Node, typename TreePath>
      void pre(Node& node, TreePath tp) const
      {
        if (node._fixed_size)
          {
            typedef typename Node::Traits::SizeType size_type;
            const size_type dim = ES::dimension;
            node._codim_used.reset();
            node._gt_used.assign(Dune::GlobalGeometryTypeIndex::size(dim),false);
            node._gt_dof_offsets.assign(Dune::GlobalGeometryTypeIndex::size(dim) * TypeTree::degree(node),0);
            node._max_local_size = 0;
          }
      }

      template<typename Node, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(Node& node, const Child& child, TreePath tp, ChildIndex childIndex) const
      {
        if (node._fixed_size)
          {
            node._codim_used |= child._codim_used;

            std::transform(node._gt_used.begin(),
                           node._gt_used.end(),
                           child._gt_used.begin(),
                           node._gt_used.begin(),
                           std::logical_or<bool>());

            node._max_local_size += child._max_local_size;

            typedef typename Node::Traits::SizeType size_type;

            const size_type per_gt_size = child._child_count > 0 ? child._child_count : 1;
            const size_type size_offset = child._child_count > 0 ? child._child_count - 1 : 0;

            for (size_type gt = 0; gt < Dune::GlobalGeometryTypeIndex::size(ES::dimension); ++gt)
              node._gt_dof_offsets[gt * TypeTree::degree(node) + childIndex] = child._gt_dof_offsets[gt * per_gt_size + size_offset];
          }
      }

      template<typename Node, typename TreePath>
      void post(Node& node, TreePath tp) const
      {
        if (node._fixed_size)
          {
            typedef typename std::vector<typename Node::Traits::SizeType>::iterator iterator;

            iterator next_gt_it = node._gt_dof_offsets.begin() + TypeTree::degree(node);
            const iterator end_it = node._gt_dof_offsets.end();

            for (iterator it = node._gt_dof_offsets.begin();
                 it != end_it;
                 it += TypeTree::degree(node), next_gt_it += TypeTree::degree(node))
              std::partial_sum(it,next_gt_it,it);
          }
      }

      update_fixed_size(const ES es_)
        : es(es_)
      {}

      ES es;

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
            node._codim_used.reset();
            node._gt_used.assign(Dune::GlobalGeometryTypeIndex::size(dim),false);
            node._gt_dof_offsets.assign(Dune::GlobalGeometryTypeIndex::size(dim) * std::max(node._child_count,static_cast<std::size_t>(1)),0);
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
    struct collect_used_geometry_types_from_cell_visitor
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp) const
      {
        if (!node._fixed_size)
          node.collect_used_geometry_types_from_cell(cell);
      }

      collect_used_geometry_types_from_cell_visitor(const Cell& cell_)
        : cell(cell_)
        , ref_el(Dune::ReferenceElements<typename Cell::Geometry::ctype,Cell::dimension>::general(cell_.type()))
      {}

      const Cell& cell;
      Dune::ReferenceElement<typename Cell::Geometry> ref_el;

    };


    template<typename ES>
    struct post_collect_used_geometry_types
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp) const
      {
        if (!node._fixed_size)
          {
            typedef typename Node::Traits::SizeType size_type;

            for (const auto& gt : es.indexSet().types())
              {
                if (node._gt_used[Dune::GlobalGeometryTypeIndex::index(gt)])
                  node._gt_entity_offsets[Dune::GlobalGeometryTypeIndex::index(gt) + 1] = es.indexSet().size(gt);
              }

            std::partial_sum(node._gt_entity_offsets.begin(),node._gt_entity_offsets.end(),node._gt_entity_offsets.begin());
            node._entity_dof_offsets.assign(node._gt_entity_offsets.back() * std::max(node._child_count,static_cast<size_type>(1)),0);
            node.setup_fixed_size_possible();
          }
      }

      template<typename Node, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(Node& node, const Child& child, TreePath tp, ChildIndex childIndex) const
      {
        if (!node._fixed_size)
          {
            node._codim_used |= child._codim_used;

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

      post_collect_used_geometry_types(const ES& es_)
        : es(es_)
      {}

      ES es;

    };


    template<typename ES>
    struct extract_per_entity_sizes_from_cell_visitor
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      static const std::size_t dim = ES::dimension;
      typedef typename ES::template Codim<0>::Entity Cell;
      typedef std::size_t size_type;

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp)
      {
        if (!node._fixed_size)
          node.extract_per_entity_sizes_from_cell(*cell,gt_sizes);
      }

      extract_per_entity_sizes_from_cell_visitor(const ES& es_)
        : es(es_)
        , cell(nullptr)
        , ref_el()
        , gt_sizes(Dune::GlobalGeometryTypeIndex::size(dim),0)
      {}

      void set_cell(const Cell& cell_)
      {
        cell = &cell_;
        ref_el = referenceElement(cell_.geometry());
      }

      ES es;
      const Cell* cell;
      Dune::ReferenceElement<typename Cell::Geometry> ref_el;
      std::vector<size_type> gt_sizes;

    };


    template<typename ES>
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
          {
            node._fixed_size_possible = true;
            node._max_local_size = 0;
          }
      }


      template<typename Node, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(Node& node, const Child& child, TreePath tp, ChildIndex childIndex) const
      {
        if (!node._fixed_size)
          {
            node._fixed_size_possible = node._fixed_size_possible && child._fixed_size;
            node._max_local_size += child._max_local_size;
          }
      }


      template<typename Node, typename TreePath>
      void post(Node& node, TreePath tp) const
      {
        if (!node._fixed_size)
          {

            typedef typename Node::Traits::SizeType size_type;
            const size_type dim = ES::dimension;

            if (node._fixed_size_possible)
              {

                for (size_type gt = 0; gt < GlobalGeometryTypeIndex::size(ES::dimension); ++gt)
                  {
                    for (size_type child_index = 0; child_index < TypeTree::degree(node); ++child_index)
                      {
                        const size_type per_gt_size = node.childOrdering(child_index)._child_count > 0 ? node.childOrdering(child_index)._child_count : 1;
                        const size_type size_offset = node.childOrdering(child_index)._child_count > 0 ? node.childOrdering(child_index)._child_count - 1 : 0;

                        node._gt_dof_offsets[gt * TypeTree::degree(node) + child_index] = node.childOrdering(child_index)._gt_dof_offsets[gt * per_gt_size + size_offset];
                      }
                  }

                typedef typename std::vector<typename Node::Traits::SizeType>::iterator iterator;

                const iterator end_it = node._gt_dof_offsets.end();

                for (iterator it = node._gt_dof_offsets.begin();
                     it != end_it;
                     it += TypeTree::degree(node))
                  std::partial_sum(it,it + TypeTree::degree(node),it);

                node._fixed_size = true;
              }
            else
              {
                typedef typename Node::Traits::SizeType size_type;

                size_type index = 0;
                for (size_type geometry_type_index = 0; geometry_type_index < GlobalGeometryTypeIndex::size(dim); ++geometry_type_index)
                  {
                    if (!node._gt_used[geometry_type_index])
                      continue;
                    const size_type entity_count = node._gt_entity_offsets[geometry_type_index+1] - node._gt_entity_offsets[geometry_type_index];
                    for (size_type entity_index = 0; entity_index < entity_count; ++entity_index)
                      {
                        size_type carry = 0;
                        for (size_type child_index = 0; child_index < TypeTree::degree(node); ++child_index)
                          node._entity_dof_offsets[index++] = (carry += node.childOrdering(child_index).size(geometry_type_index,entity_index));
                      }
                  }

              }
          }
      }

      post_extract_per_entity_sizes(const ES& es_)
        : es(es_)
      {}

      ES es;

    };


    template<typename LocalOrdering>
    class GridViewOrdering
      : public TypeTree::CompositeNode<LocalOrdering>
      , public VirtualOrderingBase<typename LocalOrdering::Traits::DOFIndex,
                                   typename LocalOrdering::Traits::ContainerIndex>
      , public OrderingBase<typename LocalOrdering::Traits::DOFIndex,
                            typename LocalOrdering::Traits::ContainerIndex>
    {
    public:
      typedef typename LocalOrdering::Traits Traits;

      static const bool has_dynamic_ordering_children = false;

      static const bool consume_tree_index = false;

    private:

      typedef TypeTree::CompositeNode<LocalOrdering> NodeT;
      typedef OrderingBase<
        typename LocalOrdering::Traits::DOFIndex,
        typename LocalOrdering::Traits::ContainerIndex
        > BaseT;

      using EntitySet = typename Traits::EntitySet;

    public:
      //! Construct ordering object
      /**
       * In general, an ordering object is not properly setup after
       * construction.  This must be done by a separate call to update().
       * This particular ordering however can be used right away.
       */
      GridViewOrdering(const typename NodeT::NodeStorage& local_ordering, bool container_blocked, typename BaseT::GFSData* gfs_data)
        : NodeT(local_ordering)
        , BaseT(*this,container_blocked,gfs_data,this)
        , _es(localOrdering().entitySet())
      {
        // make sure to switch off container blocking handling in the local ordering,
        // we already handle it in the GridViewOrdering
        localOrdering().disable_container_blocking();
      }

#ifndef DOXYGEN

// we need to override the default copy / move ctor to fix the delegate pointer, but that is
// hardly interesting to our users...

      GridViewOrdering(const GridViewOrdering& r)
        : NodeT(r.nodeStorage())
        , BaseT(r)
        , _es(r._es)
        , _gt_dof_offsets(r._gt_dof_offsets)
        , _gt_entity_offsets(r._gt_entity_offsets)
        , _entity_dof_offsets(r._entity_dof_offsets)
      {
        this->setDelegate(this);
      }

      GridViewOrdering(GridViewOrdering&& r)
        : NodeT(r.nodeStorage())
        , BaseT(std::move(r))
        , _es(std::move(r._es))
        , _gt_dof_offsets(std::move(r._gt_dof_offsets))
        , _gt_entity_offsets(std::move(r._gt_entity_offsets))
        , _entity_dof_offsets(std::move(r._entity_dof_offsets))
      {
        this->setDelegate(this);
      }

#endif // DOXYGEN

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
        mapIndex(di,ci);
      }

      typename Traits::ContainerIndex mapIndex(const typename Traits::DOFIndex& di) const
      {
        typename Traits::ContainerIndex ci;
        mapIndex(di.view(),ci);
        return ci;
      }

      void mapIndex(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {
        typedef typename Traits::SizeType size_type;
        const size_type geometry_type_index = Traits::DOFIndexAccessor::geometryType(di);
        const size_type entity_index = Traits::DOFIndexAccessor::entityIndex(di);
        localOrdering().map_local_index(geometry_type_index,entity_index,di.treeIndex(),ci);
        if (_container_blocked)
          {
            if (_fixed_size)
              {
                ci.push_back(_gt_dof_offsets[geometry_type_index] + entity_index);
              }
            else
              {
                ci.push_back(_gt_entity_offsets[geometry_type_index] + entity_index);
              }
          }
        else
          {
            if (_fixed_size)
              {
                ci.back() += _gt_dof_offsets[geometry_type_index] + entity_index * localOrdering().size(geometry_type_index,entity_index);
              }
            else
              {
                ci.back() += _entity_dof_offsets[_gt_entity_offsets[geometry_type_index] + entity_index];
              }
          }
      }

      template<typename ItIn, typename ItOut>
      void map_lfs_indices(const ItIn begin, const ItIn end, ItOut out) const
      {
        typedef typename Traits::SizeType size_type;
        if (_container_blocked)
          {
            if (_fixed_size)
              for (ItIn in = begin; in != end; ++in, ++out)
                {
                  const size_type geometry_type_index = Traits::DOFIndexAccessor::geometryType(*in);
                  const size_type entity_index = Traits::DOFIndexAccessor::entityIndex(*in);
                  out->push_back(_gt_dof_offsets[geometry_type_index] + entity_index);
                }
            else
              for (ItIn in = begin; in != end; ++in, ++out)
                {
                  const size_type geometry_type_index = Traits::DOFIndexAccessor::geometryType(*in);
                  const size_type entity_index = Traits::DOFIndexAccessor::entityIndex(*in);
                  out->push_back(_gt_entity_offsets[geometry_type_index] + entity_index);
                }
          }
        else if (_fixed_size)
          {
            for (ItIn in = begin; in != end; ++in, ++out)
              {
                const size_type geometry_type_index = Traits::DOFIndexAccessor::geometryType(*in);
                const size_type entity_index = Traits::DOFIndexAccessor::entityIndex(*in);
                out->back() += _gt_dof_offsets[geometry_type_index] + entity_index * localOrdering().size(geometry_type_index,entity_index);
              }
          }
        else
          {
            for (ItIn in = begin; in != end; ++in, ++out)
              {
                const size_type geometry_type_index = Traits::DOFIndexAccessor::geometryType(*in);
                const size_type entity_index = Traits::DOFIndexAccessor::entityIndex(*in);
                out->back() += _entity_dof_offsets[_gt_entity_offsets[geometry_type_index] + entity_index];
              }
          }
      }

      template<typename CIOutIterator>
      typename Traits::SizeType
      extract_entity_indices(const typename Traits::DOFIndex::EntityIndex& ei,
                             typename Traits::SizeType child_index,
                             CIOutIterator ci_out, const CIOutIterator ci_end) const
      {
        typedef typename Traits::SizeType size_type;

        const size_type geometry_type_index = Traits::DOFIndexAccessor::GeometryIndex::geometryType(ei);
        const size_type entity_index = Traits::DOFIndexAccessor::GeometryIndex::entityIndex(ei);

        if (_container_blocked)
          {
            if (_fixed_size)
              for (; ci_out != ci_end; ++ci_out)
                {
                  ci_out->push_back(_gt_dof_offsets[geometry_type_index] + entity_index);
                }
            else
              for (; ci_out != ci_end; ++ci_out)
                {
                  ci_out->push_back(_gt_entity_offsets[geometry_type_index] + entity_index);
                }
          }
        else if (_fixed_size)
          {
            for (; ci_out != ci_end; ++ci_out)
              {
                ci_out->back() += _gt_dof_offsets[geometry_type_index] + entity_index * localOrdering().size(geometry_type_index,entity_index);
              }
          }
        else
          {
            for (; ci_out != ci_end; ++ci_out)
              {
                ci_out->back() += _entity_dof_offsets[_gt_entity_offsets[geometry_type_index] + entity_index];
              }
          }

        // The return value is not used for non-leaf orderings.
        return 0;
      }

      void update()
      {

        typedef typename Traits::SizeType size_type;
        using ES = typename Traits::EntitySet;
        const size_type dim = ES::dimension;

        typename ES::CodimMask codims;
        codims.set(0); // we always need cells

        TypeTree::applyToTree(localOrdering(),collect_used_codims<typename ES::CodimMask>(codims));

        for (typename ES::dim_type codim = 0; codim <= ES::dimension; ++codim)
          if (codims.test(codim))
            _es.addCodim(codim);

        _es.update();

        // Do we already know that we have fixed per-GeometryType sizes?
        collect_a_priori_fixed_size fixed_size_collector;
        TypeTree::applyToTree(localOrdering(),fixed_size_collector);
        _fixed_size = localOrdering().fixedSize();

        const size_type gt_index_count = GlobalGeometryTypeIndex::size(ES::dimension);

        if (fixed_size_collector.any)
          {
            // collect used GeometryTypes
            TypeTree::applyToTree(localOrdering(),update_fixed_size<ES>(_es));
          }

        if (!fixed_size_collector.all)
          {
            TypeTree::applyToTree(localOrdering(),pre_collect_used_geometry_types(ES::dimension));

            using Element = typename ES::template Codim<0>::Entity;

            for (const auto& element : elements(_es))
              {
                TypeTree::applyToTree(localOrdering(),collect_used_geometry_types_from_cell_visitor<Element>(element));
              }
            TypeTree::applyToTree(localOrdering(),post_collect_used_geometry_types<ES>(_es));
            // allocate

            //TypeTree::applyToTree(localOrdering(),pre_extract_per_entity_sizes<GV>(_gv));
            extract_per_entity_sizes_from_cell_visitor<ES> visitor(_es);
            for (const auto& element : elements(_es))
              {
                visitor.set_cell(element);
                TypeTree::applyToTree(localOrdering(),visitor);
              }
            TypeTree::applyToTree(localOrdering(),post_extract_per_entity_sizes<ES>(_es));
          }

        _codim_used = localOrdering()._codim_used;

        if (localOrdering().fixedSize())
          {
            _fixed_size = true;
            _gt_dof_offsets.assign(gt_index_count + 1,0);

            _block_count = 0;

            _size = 0;

            for (const auto& gt : _es.indexSet().types())
              {
                const size_type gt_index = GlobalGeometryTypeIndex::index(gt);
                size_type gt_size = localOrdering().size(gt_index,0);
                const size_type gt_entity_count = _es.indexSet().size(gt);
                _size += gt_size * gt_entity_count;
                if (_container_blocked)
                  gt_size = gt_size > 0;
                _gt_dof_offsets[gt_index + 1] = gt_size * gt_entity_count;
              }

            std::partial_sum(_gt_dof_offsets.begin(),_gt_dof_offsets.end(),_gt_dof_offsets.begin());
            _block_count = _gt_dof_offsets.back();

            _codim_fixed_size.set();

          }
        else
          {
            _gt_entity_offsets.assign(gt_index_count + 1,0);

            for (const auto& gt : _es.indexSet().types())
              {
                if (!localOrdering().contains(gt))
                  continue;
                const size_type gt_index = GlobalGeometryTypeIndex::index(gt);
                _gt_entity_offsets[gt_index + 1] = _es.indexSet().size(gt);
              }

            std::partial_sum(_gt_entity_offsets.begin(),_gt_entity_offsets.end(),_gt_entity_offsets.begin());
            _entity_dof_offsets.assign(_gt_entity_offsets.back()+1,0);
            _block_count = 0;

            size_type carry = 0;
            size_type index = 0;
            for (size_type gt_index = 0; gt_index < GlobalGeometryTypeIndex::size(dim); ++gt_index)
              {
                if (!localOrdering().contains_geometry_type(gt_index))
                  continue;
                const size_type entity_count = _gt_entity_offsets[gt_index + 1] - _gt_entity_offsets[gt_index];
                for (size_type entity_index = 0; entity_index < entity_count; ++entity_index)
                  {
                    const size_type size = localOrdering().size(gt_index,entity_index);
                    _entity_dof_offsets[++index] = (carry += size);
                    _block_count += (size > 0);
                  }
              }
            _size = _entity_dof_offsets.back();

            if (!_container_blocked)
              _block_count = _size;

            _codim_fixed_size.reset();
          }

        _max_local_size = localOrdering().maxLocalSize();
      }

      using BaseT::fixedSize;

    private:

      using BaseT::_container_blocked;
      using BaseT::_fixed_size;
      using BaseT::_max_local_size;
      using BaseT::_size;
      using BaseT::_block_count;
      using BaseT::_codim_used;
      using BaseT::_codim_fixed_size;

      typename Traits::EntitySet _es;
      std::vector<typename Traits::SizeType> _gt_dof_offsets;
      std::vector<typename Traits::SizeType> _gt_entity_offsets;
      std::vector<typename Traits::SizeType> _entity_dof_offsets;

    };


   //! \} group Ordering
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_GRIDVIEWORDERING_HH
