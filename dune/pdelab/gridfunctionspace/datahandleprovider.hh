// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_DATAHANDLEPROVIDER_HH
#define DUNE_PDELAB_DATAHANDLEPROVIDER_HH

#include <vector>
#include <stack>

#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/reservedvector.hh>


namespace Dune {
  namespace PDELab {

    namespace {

      template<typename EntityIndex>
      struct get_size_for_entity
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename Ordering, typename TreePath>
        void leaf(const Ordering& ordering, TreePath tp)
        {
          _size += ordering.size(_entity_index);
        }

        get_size_for_entity(const EntityIndex& entity_index)
          : _size(0)
          , _entity_index(entity_index)
        {}

        std::size_t size() const
        {
          return _size;
        }

      private:

        std::size_t _size;
        const EntityIndex& _entity_index;

      };


      template<typename EntityIndex, typename ContainerIndex, std::size_t tree_depth>
      struct indices_for_entity
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        typedef std::size_t size_type;

        template<typename Ordering, typename Child, typename TreePath, typename ChildIndex>
        void beforeChild(const Ordering& ordering, const Child& child, TreePath tp, ChildIndex childIndex)
        {
          _stack.push(_it);
        }

        template<typename Ordering, typename TreePath>
        void leaf(const Ordering& ordering, TreePath tp)
        {
          size_type size = ordering.containerIndices(_entity_index,
                                                     tp.back(),
                                                     _it,
                                                     _end);

          _end += size;
          _it = _end;
        }

        template<typename Ordering, typename Child, typename TreePath, typename ChildIndex>
        void afterChild(const Ordering& ordering, const Child& child, TreePath tp, ChildIndex childIndex)
        {
          // pop
          ordering.containerIndices(_entity_index,
                                    childIndex,
                                    _stack.top(),
                                    _end);
          _stack.pop();
        }


        indices_for_entity(const EntityIndex& entity_index,
                           std::vector<ContainerIndex>& indices)
          : _entity_index(entity_index)
          , _indices(indices)
          , _it(indices.begin())
          , _end(indices.begin())
        {}

      private:

        typedef typename std::vector<ContainerIndex>::iterator Iterator;

        const EntityIndex& _entity_index;
        std::vector<ContainerIndex>& _indices;
        Iterator _it;
        Iterator _end;

        std::stack<Iterator,ReservedVector<Iterator,tree_depth> > _stack;
      };

    } // anonymous namespace


    template<typename GFS>
    class DataHandleProvider
    {

    public:

      typedef std::size_t size_type;

      //------------------------------
      // generic data handle interface
      //------------------------------

      //! returns true if data for this codim should be communicated
      bool dataHandleContains (int dim, int codim) const
      {
        return gfs().ordering().contains(codim);
      }

      //! returns true if size per entity of given dim and codim is a constant
      bool dataHandleFixedSize (int dim, int codim) const
      {
        return gfs().ordering().fixedSize(codim);
      }

      /*! how many objects of type DataType have to be sent for a given entity

        Note: Only the sender side needs to know this size.
      */
      template<class EntityType>
      size_type dataHandleSize (const EntityType& e) const
      {
        typedef typename GFS::Ordering Ordering;

        typedef typename Ordering::Traits::DOFIndex::EntityIndex EntityIndex;
        EntityIndex ei;

        Ordering::Traits::DOFIndexAccessor::GeometryIndex::store(
          ei,
          e.type(),
          gfs().gridView().indexSet().index(e)
        );

        get_size_for_entity<EntityIndex> get_size(ei);
        TypeTree::applyToTree(gfs().ordering(),get_size);

        return get_size.size();
      }

      //! return vector of global indices associated with the given entity
      template<class EntityType, typename ContainerIndex>
      size_type dataHandleContainerIndices (const EntityType& e,
                                            std::vector<ContainerIndex>& indices) const
      {
        typedef typename GFS::Ordering Ordering;

        dune_static_assert((is_same<ContainerIndex,typename Ordering::Traits::ContainerIndex>::value),
                           "dataHandleContainerIndices() called with invalid ContainerIndex type.");

        // Clear index state
        for (typename std::vector<ContainerIndex>::iterator it = indices.begin(),
               endit = indices.end();
             it != endit;
             ++it)
          it->clear();

        typedef typename Ordering::Traits::DOFIndex::EntityIndex EntityIndex;
        EntityIndex ei;

        Ordering::Traits::DOFIndexAccessor::GeometryIndex::store(
          ei,
          e.type(),
          gfs().gridView().indexSet().index(e)
        );

        get_size_for_entity<EntityIndex> get_size(ei);
        TypeTree::applyToTree(gfs().ordering(),get_size);

        indices.resize(get_size.size());

        indices_for_entity<
          EntityIndex,
          ContainerIndex,
          TypeTree::TreeInfo<Ordering>::depth
          > extract_indices(ei,indices);
        TypeTree::applyToTree(gfs().ordering(),extract_indices);

        return get_size.size();
      }

    private:

      const GFS& gfs() const
      {
        return static_cast<const GFS&>(*this);
      }

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_DATAHANDLEPROVIDER
