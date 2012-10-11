// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMMON_DOFINDEX_HH
#define DUNE_PDELAB_COMMON_DOFINDEX_HH

#include <dune/pdelab/common/multiindex.hh>

namespace Dune {

  namespace PDELab {

    template<typename T, std::size_t tree_n, std::size_t entity_n = 1>
    class DOFIndex
    {

    public:

      //! The maximum possible depth of the MultiIndex.
      static const std::size_t max_depth = tree_n;
      static const std::size_t entity_capacity = entity_n;

      typedef array<T,entity_capacity> EntityIndex;
      typedef MultiIndex<T,max_depth> TreeIndex;

      typedef typename TreeIndex::size_type size_type;
      typedef T value_type;

      class View
      {

        friend class DOFIndex;

      public:

        static const std::size_t max_depth = tree_n;
        static const std::size_t entity_capacity = entity_n;

        typedef const array<T,entity_n>& EntityIndex;
        typedef typename MultiIndex<T,tree_n>::View TreeIndex;

        const EntityIndex& entityIndex() const
        {
          return _entity_index_view;
        }

        const TreeIndex& treeIndex() const
        {
          return _tree_index_view;
        }

        View back_popped() const
        {
          return View(_entity_index_view,_tree_index_view.back_popped());
        }

        friend std::ostream& operator<< (std::ostream& s, const View& di)
        {
          s << "(";

          for (typename std::remove_reference<EntityIndex>::type::const_iterator it = di._entity_index_view.begin(); it != di._entity_index_view.end(); ++it)
            s << std::setw(4) << *it;

          s << " | "
            << di._tree_index_view
            << ")";
          return s;
        }

        std::size_t size() const
        {
          return _tree_index_view.size();
        };

      private:

        explicit View(const DOFIndex& dof_index)
          : _entity_index_view(dof_index._entity_index)
          , _tree_index_view(dof_index._tree_index.view())
        {}

        View(const DOFIndex& dof_index, std::size_t size)
          : _entity_index_view(dof_index._entity_index)
          , _tree_index_view(dof_index._tree_index.view(size))
        {}

        View(const EntityIndex& entity_index, const TreeIndex& tree_index)
          : _entity_index_view(entity_index)
          , _tree_index_view(tree_index)
        {}

        EntityIndex _entity_index_view;
        TreeIndex _tree_index_view;

      };

      //! Default constructor.
      DOFIndex()
      {}

      View view() const
      {
        return View(*this);
      }

      View view(std::size_t size) const
      {
        return View(*this,size);
      }

      void clear()
      {
        std::fill(_entity_index.begin(),_entity_index.end(),0);
        _tree_index.clear();
      }

      //! Returns the index of the grid entity associated with the DOF.
      EntityIndex& entityIndex()
      {
        return _entity_index;
      }

      const EntityIndex& entityIndex() const
      {
        return _entity_index;
      }

      TreeIndex& treeIndex()
      {
        return _tree_index;
      }

      const TreeIndex& treeIndex() const
      {
        return _tree_index;
      }

      //! Writes a pretty representation of the MultiIndex to the given std::ostream.
      friend std::ostream& operator<< (std::ostream& s, const DOFIndex& di)
      {
        s << "(";

        for (typename EntityIndex::const_iterator it = di._entity_index.begin(); it != di._entity_index.end(); ++it)
          s << std::setw(4) << *it;

        s << " | "
          << di._tree_index
          << ")";
        return s;
      }

      //! Tests whether two MultiIndices are equal.
      /**
       * \note Only MultiIndices of identical max_depth are comparable
       */
      bool operator== (const DOFIndex& r) const
      {
        return
          std::equal(_entity_index.begin(),_entity_index.end(),r._entity_index.begin()) &&
          _tree_index == r._tree_index;
      }

      //! Tests whether two MultiIndices are not equal.
      bool operator!= (const DOFIndex& r) const
      {
        return !(*this == r);
      }

#if 0
      bool operator< (const DOFIndex& r) const
      {
        // FIXME: think about natural ordering
        return _c.size() < _r.size();
        return std::lexicographical_compare(_c.begin(),_c.end(),r._c.begin(),r._c.end());
      }
#endif

      std::size_t size() const
      {
        return _tree_index.size();
      }

    private:

      EntityIndex _entity_index;
      TreeIndex _tree_index;

    };

    template<typename T, std::size_t n1, std::size_t n2>
    inline std::size_t hash_value(const DOFIndex<T,n1,n2>& di)
    {
      std::size_t seed = 0;
      hash_range(seed,di.entityIndex().begin(),di.entityIndex().end());
      hash_range(seed,di.treeIndex().begin(),di.treeIndex().end());
      return seed;
    }



  } // namespace PDELab
} // namespace Dune

DUNE_DEFINE_HASH(DUNE_HASH_TEMPLATE_ARGS(typename T, std::size_t n1, std::size_t n2),DUNE_HASH_TYPE(Dune::PDELab::DOFIndex<T,n1,n2>))

#endif // DUNE_PDELAB_COMMON_DOFINDEX_HH
