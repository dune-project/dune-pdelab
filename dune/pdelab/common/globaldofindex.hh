// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMMON_GLOBALDOFINDEX_HH
#define DUNE_PDELAB_COMMON_GLOBALDOFINDEX_HH

#include <dune/pdelab/common/multiindex.hh>

namespace Dune {

  namespace PDELab {


    template<typename T, std::size_t tree_n, typename ID>
    class GlobalDOFIndex
    {

    public:

      //! The maximum possible depth of the MultiIndex.
      static const std::size_t max_depth = tree_n;

      typedef ID EntityID;
      typedef MultiIndex<T,max_depth> TreeIndex;

      typedef typename TreeIndex::size_type size_type;
      typedef T value_type;

      class View
      {

        friend class GlobalDOFIndex;

      public:

        static const std::size_t max_depth = tree_n;

        typedef ID EntityID;
        typedef typename MultiIndex<T,tree_n>::View TreeIndex;

        const EntityID& entityID() const
        {
          return _entity_id;
        }

        const TreeIndex& treeIndex() const
        {
          return _tree_index_view;
        }

        View back_popped() const
        {
          return View(_entity_id,_tree_index_view.back_popped());
        }

        friend std::ostream& operator<< (std::ostream& s, const View& di)
        {
          s << "(";

          s << di._entity_id;

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

        explicit View(const GlobalDOFIndex& dof_index)
          : _entity_id(dof_index._entity_id)
          , _tree_index_view(dof_index._tree_index.view())
        {}

        View(const GlobalDOFIndex& dof_index, std::size_t size)
          : _entity_id(dof_index._entity_id)
          , _tree_index_view(dof_index._tree_index.view(size))
        {}

        View(const EntityID& entity_id, const TreeIndex& tree_index)
          : _entity_id(entity_id)
          , _tree_index_view(tree_index)
        {}

        const EntityID& _entity_id;
        TreeIndex _tree_index_view;

      };

      //! Default constructor.
      GlobalDOFIndex()
      {}

      GlobalDOFIndex(const EntityID& entity_id, const TreeIndex& tree_index)
        : _entity_id(entity_id)
        , _tree_index(tree_index)
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
        _entity_id = EntityID();
        _tree_index.clear();
      }

      //! Returns the index of the grid entity associated with the DOF.
      EntityID& entityID()
      {
        return _entity_id;
      }

      const EntityID& entityID() const
      {
        return _entity_id;
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
      friend std::ostream& operator<< (std::ostream& s, const GlobalDOFIndex& di)
      {
        s << "(";

        s << di._entity_id;

        s << " | "
          << di._tree_index
          << ")";
        return s;
      }

      //! Tests whether two MultiIndices are equal.
      /**
       * \note Only MultiIndices of identical max_depth are comparable
       */
      bool operator== (const GlobalDOFIndex& r) const
      {
        return
          _entity_id == r._entity_id && _tree_index == r._tree_index;
      }

      //! Tests whether two MultiIndices are not equal.
      bool operator!= (const GlobalDOFIndex& r) const
      {
        return !(*this == r);
      }


      std::size_t size() const
      {
        return _tree_index.size();
      }

    private:

      EntityID _entity_id;
      TreeIndex _tree_index;

    };

    template<typename T, std::size_t n, typename ID>
    inline std::size_t hash_value(const GlobalDOFIndex<T,n,ID>& di)
    {
      std::size_t seed = 0;
      std::hash<ID> id_hasher;
      hash_combine(seed,id_hasher(di.entityID()));
      hash_range(seed,di.treeIndex().begin(),di.treeIndex().end());
      return seed;
    }



  } // namespace PDELab
} // namespace Dune

DUNE_DEFINE_HASH(DUNE_HASH_TEMPLATE_ARGS(typename T, std::size_t n, typename ID),DUNE_HASH_TYPE(Dune::PDELab::GlobalDOFIndex<T,n,ID>))

#endif // DUNE_PDELAB_COMMON_GLOBALDOFINDEX_HH
