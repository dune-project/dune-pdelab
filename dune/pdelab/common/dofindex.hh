// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMMON_DOFINDEX_HH
#define DUNE_PDELAB_COMMON_DOFINDEX_HH

#include <dune/pdelab/common/multiindex.hh>
#include <dune/common/array.hh>

namespace Dune {

  namespace PDELab {

    //! A multi-index representing a degree of freedom in a GridFunctionSpace.
    /**
     * A DOFIndex provides a way for identifying degrees of freedom in a (possibly
     * nested) GridFunctionSpace without imposing any kind of ordering. For that
     * purpose, a DOFIndex identifies a degree of freedom by recording
     *
     * - the geometry type of the grid entity associated with the DOF,
     * - an index value uniquely identifying the grid entity among all grid entities
     *   of its geometry type (usually just the index value of some Grid IndexSet),
     * - a tuple of entity-local indices.
     *
     * The length of this index tuple is limited by the template parameter \p tree_n, which will
     * usually be equal to the maximum depth of the current GridFunctionSpace tree. Moreover,
     * there will never be two identical index tuples associated with the same grid entity.
     *
     * The index tuple is oriented from left to right when traversing up the tree (i.e.
     * towards the root node) and from right to left when drilling down from the root node
     * towards a leaf. For non-leaf nodes, the associated index entry identifies the child
     * GridFunctionSpace that the degree of freedom is associated with, while for leafs, it
     * provides a way to provide multiple degrees of freedom for a single grid entity (usually,
     * the index value for a leaf space will correspond to the LocalKey::index() value from
     * the finite element).
     *
     * Note that in general, the length of the index tuple will not be the same for all degrees
     * of freedom in a GridFunctionSpace. Consider the following example of a Taylor-Hood element
     * (P<sub>2</sub>-P<sub>1</sub>):
     * \dot
     * graph taylor_hood {
     * node [shape=record, style=rounded, fontname=Helvetica, fontsize=8, height=0.2, width=0.4];
     * TH [ label="Taylor-Hood"];
     * TH -- V;
     * V [ label="Velocity"];
     * TH -- P;
     * P [ label="Pressure"];
     * V -- Vx;
     * Vx [ label="x Velocity" ];
     * V -- Vy;
     * Vy [ label="y Velocity" ];
     * }
     * \enddot
     * In this case, degrees of freedom for the velocity components will have an index tuple of length
     * 3, while those related to pressure will only have an index tuple of length 2. For the Taylor-Hood
     * space given above, the DOFIndex associated with a triangle with vertex and edge indices in the
     * range {0,1,2} are
     *
     * <table>
     *  <tr>
     *   <th>GeometryType</th>
     *   <th>grid entity index</th>
     *   <th>entity-local index tuple</th>
     *  </tr>
     *  <tr>
     *   <td>Point</td>
     *   <td>0</td>
     *   <td align="right">0, 0, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Point</td>
     *   <td>1</td>
     *   <td align="right">0, 0, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Point</td>
     *   <td>2</td>
     *   <td align="right">0, 0, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Line</td>
     *   <td>0</td>
     *   <td align="right">0, 0, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Line</td>
     *   <td>1</td>
     *   <td align="right">0, 0, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Line</td>
     *   <td>2</td>
     *   <td align="right">0, 0, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Point</td>
     *   <td>0</td>
     *   <td align="right">0, 1, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Point</td>
     *   <td>1</td>
     *   <td align="right">0, 1, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Point</td>
     *   <td>2</td>
     *   <td align="right">0, 1, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Line</td>
     *   <td>0</td>
     *   <td align="right">0, 1, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Line</td>
     *   <td>1</td>
     *   <td align="right">0, 1, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Line</td>
     *   <td>2</td>
     *   <td align="right">0, 1, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Point</td>
     *   <td>0</td>
     *   <td align="right">0, 1</td>
     *  </tr>
     *  <tr>
     *   <td>Point</td>
     *   <td>1</td>
     *   <td align="right">0, 1</td>
     *  </tr>
     *  <tr>
     *   <td>Point</td>
     *   <td>2</td>
     *   <td align="right">0, 1</td>
     *  </tr>
     * </table>
     *
     * \tparam T         the type of the index entries.
     * \tparam tree_n    the maximum length of the tuple of entity-local indices.
     * \tparam entity_n  the maximum length of a tuple which represents the entity index.
     */
    template<typename T, std::size_t tree_n, std::size_t entity_n = 1>
    class DOFIndex
    {

    public:

      //! The maximum possible length of a tuple which represents the entity index.
      static const std::size_t entity_capacity = entity_n;

      //! The maximum possible length of the tuple of entity-local indices.
      static const std::size_t max_depth = tree_n;

      typedef std::array<T,entity_capacity> EntityIndex;
      typedef MultiIndex<T,max_depth> TreeIndex;

      typedef typename TreeIndex::size_type size_type;
      typedef T value_type;

      class View
      {

        friend class DOFIndex;

      public:

        static const std::size_t max_depth = tree_n;
        static const std::size_t entity_capacity = entity_n;

        typedef const std::array<T,entity_n>& EntityIndex;
        typedef typename MultiIndex<T,tree_n>::View TreeIndex;

        EntityIndex& entityIndex() const
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

        View(EntityIndex& entity_index, const TreeIndex& tree_index)
          : _entity_index_view(entity_index)
          , _tree_index_view(tree_index)
        {}

        EntityIndex _entity_index_view;
        TreeIndex _tree_index_view;

      };

      //! Default constructor.
      DOFIndex()
      {}

      explicit DOFIndex(const View& view)
        : _entity_index(view._entity_index_view)
        , _tree_index(view._tree_index_view)
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

      //! Returns the index of the grid entity associated with the DOF.
      const EntityIndex& entityIndex() const
      {
        return _entity_index;
      }

      //! Returns the tuple of entity-local indices associated with the DOF.
      TreeIndex& treeIndex()
      {
        return _tree_index;
      }

      //! Returns the tuple of entity-local indices associated with the DOF.
      const TreeIndex& treeIndex() const
      {
        return _tree_index;
      }

      //! Writes a pretty representation of the DOFIndex to the given std::ostream.
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

      //! Tests whether two DOFIndices are equal.
      /**
       * \note Only DOFIndices of identical max_depth are comparable.
       */
      bool operator== (const DOFIndex& r) const
      {
        return
          std::equal(_entity_index.begin(),_entity_index.end(),r._entity_index.begin()) &&
          _tree_index == r._tree_index;
      }

      //! Tests whether two DOFIndices are not equal.
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
