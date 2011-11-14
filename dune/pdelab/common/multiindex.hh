// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMMON_MULTIINDEX_HH
#define DUNE_PDELAB_COMMON_MULTIINDEX_HH

#include <dune/common/reservedvector.hh>

#include <dune/geometry/type.hh>

#include <algorithm>
#include <iomanip>

namespace Dune {

  namespace PDELab {

    //! A multi-index representing a degree of freedom in a GridFunctionSpace.
    /**
     * A MultiIndex provides a way for identifying degrees of freedom in a (possibly
     * nested) GridFunctionSpace without imposing any kind of ordering. For that
     * purpose, a MultiIndex identifies a degree of freedom by recording
     *
     * - the geometry type of the grid entity associated with the DOF,
     * - an index value uniquely identifying the grid entity among all grid entities
     *   of its geometry type (usually just the index value of some Grid IndexSet),
     * - a tuple of entity-local indices.
     *
     * The length of this index tuple is limited by the template parameter \ref n, which will
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
     * of freedom in a GridFunctionSpace. Consider the following example of a Taylor-Hood element:
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
     * space given above, the multiindices associated to a triangle with vertex and edge indices in the
     * range {0,1,2} are
     *
     * <table>
     *  <tr>
     *   <th>GeometryType</th>
     *   <th>entity index</th>
     *   <th>index tuple</th>
     *  </tr>
     *  <tr>
     *   <td>Point</td>
     *   <td>0</td>
     *   <td>0, 0, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Point</td>
     *   <td>1</td>
     *   <td>0, 0, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Point</td>
     *   <td>2</td>
     *   <td>0, 0, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Line</td>
     *   <td>0</td>
     *   <td>0, 0, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Line</td>
     *   <td>1</td>
     *   <td>0, 0, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Line</td>
     *   <td>2</td>
     *   <td>0, 0, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Point</td>
     *   <td>0</td>
     *   <td>0, 1, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Point</td>
     *   <td>1</td>
     *   <td>0, 1, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Point</td>
     *   <td>2</td>
     *   <td>0, 1, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Line</td>
     *   <td>0</td>
     *   <td>0, 1, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Line</td>
     *   <td>1</td>
     *   <td>0, 1, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Line</td>
     *   <td>2</td>
     *   <td>0, 1, 0</td>
     *  </tr>
     *  <tr>
     *   <td>Point</td>
     *   <td>0</td>
     *   <td>0, 1</td>
     *  </tr>
     *  <tr>
     *   <td>Point</td>
     *   <td>1</td>
     *   <td>0, 1</td>
     *  </tr>
     *  <tr>
     *   <td>Point</td>
     *   <td>2</td>
     *   <td>0, 1</td>
     *  </tr>
     * </table>
     *
     * \tparam T  the type of the index entries.
     * \tparam n  the maximum number of indices in the MultiIndex.
     */
    template<typename T, std::size_t n>
    class MultiIndex
    {

      typedef ReservedVector<T,n> container_type;

    public:

      //! The maximum possible depth of the MultiIndex.
      static const std::size_t max_depth = n;

      //! Unsigned integral type.
      typedef typename container_type::size_type size_type;

      //! The value type of the individual index entries.
      typedef typename container_type::value_type value_type;

      //! A reference to an index entry.
      typedef typename container_type::reference reference;

      //! A const reference to an index entry.
      typedef typename container_type::const_reference const_reference;

      //! Iterator over index entries.
      typedef typename container_type::iterator iterator;

      //! Const iterator over index entries.
      typedef typename container_type::const_iterator const_iterator;

      //! Default constructor.
      MultiIndex()
      {}

      //! Constructor to be called from leaf GridFunctionSpaces.
      /**
       * \param gt            the GeometryType of the grid entity associated with the DOF.
       * \param entity_index  the index of the grid entity associated with the DOF.
       * \param index         the first index of the DOF within the MultiIndex.
       */
      MultiIndex(GeometryType gt, value_type entity_index, value_type index)
        : _gt(gt)
        , _entity_index(entity_index)
      {
        push_back(index);
      }

      //! Sets the MultiIndex to a new DOF.
      /**
       * This should only be called from a leaf GridFunctionSpace as it will delete all exisiting
       * index entries.
       *
       * \param gt            the GeometryType of the grid entity associated with the DOF.
       * \param entity_index  the index of the grid entity associated with the DOF.
       * \param index         the first index of the DOF within the MultiIndex.
       */
      void set(GeometryType gt, value_type entity_index, value_type index)
      {
        _gt = gt;
        _entity_index = entity_index;
        _c.clear();
        push_back(index);
      }

      //! Returns the GeometryType of the grid entity associated with the DOF.
      Dune::GeometryType geometryType() const
      {
        return _gt;
      }

      //! Returns the index of the grid entity associated with the DOF.
      value_type entityIndex() const
      {
        return _entity_index;
      }

      //! Appends a new index entry to the MultiIndex.
      void push_back(value_type i)
      {
        _c.push_back(i);
      }

      //! Writes a pretty representation of the MultiIndex to the given std::ostream.
      friend std::ostream& operator<< (std::ostream& s, const MultiIndex& mi)
      {
        s << "("
          // entity information
          << std::setw(12) << mi._gt << " " << std::setw(4) << mi._entity_index
          << " |";
        // index tuple
        for (const_iterator it = mi._c.begin(); it != mi._c.end(); ++it)
          s << std::setw(3) << *it;
        // fill up to maximum depth for consistent formatting
        for (std::size_t i = mi._c.size(); i < mi.max_depth; ++i)
          s << "  -";
        s << ")";
        return s;
      }

      //! Tests whether two MultiIndices are equal.
      /**
       * \note Only MultiIndices of identical max_depth are comparable
       */
      bool operator== (const MultiIndex& r) const
      {
        return
          _gt == r._gt &&
          _entity_index = r._entity_index &&
          std::equal(_c.begin(),_c.end(),r._c.begin());
      }

      //! Tests whether two MultiIndices are not equal.
      bool operator!= (const MultiIndex& r) const
      {
        return !(*this == r);
      }

#if 0
      bool operator< (const MultiIndex& r) const
      {
        // FIXME: think about natural ordering
        return _c.size() < _r.size();
        return std::lexicographical_compare(_c.begin(),_c.end(),r._c.begin(),r._c.end());
      }
#endif

      //! The size (number of index entries) of the MultiIndex.
      size_type size() const
      {
        return _c.size();
      }

      //! Returns a const_iterator pointing to the beginning of the MultiIndex.
      const_iterator begin() const
      {
        return _c.begin();
      }

      //! Returns a const_iterator pointing to the end of the MultiIndex.
      const_iterator end() const
      {
        return _c.end();
      }

      //! Returns an iterator pointing to the beginning of the MultiIndex.
      iterator begin()
      {
        return _c.begin();
      }

      //! Returns an iterator pointing to the end of the MultiIndex.
      iterator end()
      {
        return _c.end();
      }

      //! Returns a reference ot the i-th index entry.
      reference operator[](size_type i)
      {
        return _c[i];
      }

      //! Returns a const reference ot the i-th index entry.
      const_reference operator[](size_type i) const
      {
        return _c[i];
      }

    private:

      GeometryType _gt;
      value_type _entity_index;
      container_type _c;

    };

  } // namespace PDELab
} // namespace Dune


#endif // DUNE_PDELAB_COMMON_MULTIINDEX_HH
