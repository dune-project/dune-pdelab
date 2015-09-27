#ifndef DUNE_PDELAB_COMMON_ELEMENTMAPPER_HH
#define DUNE_PDELAB_COMMON_ELEMENTMAPPER_HH

#include <vector>
#include <algorithm>
#include <numeric>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>
#include <dune/grid/common/capabilities.hh>

namespace Dune {

  namespace PDELab {

#ifndef DOXYGEN

    // implementation for mixed grids
    template<typename GV, bool has_single_cell_type>
    class ElementMapperBase
    {

    protected:

      typedef typename GV::template Codim<0>::Entity Element;
      typedef std::size_t size_type;

    private:

      static const size_type dim = GV::dimension;
      typedef typename GV::IndexSet IndexSet;

    protected:

      void update()
      {
        // clear old values
        std::fill(_gt_offsets.begin(),_gt_offsets.end(),0);

        // extract per-GeometryType sizes in codim 0
        for (auto gt : _index_set.types(0))
          {
            _gt_offsets[LocalGeometryTypeIndex::index(gt) + 1] = _index_set.size(gt);
          }

        // convert to offsets
        std::partial_sum(_gt_offsets.begin(),_gt_offsets.end(),_gt_offsets.begin());
      }

      size_type map(const Element& e) const
      {
        return _gt_offsets[LocalGeometryTypeIndex::index(e.type())] + _index_set.index(e);
      }

      ElementMapperBase(const GV& gv)
        : _gt_offsets(LocalGeometryTypeIndex::size(dim) + 1)
        , _index_set(gv.indexSet())
      {
        update();
      }

    private:

      std::vector<size_type> _gt_offsets;
      const IndexSet& _index_set;

    };

    // implementation for grids with a single codim 0 geometry type
    template<typename GV>
    class ElementMapperBase<GV,true>
    {

    protected:

      typedef typename GV::template Codim<0>::Entity Element;
      typedef typename GV::IndexSet IndexSet;
      typedef std::size_t size_type;

      void update()
      {}

      size_type map(const Element& e) const
      {
        return _index_set.index(e);
      }

      ElementMapperBase(const GV& gv)
        : _index_set(gv.indexSet())
      {}

    private:

      const IndexSet& _index_set;

    };

#endif // DOXYGEN


    //! Class providing a consecutive index for codim 0 entities of a GridView.
    /**
     * ElementMapper yields a unique, consecutive, zero-based index for all
     * codim 0 entities of a GridView. Conceptually, this can also be achieved
     * using a Mapper from dune-grid, but this implementation has better performance
     * as it avoids looking up the GeometryType in a std::map. For the common case
     * of grids with a single codim 0 GeometryType, ElementMapper is specialized
     * to just look up the cell index on the IndexSet of the GridView.
     *
     * \tparam GV  The type of the GridView to operate on.
     */
    template<typename GV>
    class ElementMapper
      : public ElementMapperBase<GV,
                                 Dune::Capabilities::hasSingleGeometryType<
                                   typename GV::Grid
                                   >::v
                                 >
    {

      typedef ElementMapperBase<
        GV,
        Dune::Capabilities::hasSingleGeometryType<
          typename GV::Grid
          >::v
        > BaseT;

    public:

      //! The type of the returned index.
      typedef typename BaseT::size_type size_type;

      //! The type of the codim 0 entities of the GridView.
      typedef typename BaseT::Element Element;

      //! Construct a CellIndexProvider for the given GridView.
      /**
       * \param gv  the GridView to operate on.
       */
      ElementMapper(const GV& gv)
        : BaseT(gv)
      {}

      //! Return the index of the given element.
      /**
       * Returns an index for the given codim 0 entity that is guaranteed to be
       * zero-based, consecutive and unique across all codim 0 entities
       * of the GridView.
       *
       * \param e  The codim 0 entity for which to calculate an index.
       * \return   The index of the entity.
       */
      size_type map(const Element& e) const
      {
        return BaseT::map(e);
      }

    };

  } // namespace PDELab

} // namespace Dune

#endif // DUNE_PDELAB_COMMON_ELEMENTMAPPER_HH
