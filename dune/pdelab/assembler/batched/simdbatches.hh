// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ASSEMBLER_BATCHED_SIMDBATCHES_HH
#define DUNE_PDELAB_ASSEMBLER_BATCHED_SIMDBATCHES_HH

#include <vendor/vectorclass/vectorclass.h>

#include <dune/pdelab/assembler/batchedassembler.hh>

namespace Dune::PDELab {

  namespace SIMD {

    template<std::size_t dim>
    constexpr std::array<int,2> make_corner_indices(index_constant<dim>,std::true_type)
    {
      return {0, (1 << dim) - 1};
    }

    template<std::size_t dim>
    constexpr std::array<int,(1<<dim)> make_corner_indices(index_constant<dim>,std::false_type)
    {
      std::array<int,(1<<dim)> r;
      std::iota(begin(r),end(r),0);
      return r;
    }

    template<std::size_t dim, bool axis_parallel>
    const inline std::array<int,(axis_parallel ? 2 : 1 << dim)> corner_indices = make_corner_indices(index_constant<dim>{},std::bool_constant<axis_parallel>{});

    template<typename ES, int size_, bool axis_parallel>
    class CellBatch
      : public Batch
    {
    public:

      using EntitySet       = ES;
      using CoordinateField = typename ES::ctype;
      using Cell      = typename ES::template Codim<0>::Entity;
      using Index     = typename ES::Traits::Index;

      CellBatch(int variant = 0)
        : Batch(BatchType::cell,variant)
      {
        reset();
      }

      void store(int i, const EntitySet& es, const Cell& cell)
      {
        indices[i] = es.indexSet().index(cell);
        auto geo = cell.geometry();
        int k = 0;
        for (auto j : corner_indices<ES::dimension,axis_parallel>)
        {
          auto c = geo.corner(j);
          for (std::size_t d = 0 ; d < ES::dimension ; ++d)
            corners[j][d][i] = c[d];
          ++k;
        }
      }

      Index index(const EntitySet& es, int i) const
      {
        return indices[i];
      }

      static constexpr int size()
      {
        return size_;
      }

      void reset()
      {
        std::fill(begin(indices),end(indices),ES::IndexSet::invalidIndex());
      }

      std::array<Index,size_> indices;
      std::array<std::array<std::array<CoordinateField,size_>,ES::dimension>,axis_parallel ? 2 : (1 << ES::dimension)> corners;

    };

    template<typename ES, int size_, BatchType type_, bool axis_parallel>
    class IntersectionBatch
      : public Batch
    {
    public:

      using EntitySet    = ES;
      using CoordinateField = typename ES::ctype;
      using Intersection = typename ES::Intersection;
      using Index        = typename ES::Traits::Index;

      IntersectionBatch(int variant = 0)
        : Batch(type_,variant)
      {
        reset();
      }

      static constexpr int size()
      {
        return size_;
      }

      template<typename ES2>
      std::enable_if_t<std::is_same_v<ES2,EntitySet> and (type_ == BatchType::skeleton or type_ == BatchType::periodic)>
      store(int i, const ES2& es, const Intersection& intersection, Index intersection_index)
      {
        auto inside  = intersection.inside();
        auto outside = intersection.outside();
        inside_indices[i]  = es.indexSet().index(inside);
        outside_indices[i] = es.indexSet().index(outside);
        auto inside_geo  = inside.geometry();
        auto outside_geo = outside.geometry();
        int k = 0;
        for (auto j : corner_indices<ES::dimension,axis_parallel>)
        {
          auto inside_c  = inside_geo.corner(j);
          auto outside_c = outside_geo.corner(j);
          for (std::size_t d = 0 ; d < ES::dimension ; ++d)
          {
            inside_corners[j][d][i]  = inside_c[d];
            outside_corners[j][d][i] = outside_c[d];
          }
          ++k;
        }
      }

      template<typename ES2>
      std::enable_if_t<std::is_same_v<ES2,EntitySet> and (type_ == BatchType::boundary or type_ == BatchType::processor)>
      store(int i, const ES2& es, const Intersection& intersection, Index intersection_index)
      {
        auto inside = intersection.inside();
        inside_indices[i]  = es.indexSet().index(inside);
        outside_indices[i] = ES::IndexSet::invalidIndex();
        auto inside_geo = inside.geometry();
        int k = 0;
        for (auto j : corner_indices<ES::dimension,axis_parallel>)
        {
          auto inside_c = inside_geo.corner(j);
          for (std::size_t d = 0 ; d < ES::dimension ; ++d)
            inside_corners[j][d][i] = inside_c[d];
          ++k;
        }
      }

      Index insideIndex(const EntitySet& es, int i) const
      {
        return inside_indices[i];
      }

      Index outsideIndex(const EntitySet& es, int i) const
      {
        return outside_indices[i];
      }

      void reset()
      {
        std::fill(begin(inside_indices),end(inside_indices),ES::IndexSet::invalidIndex());
        std::fill(begin(outside_indices),end(outside_indices),ES::IndexSet::invalidIndex());
      }

      std::array<Index,size_> inside_indices;
      std::array<Index,size_> outside_indices;

      std::array<std::array<std::array<CoordinateField,size_>,ES::dimension>,axis_parallel ? 2 : (1 << ES::dimension)> inside_corners;
      std::array<std::array<std::array<CoordinateField,size_>,ES::dimension>,axis_parallel ? 2 : (1 << ES::dimension)> outside_corners;

    };

  }

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_ASSEMBLER_BATCHED_SIMDBATCHES_HH
