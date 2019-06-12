// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ASSEMBLER_BATCHEDASSEMBLER_HH
#define DUNE_PDELAB_ASSEMBLER_BATCHEDASSEMBLER_HH

#include <algorithm>
#include <array>
#include <tuple>
#include <type_traits>

#include <dune/common/iteratorrange.hh>

#include <dune/pdelab/common/intersectiontype.hh>
#include <dune/pdelab/logging/logger.hh>

namespace Dune::PDELab {

  enum class BatchType
  {
    cell,
    boundary,
    processor,
    skeleton,
    periodic
  };

  class Batch
  {

    BatchType _type;
    int _variant;

  public:

    BatchType type() const
    {
      return _type;
    }

    int variant() const
    {
      return _variant;
    }

    void setVariant(int variant)
    {
      _variant = variant;
    }

  protected:

    Batch(BatchType type, int variant)
      : _type(type)
      , _variant(variant)
    {}

  };

  namespace detail {

    struct NoMask
    {
      static constexpr bool active()
      {
        return false;
      }

      template<typename Batch>
      static Batch& bind(Batch& batch)
      {
        return batch;
      }

      template<typename Batch>
      static Batch& unbind(Batch& batch)
      {
        return batch;
      }

    };

    template<typename T, std::size_t size>
    class PaddedRange
    {
    public:

      template<typename Range>
      PaddedRange(Range& range, std::size_t payload_size)
      {
        using std::min;
        for (std::size_t i = 0 ; i < size ; ++i)
          _range[i] = &range[min(i,payload_size-1)];
      }

      struct iterator
      {

        T& operator*()
        {
          return **_it;
        }

        iterator& operator++()
        {
          ++_it;
          return *this;
        }

        bool operator==(const iterator& r)
        {
          return _it == r._it;
        }

        bool operator!=(const iterator& r)
        {
          return _it != r._it;
        }

        typename std::array<T*,size>::iterator _it;
      };

      iterator begin()
      {
        return {begin(_range)};
      }

      iterator end()
      {
        return {end(_range)};
      }

    private:
      std::array<T*,size> _range;
    };

    class PostfixMask
    {
    public:

      PostfixMask(int size)
        : _size(size)
      {}

      static constexpr bool active()
      {
        return true;
      }

      template<typename Batch>
      auto& bind(Batch& batch)
      {
        return PaddedRange{batch,_size};
      }

      template<typename Batch>
      auto& unbind(Batch& batch)
      {
        return IteratorRange(begin(batch),begin(batch) + _size);
      }

    private:

      int _size;

    };

    template<typename ES, int size_>
    class CellBatch
      : public Batch
    {
    public:

      using EntitySet = ES;
      using Cell      = typename ES::template Codim<0>::Entity;
      using Index     = typename ES::Traits::Index;

      CellBatch(int variant = 0)
        : Batch(BatchType::cell,variant)
      {}

      static constexpr int size()
      {
        return size_;
      }

      const Cell& cell(int i) const
      {
        return std::get<0>(_cells[i]);
      }

      void store(int i, const EntitySet& es, const Cell& cell)
      {
        _cells[i] = std::make_tuple(cell,es.indexSet().index(cell),es.indexSet().uniqueIndex(cell));
      }

      Index index(const EntitySet& es, int i) const
      {
        return std::get<1>(_cells[i]);
      }

      Index uniqueIndex(int i) const
      {
        return std::get<2>(_cells[i]);
      }

      void reset()
      {
        for (auto& cell : _cells)
          cell = std::make_tuple(Cell(),ES::IndexSet::invalidIndex(),ES::IndexSet::invalidIndex());
      }

      std::array<std::tuple<Cell,Index,Index>,size_> _cells;

    };

    template<typename ES, int size_, BatchType type_>
    class IntersectionBatch
      : public Batch
    {
    public:

      using EntitySet    = ES;
      using Intersection = typename ES::Intersection;
      using Index        = typename ES::Traits::Index;

      IntersectionBatch(int variant = 0)
        : Batch(type_,variant)
      {}

      static constexpr int size()
      {
        return size_;
      }

      void store(int i, const EntitySet& es, const Intersection& intersection, int index)
      {
        _intersections[i] = std::make_tuple(intersection,index);
      }

      const Intersection& intersection(int i) const
      {
        return std::get<0>(_intersections[i]);
      }

      Index index(int i) const
      {
        return std::get<1>(_intersections[i]);
      }

      Index insideIndex(const EntitySet& es, int i) const
      {
        return es.index(_intersections[i].inside());
      }

      Index outsideIndex(const EntitySet& es, int i) const
      {
        return es.index(_intersections[i].outside());
      }

      void reset()
      {
        for (auto& intersection : _intersections)
          intersection = std::make_tuple(Intersection(),ES::IndexSet::invalidIndex());
      }

      std::array<std::tuple<Intersection,Index>,size_> _intersections;

    };

  }

  template<typename ES, typename BC>
  class BatchedAssembler
  {

  public:

    using EntitySet   = ES;
    using BatchConfig = BC;

    const EntitySet& entitySet() const
    {
      return _es;
    }

    static constexpr int batchSize()
    {
      return BatchConfig::batchSize();
    }

    const BatchConfig& batchConfig() const
    {
      return *_bc;
    }

    using CellBatch      = typename BatchConfig::CellBatch;
    using BoundaryBatch  = typename BatchConfig::BoundaryBatch;
    using ProcessorBatch = typename BatchConfig::ProcessorBatch;
    using SkeletonBatch  = typename BatchConfig::SkeletonBatch;
    using PeriodicBatch  = typename BatchConfig::PeriodicBatch;


    BatchedAssembler(const EntitySet& es, const BC& batch_controller, const Logger& log)
      : _es(es)
      , _bc(&batch_controller)
      , _log(log)
    {
      createBatches();
    }

    template<typename Engine>
    decltype(auto) assemble(Engine&& engine, std::enable_if_t<std::is_same_v<Engine,std::decay_t<Engine>>,int> = 0) const
    {
      return assemble(engine);
    }


    void createBatches()
    {

      _log.info("Starting batch creation, batchSize={}"_fmt,batchSize());

      constexpr bool visit_periodic_intersections = BC::visitPeriodicIntersections();
      constexpr bool visit_skeleton_intersections = BC::visitSkeletonIntersections();
      constexpr bool visit_boundary_intersections = BC::visitBoundaryIntersections();
      constexpr bool visit_processor_intersections = BC::visitProcessorIntersections();
      constexpr bool visit_intersections =
                  visit_periodic_intersections or
                  visit_skeleton_intersections or
                  visit_boundary_intersections or
                  visit_processor_intersections;
      constexpr bool intersections_two_sided [[maybe_unused]] = BC::intersectionsTwoSided();

      using Index = typename EntitySet::IndexSet::Index;

      auto& entity_set = _es;
      auto& index_set = _es.indexSet();
      auto& batch_config = batchConfig();

      _batches.clear();
      _partial_batches.clear();
      for (auto& variant : _cell_batches)
        variant.clear();
      for (auto& variant : _boundary_batches)
        variant.clear();
      for (auto& variant : _processor_batches)
        variant.clear();
      for (auto& variant : _skeleton_batches)
        variant.clear();
      for (auto& variant : _periodic_batches)
        variant.clear();

      {
        int i = 0;
        for (auto& variant : _tail_cell_batches)
        {
          for (auto& [batch,size] : variant)
          {
            batch.setVariant(i);
            size = 0;
          }
          ++i;
        }

        i = 0;
        for (auto& variant : _tail_boundary_batches)
        {
          for (auto& [batch,size] : variant)
          {
            batch.setVariant(i);
            size = 0;
          }
          ++i;
        }

        i = 0;
        for (auto& variant : _tail_processor_batches)
        {
          for (auto& [batch,size] : variant)
          {
            batch.setVariant(i);
            size = 0;
          }
          ++i;
        }

        i = 0;
        for (auto& variant : _tail_skeleton_batches)
        {
          for (auto& [batch,size] : variant)
          {
            batch.setVariant(i);
            size = 0;
          }
          ++i;
        }

        i = 0;
        for (auto& variant : _tail_periodic_batches)
        {
          for (auto& [batch,size] : variant)
          {
            batch.setVariant(i);
            size = 0;
          }
          ++i;
        }

      }

      std::vector<std::pair<BatchType,int>> batch_order;
      int cell_batch_count = 0;
      int boundary_batch_count = 0;
      int processor_batch_count = 0;
      int skeleton_batch_count = 0;
      int periodic_batch_count = 0;

      for (const auto& element : elements(entity_set))
      {

        auto entity_index = index_set.index(element);
        auto unique_index = index_set.uniqueIndex(element);

        _log.trace("At element {}"_fmt,entity_index);

        if (not batch_config.skipCell(element,entity_index))
        {
          bool batched = false;
          int variant = batch_config.volumeVariant(element,entity_index,unique_index);
          auto& colors = _tail_cell_batches[variant];
          int color = 0;
          for (auto& [batch,size] : colors)
          {
            if (batch_config.volumeCompatible(element,entity_index,batch,size))
            {
              batch.store(size++,entity_set,element);
              batched = true;
              if (size == batchSize())
              {
                _cell_batches[variant].push_back(batch);
                batch_order.emplace_back(BatchType::cell,variant);
                batch.reset();
                size = 0;
                ++cell_batch_count;
                _log.trace("Completed cell batch: variant={}, color={}"_fmt,variant,color);
              }
              break;
            }
            ++color;
          }

          if (not batched)
            DUNE_THROW(Exception,"insufficient number of cell colors: " << BatchConfig::maxColors());

          if constexpr (visit_intersections)
          {
            Index intersection_index = 0;

            for (const auto& intersection : intersections(_es,element))
            {
              auto intersection_data = classifyIntersection(entity_set,intersection);
              auto intersection_type = std::get<0>(intersection_data);
              auto& outside_element = std::get<1>(intersection_data);

              switch(intersection_type)
              {
              case IntersectionType::skeleton:
                if constexpr (visit_skeleton_intersections)
                {
                  auto entity_idn = index_set.index(outside_element);
                  auto unique_idn = index_set.uniqueIndex(outside_element);
                  // The final condition in here makes sure that the intersection will be visited even if the engine decides to skip
                  // the outside element
                  bool visit_face = intersections_two_sided or unique_index < unique_idn or batch_config.skipCell(outside_element,entity_idn);

                  if (visit_face)
                  {

                    int variant = batch_config.skeletonVariant(intersection,intersection_index,element,entity_index,unique_index,outside_element,entity_idn,unique_idn);
                    auto& colors = _tail_skeleton_batches[variant];
                    bool batched = false;
                    int color = 0;
                    for (auto& [batch,size] : colors)
                    {
                      if (batch_config.skeletonCompatible(intersection,intersection_index,element,entity_index,unique_index,outside_element,entity_idn,unique_idn,batch,size))
                      {
                        batch.store(size++,entity_set,intersection,intersection_index);
                        batched = true;
                        if (size == batchSize())
                        {
                          _skeleton_batches[variant].push_back(std::move(batch));
                          batch_order.emplace_back(BatchType::skeleton,variant);
                          batch.reset();
                          size = 0;
                          ++skeleton_batch_count;
                          _log.trace("Completed skeleton batch: variant={}, color={}"_fmt,variant,color);
                        }
                        break;
                      }
                      ++color;
                    }

                    if (not batched)
                      DUNE_THROW(Exception,"insufficient number of skeleton intersection colors: " << BatchConfig::maxColors());
                  }
                }
                break;
              case IntersectionType::periodic:
                if (visit_periodic_intersections)
                {
                  auto entity_idn = index_set.index(outside_element);
                  auto unique_idn = index_set.uniqueIndex(outside_element);
                  bool visit_face = intersections_two_sided or unique_index < unique_idn or batch_config.skipCell(outside_element,entity_idn);

                  if (visit_face)
                  {

                    int variant = batch_config.periodicVariant(intersection,intersection_index,element,entity_index,unique_index,outside_element,entity_idn,unique_idn);
                    auto& colors = _tail_periodic_batches[variant];
                    bool batched = false;
                    int color = 0;
                    for (auto& [batch,size] : colors)
                    {
                      if (batch_config.periodicCompatible(intersection,intersection_index,element,entity_index,unique_index,outside_element,entity_idn,unique_idn,batch,size))
                      {
                        batch.store(size++,entity_set,intersection,intersection_index);
                        batched = true;
                        if (size == batchSize())
                        {
                          _periodic_batches[variant].push_back(std::move(batch));
                          batch_order.emplace_back(BatchType::periodic,variant);
                          batch.reset();
                          size = 0;
                          ++periodic_batch_count;
                          _log.trace("Completed periodic batch: variant={}, color={}"_fmt,variant,color);
                        }
                        break;
                      }
                      ++color;
                    }

                    if (not batched)
                      DUNE_THROW(Exception,"insufficient number of periodic intersection colors: " << BatchConfig::maxColors());
                  }
                }
                break;
              case IntersectionType::boundary:
                if (visit_boundary_intersections)
                {

                  int variant = batch_config.boundaryVariant(intersection,intersection_index,element,entity_index,unique_index);
                  auto& colors = _tail_boundary_batches[variant];
                  bool batched = false;
                  int color = 0;
                  for (auto& [batch,size] : colors)
                  {
                    if (batch_config.boundaryCompatible(intersection,intersection_index,element,entity_index,unique_index,batch,size))
                    {
                      batch.store(size++,entity_set,intersection,intersection_index);
                      batched = true;
                      if (size == batchSize())
                      {
                        _boundary_batches[variant].push_back(std::move(batch));
                        batch_order.emplace_back(BatchType::boundary,variant);
                        batch.reset();
                        size = 0;
                        ++boundary_batch_count;
                        _log.trace("Completed boundary batch: variant={}, color={}"_fmt,variant,color);
                      }
                      break;
                    }
                    ++color;
                  }

                  if (not batched)
                    DUNE_THROW(Exception,"insufficient number of boundary intersection colors: " << BatchConfig::maxColors());

                }
                break;
              case IntersectionType::processor:
                if (visit_processor_intersections)
                {

                  int variant = batch_config.processorVariant(intersection,intersection_index,element,entity_index,unique_index);
                  auto& colors = _tail_processor_batches[variant];
                  bool batched = false;
                  int color;
                  for (auto& [batch,size] : colors)
                  {
                    if (batch_config.processorCompatible(intersection,intersection_index,element,entity_index,unique_index,batch,size))
                    {
                      batch.store(size++,entity_set,intersection,intersection_index);
                      batched = true;
                      if (size == batchSize())
                      {
                        _processor_batches[variant].push_back(std::move(batch));
                        batch_order.emplace_back(BatchType::processor,variant);
                        batch.reset();
                        size = 0;
                        ++processor_batch_count;
                        _log.trace("Completed processor batch: variant={}, color={}"_fmt,variant,color);
                      }
                      break;
                    }
                    ++color;
                  }

                  if (not batched)
                    DUNE_THROW(Exception,"insufficient number of processor intersection colors: " << BatchConfig::maxColors());

                }
                break;
              default:
                DUNE_THROW(Dune::Exception,"Encountered invalid intersection type during assembly");
              }
              ++intersection_index;
            }

          }

        }
      }

      {
        std::array<typename std::vector<CellBatch>::iterator,BatchConfig::volumeVariants()> cell_batches;
        std::array<typename std::vector<BoundaryBatch>::iterator,BatchConfig::boundaryVariants()> boundary_batches;
        std::array<typename std::vector<ProcessorBatch>::iterator,BatchConfig::processorVariants()> processor_batches;
        std::array<typename std::vector<SkeletonBatch>::iterator,BatchConfig::skeletonVariants()> skeleton_batches;
        std::array<typename std::vector<PeriodicBatch>::iterator,BatchConfig::periodicVariants()> periodic_batches;

        std::transform(begin(_cell_batches),end(_cell_batches),begin(cell_batches),[](auto& v){ return begin(v);});
        std::transform(begin(_boundary_batches),end(_boundary_batches),begin(boundary_batches),[](auto& v){ return begin(v);});
        std::transform(begin(_processor_batches),end(_processor_batches),begin(processor_batches),[](auto& v){ return begin(v);});
        std::transform(begin(_skeleton_batches),end(_skeleton_batches),begin(skeleton_batches),[](auto& v){ return begin(v);});
        std::transform(begin(_periodic_batches),end(_periodic_batches),begin(periodic_batches),[](auto& v){ return begin(v);});

        for (auto& [type,variant] : batch_order)
          switch (type)
          {
          case BatchType::cell:
            _batches.push_back(&(*cell_batches[variant]++));
            break;
          case BatchType::boundary:
            _batches.push_back(&(*boundary_batches[variant]++));
            break;
          case BatchType::processor:
            _batches.push_back(&(*processor_batches[variant]++));
            break;
          case BatchType::skeleton:
            _batches.push_back(&(*skeleton_batches[variant]++));
            break;
          case BatchType::periodic:
            _batches.push_back(&(*periodic_batches[variant]++));
            break;
          }
      }

      int variant = 0;
      for (auto& colors : _tail_cell_batches)
      {
        int color = 0;
        for (auto& [batch,size] : colors)
        {
          if (size > 0)
          {
            _partial_batches.emplace_back(&batch,size);
            _log.trace("Partial cell batch: variant={}, color={}, size={}"_fmt,variant,color,size);
          }
          ++color;
        }
        ++variant;
      }

      variant = 0;
      for (auto& colors : _tail_boundary_batches)
      {
        int color = 0;
        for (auto& [batch,size] : colors)
        {
          if (size > 0)
          {
            _partial_batches.emplace_back(&batch,size);
            _log.trace("Partial boundary batch: variant={}, color={}, size={}"_fmt,variant,color,size);
          }
          ++color;
        }
        ++variant;
      }

      variant = 0;
      for (auto& colors : _tail_processor_batches)
      {
        int color = 0;
        for (auto& [batch,size] : colors)
        {
          if (size > 0)
          {
            _partial_batches.emplace_back(&batch,size);
            _log.trace("Partial processor batch: variant={}, color={}, size={}"_fmt,variant,color,size);
          }
          ++color;
        }
        ++variant;
      }

      variant = 0;
      for (auto& colors : _tail_skeleton_batches)
      {
        int color = 0;
        for (auto& [batch,size] : colors)
        {
          if (size > 0)
          {
            _partial_batches.emplace_back(&batch,size);
            _log.trace("Partial skeleton batch: variant={}, color={}, size={}"_fmt,variant,color,size);
          }
          ++color;
        }
        ++variant;
      }

      variant = 0;
      for (auto& colors : _tail_periodic_batches)
      {
        int color = 0;
        for (auto& [batch,size] : colors)
        {
          if (size > 0)
          {
            _partial_batches.emplace_back(&batch,size);
            _log.trace("Partial periodic batch: variant={}, color={}, size={}"_fmt,variant,color,size);
          }
          ++color;
        }
        ++variant;
      }

      _log.info("Completed batch creation: batches={}, partial_batches={}"_fmt,_batches.size(),_partial_batches.size());
      _log.detail(
        "cell_batches={}, boundary_batches={}, processor_batches={}, skeleton_batches={}, periodic_batches={}"_fmt,
        cell_batch_count,
        boundary_batch_count,
        processor_batch_count,
        skeleton_batch_count,
        periodic_batch_count
        );

    }

    using Element = typename ES::template Codim<0>::Entity;
    using Index   = typename ES::Traits::Index;

    struct Cache
    {
      std::array<Element,batchSize()> inside_elements;
      std::array<Index,batchSize()> inside_indices;
      std::array<Index,batchSize()> inside_unique_indices;
      std::array<Element,batchSize()> outside_elements;
      std::array<Index,batchSize()> outside_indices;
      std::array<Index,batchSize()> outside_unique_indices;
      Element invalid_element;
      Index invalid_index = ES::IndexSet::invalidIndex();
    };

    template<typename IndexSet, typename Engine, typename Context, typename Mask = detail::NoMask>
    void visit(CellBatch& batch, const IndexSet& index_set, Engine& engine, Context& ctx, Cache& cache, Mask mask = Mask{})
    {
      int i = 0;
      for (auto& [element,entity_index,unique_index] : mask.bind(batch))
      {
        ctx[i++].bind(element,entity_index,unique_index);
      }

      if constexpr (Mask::active())
        engine.volume(ctx,mask());
      else
        engine.volume(ctx);

      i = 0;
      for (auto& [element,entity_index,unique_index] : mask.unbind(batch))
      {
        ctx[i++].unbind(element,entity_index,unique_index);
      }
    }

    template<typename IndexSet, typename Engine, typename Context, typename Mask = detail::NoMask>
    void visit(BoundaryBatch& batch, const IndexSet& index_set, Engine& engine, Context& ctx, Cache& cache, Mask mask = Mask{})
    {
      int i = 0;
      for (auto& [intersection,intersection_index] : mask.bind(batch))
      {
        cache.inside_elements[i] = intersection.inside();
        cache.inside_indices[i] = index_set.index(cache.inside_elements[i]);
        cache.inside_unique_indices[i] = index_set.uniqueIndex(cache.inside_elements[i]);

        ctx[i].bind(cache.inside_elements[i],cache.inside_indices[i],cache.inside_unique_indices[i]);

        ctx[i++].bind(
          std::integral_constant<IntersectionType,IntersectionType::boundary>{},
          intersection,
          intersection_index,
          cache.invalid_element,
          cache.invalid_index,
          cache.invalid_index
          );
      }

      if constexpr (Mask::active())
        engine.boundary(ctx,mask());
      else
        engine.boundary(ctx);

      i = 0;
      for (auto& [intersection,intersection_index] : mask.unbind(batch))
      {
        ctx.unbind(
          std::integral_constant<IntersectionType,IntersectionType::boundary>{},
          intersection,
          intersection_index,
          cache.invalid_element,
          cache.invalid_index,
          cache.invalid_index
          );

        ctx[i].unbind(cache.inside_elements[i],cache.inside_indices[i],cache.inside_unique_indices[i]);
        ++i;
      }
    }

    template<typename IndexSet, typename Engine, typename Context, typename Mask = detail::NoMask>
    void visit(ProcessorBatch& batch, const IndexSet& index_set, Engine& engine, Context& ctx, Cache& cache, Mask mask = Mask{})
    {
      int i = 0;
      for (auto& [intersection,intersection_index] : mask.bind(batch))
      {
        cache.inside_elements[i] = intersection.inside();
        cache.inside_indices[i] = index_set.index(cache.inside_elements[i]);
        cache.inside_unique_indices[i] = index_set.uniqueIndex(cache.inside_elements[i]);

        ctx[i].bind(cache.inside_elements[i],cache.inside_indices[i],cache.inside_unique_indices[i]);

        ctx[i++].bind(
          std::integral_constant<IntersectionType,IntersectionType::processor>{},
          intersection,
          intersection_index,
          cache.invalid_element,
          cache.invalid_index,
          cache.invalid_index
          );
      }

      if constexpr (Mask::active())
        engine.processor(ctx,mask());
      else
        engine.processor(ctx);

      i = 0;
      for (auto& [intersection,intersection_index] : mask.unbind(batch))
      {
        ctx.unbind(
          std::integral_constant<IntersectionType,IntersectionType::processor>{},
          intersection,
          intersection_index,
          cache.invalid_element,
          cache.invalid_index,
          cache.invalid_index
          );

        ctx[i].unbind(cache.inside_elements[i],cache.inside_indices[i],cache.inside_unique_indices[i]);
        ++i;
      }
    }

    template<typename IndexSet, typename Engine, typename Context, typename Mask = detail::NoMask>
    void visit(SkeletonBatch& batch, const IndexSet& index_set, Engine& engine, Context& ctx, Cache& cache, Mask mask = Mask{})
    {
      int i = 0;
      for (auto& [intersection,intersection_index] : mask.bind(batch))
      {
        cache.inside_elements[i] = intersection.inside();
        cache.inside_indices[i] = index_set.index(cache.inside_elements[i]);
        cache.inside_unique_indices[i] = index_set.uniqueIndex(cache.inside_elements[i]);

        cache.outside_elements[i] = intersection.outside();
        cache.outside_indices[i] = index_set.index(cache.outside_elements[i]);
        cache.outside_unique_indices[i] = index_set.uniqueIndex(cache.outside_elements[i]);

        ctx[i].bind(cache.inside_elements[i],cache.inside_indices[i],cache.inside_unique_indices[i]);

        ctx[i].bind(
          std::integral_constant<IntersectionType,IntersectionType::skeleton>{},
          intersection,
          intersection_index,
          cache.outside_elements[i],
          cache.outside_indices[i],
          cache.outside_unique_indices[i]
          );

        ++i;
      }

      if constexpr (Mask::active())
        engine.skeleton(ctx,mask());
      else
        engine.skeleton(ctx);

      i = 0;
      for (auto& [intersection,intersection_index] : mask.unbind(batch))
      {

        ctx[i].unbind(
          std::integral_constant<IntersectionType,IntersectionType::skeleton>{},
          intersection,
          intersection_index,
          cache.outside_elements[i],
          cache.outside_indices[i],
          cache.outside_unique_indices[i]
          );

        ctx[i].unbind(cache.inside_elements[i],cache.inside_indices[i],cache.inside_unique_indices[i]);
        ++i;
      }
    }

    template<typename IndexSet, typename Engine, typename Context, typename Mask = detail::NoMask>
    void visit(PeriodicBatch& batch, const IndexSet& index_set, Engine& engine, Context& ctx, Cache& cache, Mask mask = Mask{})
    {
      int i = 0;
      for (auto& [intersection,intersection_index] : mask.bind(batch))
      {
        cache.inside_elements[i] = intersection.inside();
        cache.inside_indices[i] = index_set.index(cache.inside_elements[i]);
        cache.inside_unique_indices[i] = index_set.uniqueIndex(cache.inside_elements[i]);

        cache.outside_elements[i] = intersection.outside();
        cache.outside_indices[i] = index_set.index(cache.outside_elements[i]);
        cache.outside_unique_indices[i] = index_set.uniqueIndex(cache.outside_elements[i]);

        ctx[i].bind(cache.inside_elements[i],cache.inside_indices[i],cache.inside_unique_indices[i]);

        ctx[i].bind(
          std::integral_constant<IntersectionType,IntersectionType::periodic>{},
          intersection,
          intersection_index,
          cache.outside_elements[i],
          cache.outside_indices[i],
          cache.outside_unique_indices[i]
          );

        ++i;
      }

      if constexpr (Mask::active())
        engine.periodic(ctx,mask());
      else
        engine.periodic(ctx);

      i = 0;
      for (auto& [intersection,intersection_index] : mask.unbind(batch))
      {

        ctx[i].unbind(
          std::integral_constant<IntersectionType,IntersectionType::periodic>{},
          intersection,
          intersection_index,
          cache.outside_elements[i],
          cache.outside_indices[i],
          cache.outside_unique_indices[i]
          );

        ctx[i].unbind(cache.inside_elements[i],cache.inside_indices[i],cache.inside_unique_indices[i]);
        ++i;
      }
    }

    template<typename Engine>
    decltype(auto) assemble(Engine& engine) const
    {
      constexpr bool visit_periodic_intersections = Engine::visitPeriodicIntersections();
      constexpr bool visit_skeleton_intersections = Engine::visitSkeletonIntersections();
      constexpr bool visit_boundary_intersections = Engine::visitBoundaryIntersections();
      constexpr bool visit_processor_intersections = Engine::visitProcessorIntersections();
      constexpr bool visit_intersections =
      visit_periodic_intersections or
      visit_skeleton_intersections or
      visit_boundary_intersections or
      visit_processor_intersections;
      constexpr bool intersections_two_sided [[maybe_unused]] = Engine::intersectionsTwoSided();

      using Index = typename EntitySet::IndexSet::Index;

      auto &entity_set = _es;
      auto &index_set = _es.indexSet();

      auto ctx = makeBatchedContext([&](){makeContext(engine.context(*this));},batchConfig());
      ctx.setup();

      engine.start(ctx);
      Cache cache;

      for (auto batch_ptr : _batches)
      {
        switch(batch_ptr->type())
        {
        case BatchType::cell:
          visit(static_cast<CellBatch&>(*batch_ptr),index_set,engine,ctx,cache);
          break;
        case BatchType::boundary:
          visit(static_cast<BoundaryBatch&>(*batch_ptr),index_set,engine,ctx,cache);
          break;
        case BatchType::processor:
          visit(static_cast<ProcessorBatch&>(*batch_ptr),index_set,engine,ctx,cache);
          break;
        case BatchType::skeleton:
          visit(static_cast<SkeletonBatch&>(*batch_ptr),index_set,engine,ctx,cache);
          break;
        case BatchType::periodic:
          visit(static_cast<PeriodicBatch&>(*batch_ptr),index_set,engine,ctx,cache);
          break;
        default:
          std::abort();
        }
      }

      for (auto& [batch_ptr,size] : _partial_batches)
      {
        switch(batch_ptr->type())
        {
        case BatchType::cell:
          visit(static_cast<CellBatch&>(*batch_ptr),index_set,engine,ctx,cache,detail::PostfixMask{size});
          break;
        case BatchType::boundary:
          visit(static_cast<BoundaryBatch&>(*batch_ptr),index_set,engine,ctx,cache,detail::PostfixMask{size});
          break;
        case BatchType::processor:
          visit(static_cast<ProcessorBatch&>(*batch_ptr),index_set,engine,ctx,cache,detail::PostfixMask{size});
          break;
        case BatchType::skeleton:
          visit(static_cast<SkeletonBatch&>(*batch_ptr),index_set,engine,ctx,cache,detail::PostfixMask{size});
          break;
        case BatchType::periodic:
          visit(static_cast<PeriodicBatch&>(*batch_ptr),index_set,engine,ctx,cache,detail::PostfixMask{size});
          break;
        default:
          std::abort();
        }
      }

      engine.finish(ctx);

      return engine.result(ctx);
    }

  private:

    EntitySet _es;
    const BatchConfig* _bc;
    Logger _log;

    std::vector<Batch*> _batches;
    std::vector<std::pair<Batch*,int>> _partial_batches;

    std::array<std::vector<CellBatch>,BatchConfig::volumeVariants()> _cell_batches;
    std::array<std::vector<BoundaryBatch>,BatchConfig::boundaryVariants()> _boundary_batches;
    std::array<std::vector<ProcessorBatch>,BatchConfig::processorVariants()> _processor_batches;
    std::array<std::vector<SkeletonBatch>,BatchConfig::skeletonVariants()> _skeleton_batches;
    std::array<std::vector<PeriodicBatch>,BatchConfig::periodicVariants()> _periodic_batches;

    std::array<std::array<std::pair<CellBatch,int>,BatchConfig::maxColors()>,BatchConfig::volumeVariants()> _tail_cell_batches;
    std::array<std::array<std::pair<BoundaryBatch,int>,BatchConfig::maxColors()>,BatchConfig::boundaryVariants()> _tail_boundary_batches;
    std::array<std::array<std::pair<ProcessorBatch,int>,BatchConfig::maxColors()>,BatchConfig::processorVariants()> _tail_processor_batches;
    std::array<std::array<std::pair<SkeletonBatch,int>,BatchConfig::maxColors()>,BatchConfig::skeletonVariants()> _tail_skeleton_batches;
    std::array<std::array<std::pair<PeriodicBatch,int>,BatchConfig::maxColors()>,BatchConfig::periodicVariants()> _tail_periodic_batches;

  };

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_ASSEMBLER_BATCHEDASSEMBLER_HH
