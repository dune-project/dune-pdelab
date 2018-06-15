// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ASSEMBLER_ASSEMBLER_HH
#define DUNE_PDELAB_ASSEMBLER_ASSEMBLER_HH

namespace Dune {
  namespace PDELab {

    template<typename ES>
    class Assembler
    {

    public:

      using EntitySet = ES;

      const EntitySet& entitySet() const
      {
        return _es;
      }

      Assembler(const EntitySet& es)
        : _es(es)
      {}

      template<typename Engine>
      decltype(auto) assemble(Engine& engine) const
      {
        constexpr bool visit_periodic_intersections = engine.visitPeriodicIntersections();
        constexpr bool visit_skeleton_intersections = engine.visitSkeletonIntersections();
        constexpr bool visit_boundary_intersections = engine.visitBoundaryIntersections();
        constexpr bool visit_processor_intersections = engine.visitProcessorIntersections();
        constexpr bool visit_intersections =
          visit_periodic_intersections or
          visit_skeleton_intersections or
          visit_boundary_intersections or
          visit_processor_intersections;
        constexpr bool intersections_two_sided [[maybe_unused]] = engine.intersectionsTwoSided();

        auto &entity_set = _es;
        auto &index_set = _es.indexSet();

        typename EntitySet::template Codim<0>::Entity invalid_element;
        constexpr auto invalid_index [[maybe_unused]] = index_set.invalidIndex();

        auto ctx = engine.context(*this);
        ctx.setup(engine);

        engine.start(ctx);

        for (const auto& element : elements(entity_set))
          {

            auto entity_index = index_set.index(element);
            auto unique_index = index_set.uniqueIndex(element);

            if (not engine.skipCell(ctx,element,entity_index))
              {

                engine.startCell(ctx,element,entity_index);

                ctx.bind(element,entity_index,unique_index);

                engine.volume(ctx);

                if constexpr (visit_intersections)
                  {
                    engine.startIntersections(ctx);

                    std::size_t intersection_index = 0;

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
                                bool visit_face = intersections_two_sided or unique_index < unique_idn or engine.skipCell(ctx,outside_element,entity_idn);

                                if (visit_face)
                                  {
                                    ctx.bind(
                                      std::integral_constant<IntersectionType,IntersectionType::skeleton>{},
                                      intersection,
                                      intersection_index,
                                      outside_element,
                                      entity_idn,
                                      unique_idn
                                      );

                                    engine.skeleton(ctx);

                                    ctx.unbind(
                                      std::integral_constant<IntersectionType,IntersectionType::skeleton>{},
                                      intersection,
                                      intersection_index,
                                      outside_element,
                                      entity_idn,
                                      unique_idn
                                      );
                                  }
                              }
                            break;
                          case IntersectionType::periodic:
                            if (visit_periodic_intersections)
                              {
                                auto entity_idn = index_set.index(outside_element);
                                auto unique_idn = index_set.uniqueIndex(outside_element);
                                bool visit_face = intersections_two_sided or unique_index < unique_idn or engine.skipCell(ctx,outside_element,entity_idn);

                                if (visit_face)
                                  {
                                    ctx.bind(
                                      std::integral_constant<IntersectionType,IntersectionType::periodic>{},
                                      intersection,
                                      intersection_index,
                                      outside_element,
                                      entity_idn,
                                      unique_idn
                                      );

                                    engine.periodic(ctx);

                                    ctx.unbind(
                                      std::integral_constant<IntersectionType,IntersectionType::periodic>{},
                                      intersection,
                                      intersection_index,
                                      outside_element,
                                      entity_idn,
                                      unique_idn
                                      );
                                  }
                              }
                            break;
                          case IntersectionType::boundary:
                            if (visit_boundary_intersections)
                              {
                                ctx.bind(
                                  std::integral_constant<IntersectionType,IntersectionType::boundary>{},
                                  intersection,
                                  intersection_index,
                                  invalid_element,
                                  invalid_index,
                                  invalid_index
                                  );

                                engine.boundary(ctx);

                                ctx.unbind(
                                  std::integral_constant<IntersectionType,IntersectionType::boundary>{},
                                  intersection,
                                  intersection_index,
                                  invalid_element,
                                  invalid_index,
                                  invalid_index
                                  );
                              }
                            break;
                          case IntersectionType::processor:
                            if (visit_processor_intersections)
                              {
                                ctx.bind(
                                  std::integral_constant<IntersectionType,IntersectionType::processor>{},
                                  intersection,
                                  intersection_index,
                                  invalid_element,
                                  invalid_index,
                                  invalid_index
                                  );

                                engine.processor(ctx);

                                ctx.unbind(
                                  std::integral_constant<IntersectionType,IntersectionType::processor>{},
                                  intersection,
                                  intersection_index,
                                  invalid_element,
                                  invalid_index,
                                  invalid_index
                                  );
                              }
                            break;
                          default:
                            DUNE_THROW(Dune::Exception,"Encountered invalid intersection type during assembly");
                          }
                        ++intersection_index;
                      }

                    engine.finishIntersections(ctx);

                  }

                engine.volumePostIntersections(ctx);

                ctx.unbind(element,entity_index,unique_index);

                engine.finishCell(ctx,element,entity_index);
              }
          }
        engine.finish(ctx);

        return engine.result(ctx);
      }

    private:

      EntitySet _es;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_ASSEMBLER_HH
