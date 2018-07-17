// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_PATTERNDATA_HH
#define DUNE_PDELAB_ASSEMBLER_PATTERNDATA_HH

#include <cassert>
#include <type_traits>

#include <dune/pdelab/common/intersectiontype.hh>
#include <dune/pdelab/assembler/utility.hh>
#include <dune/pdelab/assembler/vectordata.hh>

namespace Dune {
  namespace PDELab {

    template<typename Context>
    class CellPatternData
      : public Context
    {

    public:

      using Context::unbind;
      using Pattern = LocalSparsityPattern;

      Context* unbind(const typename Context::Entity&, typename Context::Index, typename Context::Index)
      {
        Context::engine().scatterPattern(_pattern,Context::test().cache(),Context::trial().cache());
        _pattern.clear();
        return this;
      }

      Pattern& pattern()
      {
        return _pattern;
      }

      CellPatternData(Context&& context)
        : Context(std::move(context))
      {}

    private:

      Pattern _pattern;

    };

    template<typename Context>
    auto cellPatternData(Context&& context)
    {
      return CellPatternData<Context>(std::move(context));
    }


    template<typename Context>
    class SkeletonPatternData
      : public Context
    {

    public:

      SkeletonPatternData(Context&& context)
        : Context(std::move(context))
      {}

      using Pattern = LocalSparsityPattern;

      Pattern _pattern_io;
      Pattern _pattern_oi;

      using Context::pattern;

      template<typename Inside, typename Outside>
      std::enable_if_t<
        std::is_same_v<Inside,typename Context::Inside> and std::is_same_v<Outside,typename Context::Outside>,
        Pattern&
        >
      pattern(const Inside&, const Outside&)
      {
        return _pattern_io;
      }

      template<typename Outside, typename Inside>
      std::enable_if_t<
        std::is_same_v<Outside,typename Context::Outside> and std::is_same_v<Inside,typename Context::Inside>,
        Pattern&
        >
      pattern(const Outside&, const Inside&)
      {
        return _pattern_oi;
      }

      using Context::unbind;

      template<typename IntersectionType_>
      Context* unbind(
        IntersectionType_ type,
        const typename Context::IntersectionDomain::Intersection&,
        typename Context::IntersectionDomain::Index,
        const typename Context::IntersectionDomain::Entity& entity,
        typename Context::IntersectionDomain::Index entity_index,
        typename Context::IntersectionDomain::Index unique_index
        )
      {
        if (type == IntersectionType::skeleton or type == IntersectionType::periodic)
          {
            Context::engine().scatterPattern(_pattern_io,Context::inside().test().cache(),Context::outside().trial().cache());
            Context::engine().scatterPattern(_pattern_oi,Context::outside().test().cache(),Context::inside().trial().cache());

            _pattern_io.clear();
            _pattern_oi.clear();
          }
        return this;
      }

    };

    template<typename Context>
    auto skeletonPatternData(Context&& context)
    {
      return SkeletonPatternData<Context>(std::move(context));
    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_PATTERNDATA_HH
