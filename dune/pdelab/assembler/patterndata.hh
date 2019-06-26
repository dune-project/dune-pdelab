// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_PATTERNDATA_HH
#define DUNE_PDELAB_ASSEMBLER_PATTERNDATA_HH

#include <cassert>
#include <type_traits>

#include <dune/pdelab/common/intersectiontype.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/assembler/utility.hh>
#include <dune/pdelab/assembler/vectordata.hh>

namespace Dune::PDELab::Experimental {

  template<typename Index_>
  struct PatternLink
  {
    using Index = Index_;

    Index i = 0;
    Index j = 0;

  };

  template<typename Index_>
  class LocalPattern
  {

    using Storage    = std::vector<PatternLink<Index_>>;

  public:

    using Index      = Index_;
    using Link       = typename Storage::value_type;
    using value_type = typename Storage::value_type;
    using size_type  = std::size_t;
    using iterator   = typename Storage::const_iterator;

    iterator begin() const
    {
      return _links.begin();
    }

    iterator end() const
    {
      return _links.end();
    }

    size_type size() const
    {
      return _links.size();
    }

    void clear()
    {
      _links.clear();
    }

    template<typename LFSV, typename LFSU>
    void addLink(const LFSV& lfsv, Index i, const LFSU& lfsu, Index j)
    {
      _links.push_back({Index(lfsv.localIndex(i)),Index(lfsu.localIndex(j))});
    }

  private:

    Storage _links;

  };

  template<typename Context>
  class ElementPatternData
    : public Context
  {

  public:

    using Context::unbind;
    using Pattern = LocalPattern<std::uint32_t>;

    Context* unbind(const typename Context::Entity&, typename Context::Index, typename Context::Index)
    {
      Context::engine().scatterPattern(*this,_pattern,Context::test().cache(),Context::trial().cache());
      _pattern.clear();
      return this;
    }

    Pattern& pattern()
    {
      return _pattern;
    }

    ElementPatternData(Context&& context)
      : Context(std::move(context))
    {}

  private:

    Pattern _pattern;

  };

  template<typename Context>
  auto elementPatternData(Context&& context)
  {
    return ElementPatternData<Context>(std::move(context));
  }


  template<typename Context>
  class SkeletonPatternData
    : public Context
  {

  public:

    SkeletonPatternData(Context&& context)
      : Context(std::move(context))
    {}

    using Pattern = LocalPattern<std::uint32_t>;

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
          Context::engine().scatterPattern(*this,_pattern_io,Context::inside().test().cache(),Context::outside().trial().cache());
          Context::engine().scatterPattern(*this,_pattern_oi,Context::outside().test().cache(),Context::inside().trial().cache());

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

} // namespace Dune::PDELab::Experimental

#endif // DUNE_PDELAB_ASSEMBLER_PATTERNDATA_HH
