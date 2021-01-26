// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_BLOCKING_HH
#define DUNE_PDELAB_BACKEND_ISTL_BLOCKING_HH

#include <tuple>
#include <utility>
#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/version.hh>
#include <dune/common/std/type_traits.hh>

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

namespace Dune {
namespace PDELab {
namespace Blocking {
  // Lightweight representation of (hierarchic) block structure
  namespace tag
  {
    struct Unknown {};

    /// \brief Blocking structure corresponding to IndexMergingStrategies \ref FlatLexicographic and \ref FlatInterleaved
    struct flat
    {
      template <class Index>
      flat operator[](const Index&) const { return {}; }
    };

    /// \brief Blocking structure corresponding to IndexMergingStrategy \ref BlockedInterleaved
    template <std::size_t N>
    struct leafBlocked
    {
      template <class Index>
      flat operator[](const Index&) const { return {}; }
    };

    /// \brief Blocking structure corresponding to IndexMergingStrategy \ref BlockedLexicographic
    template <class Tag0, class... Tags>
    struct blocked
    {
      template <std::size_t i>
      auto operator[](index_constant<i>) const { return std::get<i>(t); }
      Tag0 operator[](std::size_t /*i*/) const { return {}; }

      std::tuple<Tag0,Tags...> t;
    };

  } // end namespace Blocking

  // forward declaration
  template <class GlobalBasis>
  struct Blocking;

  /// \brief Extract blocking structure of GlobalBasis
  template <class GlobalBasis>
  using Blocking_t = typename Blocking<GlobalBasis>::type;


  namespace Concept
  {
    template <class T>
    struct Blocked
        : std::false_type {};

    template <class... Tags>
    struct Blocked<tag::blocked<Tags...>>
        : std::true_type {};

    template <class T>
    constexpr bool isBlocked() { return Blocked<T>::value; }

  } // end namespace Concept


  namespace Impl_
  {
    // Handle non-leaf nodes
    template <template <class...> class Tag, class PreBasis>
    struct CompositeBlocking
    {
      using type = tag::flat;
    };

    // Flat index-merging strategy
    template <class PreBasis, class IndexMergingStrategy>
    struct BlockingImpl
    {
      // Collapse flat tags to flat, if hierarchy is all flat, otherwise
      // to tag::unknown
      template <class... T>
      using FlatTag_t = std::conditional_t<Std::conjunction<std::is_same<tag::flat,T>...>::value,
        tag::flat, tag::Unknown>;

      using type = typename Impl_::CompositeBlocking<FlatTag_t,PreBasis>::type;
    };

    // Block index-merging strategy
    template <class PB>
    struct BlockingImpl<PB, Functions::BasisFactory::BlockedLexicographic>
    {
      using type = typename Impl_::CompositeBlocking<tag::blocked, PB>::type;
    };

    // Leaf-Block index-merging strategy
    template <class PB>
    struct BlockingImpl<PB, Functions::BasisFactory::BlockedInterleaved>
    {
#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,7)
      using RootTreePath = decltype(Dune::TypeTree::hybridTreePath());
      using Node = typename PB::template Node<RootTreePath>;
#else
      using Node = typename PB::Node;
#endif
      using type = tag::leafBlocked<Node::CHILDREN>;
    };

    // Handle leaf-nodes in basis-tree
    template <class PreBasis, class = void>
    struct BlockingPreBasisImpl
    {
      using type = tag::flat;
    };

    // Handle all nodes with an index-merging trategy
    template <class PB>
    struct BlockingPreBasisImpl<PB, void_t<typename PB::IndexMergingStrategy>>
    {
      using type = typename BlockingImpl<PB, typename PB::IndexMergingStrategy>::type;
    };

    /// \brief Extract block structure of PreBasis
    template <class PreBasis>
    using BlockingPreBasis_t = typename Impl_::BlockingPreBasisImpl<PreBasis>::type;

    // defaultbasis
    template <class Basis, class = void>
    struct PreBasisImpl
    {
      using type = typename Basis::PreBasis;
    };

    // subspacebasis
    template <class Basis>
    struct PreBasisImpl<Basis, void_t<typename Basis::RootBasis>>
    {
      using type = typename Basis::RootBasis::PreBasis;
    };

    template <template <class...> class Builder, class MI, class IMS, class... SF>
    struct CompositeBlocking<Builder, Functions::CompositePreBasis<MI, IMS, SF...>>
    {
      using type = Builder<BlockingPreBasis_t<SF>...>;
    };

    template <template <class...> class Builder, class MI, class IMS, class SF, std::size_t C>
    struct CompositeBlocking<Builder, Functions::PowerPreBasis<MI, IMS, SF, C>>
    {
      template <std::size_t I> using SubFactories = SF;
      template <class Idx>        struct Expand;
      template <std::size_t... I> struct Expand<std::index_sequence<I...>>
      {
        using type = Builder<BlockingPreBasis_t<SubFactories<I>>...>;
      };

      using type = typename Expand<std::make_index_sequence<C>>::type;
    };

  } // end namespace Impl_

  template <class GlobalBasis>
  struct Blocking
  {
    using type = Impl_::BlockingPreBasis_t<typename Impl_::PreBasisImpl<GlobalBasis>::type>;
  };

}}} // end namespace Dune::PDELab::Blocking

#endif // DUNE_PDELAB_BACKEND_ISTL_BLOCKING_HH
