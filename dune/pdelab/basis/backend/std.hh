#ifndef DUNE_PDELAB_BASIS_STD_BACKEND_HH
#define DUNE_PDELAB_BASIS_STD_BACKEND_HH

#include <dune/common/typetraits.hh>
#include <dune/common/tuplevector.hh>

#include <vector>
#include <array>
#include <tuple>
#include <concepts>
#include <type_traits>

namespace Dune::PDELab {

  struct STDBaseBackend {
    template<class T>
    static auto makeVector() {
      return std::vector<T>{};
    }

    template<class T, std::size_t k>
    static auto makeArray() {
      return std::array<T, k>{};
    }

    template<class... T>
    static auto makeTuple() {
      // an std tuple with integral constant bracket accessors
      return Dune::TupleVector<T...>{};
    }

    template<class LeafSpaceNode>
    static auto makeField() {
      static_assert(AlwaysFalse<LeafSpaceNode>{}, "Not Implemented");
    }

    template<class C>
    struct block_type {
      using type = typename C::value_type;
    };

    template<class T0, class... T>
    struct block_type<Dune::TupleVector<T0, T...>> {
      static_assert((std::common_with<T0,T> && ...), "Tuple arguments do not have a common block type");
      using type = std::common_type_t<T0,T...>;
    };
  };

  template<class FieldType>
  struct STDUniformBackend : public STDBaseBackend {

    template<class LeafSpaceNode>
    static auto makeField() {
      return FieldType{};
    }

  };


  struct STDRangeBackend : public STDBaseBackend {

    template<class LeafSpaceNode>
    static auto makeField() {
      return typename LeafSpaceNode::Traits::FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType{}; // omg
    }
  };

  // makes std container with the same field everywhere
  template<class Field>
  static constexpr auto std_uniform_backend(Field) {return STDUniformBackend<Field>{};}

  // makes std container with doubles everywhere
  static constexpr auto std_backend = std_uniform_backend(double{});

  // makes std container with the range of the finite element
  static constexpr STDRangeBackend std_range_field_backend{};
}
 // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_BASIS_STD_BACKEND_HH
