#ifndef DUNE_PDELAB_BASIS_ISTL_BACKEND_HH
#define DUNE_PDELAB_BASIS_ISTL_BACKEND_HH

#include <dune/istl/bvector.hh>

#include <dune/common/fvector.hh>
#include <dune/common/tuplevector.hh>

#include <concepts>
#include <type_traits>

namespace Dune::PDELab::inline Experimental {

  struct ISTLBaseBackend {
    template<class T>
    static auto makeVector() {
      return Dune::BlockVector<T>{};
    }

    template<class T, std::size_t k>
    static auto makeArray() {
      return Dune::FieldVector<T, k>{};
    }

    template<class... T>
    static auto makeTuple() {
      return Dune::TupleVector<T...>{};
    }

    template<class LeafSpaceNode>
    static auto makeField() {
      static_assert(AlwaysFalse<LeafSpaceNode>{}, "Not Implemented");
    }

    template<class C>
    struct block_type {
      using type = typename C::block_type;
    };

    template<class T0, class... T>
    struct block_type<Dune::TupleVector<T0, T...>> {
      static_assert((std::common_with<T0,T> && ...), "Tuple arguments do not have a common block type");
      using type = std::common_type_t<T0,T...>;
    };
  };

  template<class FieldType>
  struct ISTLUniformBackend : public ISTLBaseBackend {

    template<class LeafSpaceNode>
    static auto makeField() {
      return FieldType{};
    }
  };

  struct ISTLRangeBackend : public ISTLBaseBackend {

    template<class LeafSpaceNode>
    static auto makeField() {
      return typename LeafSpaceNode::Traits::FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType{}; // omg
    }
  };

  // makes istl container with the same field everywhere
  template<class Field>
  static constexpr auto istl_uniform_backend(Field) {return ISTLUniformBackend<Field>{};}

  // makes istl container with doubles everywhere
  static constexpr auto istl_backend = istl_uniform_backend(double{});

  // makes istl container with the range of the finite element
  static constexpr ISTLRangeBackend istl_range_field_backend{};

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_BASIS_ISTL_BACKEND_HH
