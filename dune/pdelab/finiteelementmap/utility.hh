#ifndef DUNE_PDELAB_FINITEELEMENTMAP_UTILITY_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_UTILITY_HH

#include <cstddef>

#include <dune/common/keywords.hh>
#include <dune/geometry/type.hh>

namespace Dune {
  namespace PDELab {

    //! Metafunction that returns the type of FEM::size() iff that function is static.
    /**
     * This metafunction is mostly for detecting the nature of FEM::size() with if_detected et. al.
     */
    template<typename FEM>
    using StaticFEMSize = decltype(FEM::size(GeometryTypes::vertex));

#ifndef DOXYGEN

    namespace Impl {

      // This function iterates over all geometry types up to the dimension of the finite element map
      // and returns the value of FEM::size(gt) iff that number is constant for all geometry types for
      // which the returned size is > 0. Otherwise it returns 0. As this only works if FEM::size() is
      // static, we use the additional argument to provide a separate overload if it is not.
      // Note that as there is no way to easily construct the set of "valid" geometry types for a
      // given dimension, we manually iterate over all possible topology ids. This creates weird
      // geometry types, but we just assume that FEM::size() will return 0 for invalid ones.
      template<typename FEM>
      constexpr std::size_t _femBlockSize(std::true_type)
      {
        constexpr int dim = FEM::dimension;
        for (int d = 0 ; d <= dim ; ++d)
          {
            std::size_t size = FEM::size(GeometryTypes::none(d));
            if (size > 0)
              return size;
            for (unsigned int topology_id = 0 ; topology_id < (1 << dim) ; ++topology_id)
              {
                std::size_t size = FEM::size(GeometryType(topology_id,d));
                if (size > 0)
                  return size;
              }
          }
        return 0;
      }

      // fallback version if `FEM::size()` is an instance method.
      template<typename FEM>
      constexpr std::size_t _femBlockSize(std::false_type)
      {
        return 0;
      }

      template<typename FEM>
      using Dimension = std::integral_constant<decltype(FEM::dimension),FEM::dimension>;

      template<typename FEM>
      constexpr std::size_t femBlockSizeCheckFEMInterface(std::true_type)
      {
        return _femBlockSize<FEM>(Std::is_detected<StaticFEMSize,FEM>());
      }

      template<typename FEM>
      DUNE_DEPRECATED_MSG("Your finite element map does not export the dimension. After the release of PDELab 2.6, this will cause compilation failures.")
      constexpr std::size_t femBlockSizeCheckFEMInterface(std::false_type)
      {
        return 0;
      }

    } // namespace Impl

#endif // DOXYGEN

    //! Returns the block size for FEM if available, 0 otherwise.
    /**
     * The block size is given by `FEM::size(gt)` iff that function returns a single value for
     * all GeometryTypes `gt` for which `FEM::size(gt) > 0`. This requires that the value is
     * the same for all instances of FEM and that FEM attaches the same number of DOFs to all
     * GeometryTypes to which it attaches DOFs at all.
     *
     * If the above condition does not hold, the function returns 0 instead.
     */
    template<typename FEM>
    constexpr std::size_t finiteElementMapBlockSize()
    {
      return Impl::femBlockSizeCheckFEMInterface<FEM>(Std::is_detected<Impl::Dimension,FEM>());
    }

    //! An alias template that encapsulates the result of `finiteElementMapBlockSize<FEM>()` in an integral constant.
    template<typename FEM>
    using FiniteElementMapBlockSize = std::integral_constant<std::size_t,finiteElementMapBlockSize<FEM>()>;

  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_UTILITY_HH
