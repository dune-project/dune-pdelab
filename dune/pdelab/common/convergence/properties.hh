#ifndef DUNE_PDELAB_COMMON_CONVERGENCE_PROPERTIES_HH
#define DUNE_PDELAB_COMMON_CONVERGENCE_PROPERTIES_HH

#include <dune/pdelab/common/convergence/reason.hh>
#include <dune/pdelab/common/property_tree.hh>

#include <dune/common/float_cmp.hh>

#include <concepts>
#include <span>
#include <string>
#include <optional>

namespace Dune::PDELab::inline Experimental::Convergence
{

  static Property makeRelativeTolerance(std::optional<std::array<double,2>> valid_range = std::array{0.,1.}) {
    return Property{
      /*name*/   "relative_tolerance",
      /*doc*/    "Minimum relative tolerance ||r0||/||ri|| to accept convergence",
      /*setter*/ [=](const Property& ppt){
        if (not ppt.has_value()) return;
        [[maybe_unused]] double rel_tolerance = unwrap_property_ref<const double>(ppt);
        if (valid_range) {
          auto [min_val, max_val] = valid_range.value();
          if (Dune::FloatCmp::lt(rel_tolerance, min_val))
            throw ppt.bad_cast_exception(ppt.type(), "Minimum relative tolerance has to be bigger or equal than " + std::to_string(min_val));
          if (Dune::FloatCmp::gt(rel_tolerance, max_val))
            throw ppt.bad_cast_exception(ppt.type(), "Minimum relative tolerance has to be less or equal than " + std::to_string(max_val));
        }
      }} = /*default value*/ 1e-8;
  }

  template<class T>
  static Property makeAbsoluteTolerance(std::optional<T> abs_tolerance = std::nullopt) {
    Property ppt{
      /*name*/   "absolute_tolerance",
      /*doc*/    "Minimum absolute tolerance ||ri|| to accept convergence",
      [](const Property& ppt){
        if (not ppt.has_value()) return;
        [[maybe_unused]] double rel_tolerance = unwrap_property_ref<const double>(ppt);
      }
    };
    if (abs_tolerance)
      ppt = abs_tolerance.value();
    return ppt;
  }


  static Property makeIterationRange() {
    return Property{
      /*name*/   "iteration_range",
      /*doc*/    "Iteration range [min_it, max_it] allowed to assert convergence",
      /*setter*/ [](const Property& ppt){
        if (not ppt.has_value()) return;
        auto& array = ppt.as_vector();
        if (array.size() != 2)
          throw ppt.bad_cast_exception(ppt.type(), "Iteration range most have two entries");
        std::size_t min_its = unwrap_property_ref<const std::size_t>(array[0]);
        std::size_t max_its = unwrap_property_ref<const std::size_t>(array[1]);
        if (min_its >= max_its)
          throw ppt.bad_cast_exception(ppt.type(), "Iteration range is invalid");
      }
    } = /*default value*/ {std::size_t{0}, std::size_t{45}};
  }

} // namespace Dune::PDELab::inline Experimental::Convergence

#endif // DUNE_PDELAB_COMMON_CONVERGENCE_PROPERTIES_HH
