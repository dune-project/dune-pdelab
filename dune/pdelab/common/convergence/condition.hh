#ifndef DUNE_PDELAB_COMMON_CONVERGENCE_CONDITION_HH
#define DUNE_PDELAB_COMMON_CONVERGENCE_CONDITION_HH

#include <dune/pdelab/common/convergence/reason.hh>
#include <dune/pdelab/common/convergence/properties.hh>
#include <dune/pdelab/common/property_tree.hh>

#include <concepts>
#include <span>
#include <string>

namespace Dune::PDELab::inline Experimental::Convergence {

  /**
   * @brief Interface for classes that implement a convergence condition for a sequence of residuals
   *
   * @tparam Residual  The type to use to check convergence
   */
  template<class Residual>
  class Condition : public PropertyTree {

  public:

    virtual ~Condition() = default;

    /**
     * @brief Evaluates a convergence reason for a sequence of defects
     *
     * @return Reason  Enumerator describing the reason of (not) convergence
     */
    virtual Reason evaluate(std::span<const Residual>) const = 0;
  };


  /**
   * @brief Default implementations for a convergence condition of a sequence of defects
   * @details This class allows iteration based on an iteration range passed by the property tree
   * If an absoulte or relative tolerance is achevied during those iterations, a
   * successful convergence reason is returned.
   *
   * @tparam Residual  The type to use to check convergence
   */
  template<class Residual>
  requires std::regular<Residual> && std::totally_ordered<Residual>
  class DefaultCondition : public Condition<Residual> {
  public:

    /**
     * @brief Construct a new default convergence condition object
     * @details The properties "iteration_range" and "relative_tolerance" are default
     * initialized whereas the field "absolute_tolerance" needs to be manually set up
     * since it's a problem dependent parameter.
     */
    DefaultCondition(std::optional<Residual> abs_tolerance = std::nullopt)
      : Condition<Residual>{}
    {
      PropertyTree& properties = *this;
      properties["iteration_range"] = makeIterationRange();
      properties["relative_tolerance"] = makeRelativeTolerance();
      properties["absolute_tolerance"] = makeAbsoluteTolerance<Residual>(abs_tolerance);
    }

    /**
     * @brief Evaluates a convergence reason for a sequence of defects
     * @details The following set of ordered conditions are evaluated to determine
     * a convergence reason:
     *
     * * If the number of iterations is below the iteration range, the method is Iterating
     * * If an absolute tolerance is set and the last defect is below, the method is ConvergedByAbsoluteTolerance
     * * If the ratio between first and last defect is below relative tolerance, the method is ConvergedByRelativeTolerance
     * * If the number of iterations is above the iteration range, the method is DivergedByIterations
     * * Otherwise, the method is Iterating
     *
     * @note The first defect in the sequence counts as the 0-th iteration and may
     * be used to evaluate the tolerance condition if the iteration range
     * starts from 0. Respectively, an iteration range with a minimum higher than 0
     * forces at least one iteration to happen.
     *
     * @param defects  A sequence of defects, one defect per iteration (including 0-th iteration)
     * @return Reason  Enumerator describing the reason of (not) convergence
     */
    Reason evaluate(std::span<const Residual> defects) const override final {
      const PropertyTree& properties = *this;
      const auto& it_range = properties.get("iteration_range").as_vector();
      std::size_t min_iterations = unwrap_property_ref<const std::size_t>(it_range[0]);
      std::size_t max_iterations = unwrap_property_ref<const std::size_t>(it_range[1]);

      if (defects.size() == 0)
        return Reason::Iterating;

      std::size_t it = defects.size() - 1;
      const auto& first_defect = defects.front();
      const auto& current_defect = defects.back();

      if (it < min_iterations)
        return Reason::Iterating;

      if constexpr (std::numeric_limits<Residual>::has_infinity)
        if (std::isinf(current_defect))
          return Reason::DivergedByNanOrInf;

      if constexpr (std::numeric_limits<Residual>::has_quiet_NaN or std::numeric_limits<Residual>::has_signaling_NaN)
        if (std::isnan(current_defect))
          return Reason::DivergedByNanOrInf;

      const auto& abs_tolerance = properties["absolute_tolerance"];
      if (abs_tolerance.has_value() and current_defect <= unwrap_property_ref<const Residual>(abs_tolerance))
        return Reason::ConvergedByAbsoluteTolerance;
      if (current_defect < Residual{properties.template get<double>("relative_tolerance")*first_defect})
        return Reason::ConvergedByRelativeTolerance;
      if (it > max_iterations)
        return Reason::DivergedByIterations;
      return Reason::Iterating;
    }
  };

} // namespace Dune::PDELab::inline Experimental::Convergence

#endif // DUNE_PDELAB_COMMON_CONVERGENCE_CONDITION_HH
