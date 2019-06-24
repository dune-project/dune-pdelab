#ifndef DUNE_PDELAB_LOCALOPERATOR_INTROSPECTION_HH
#define DUNE_PDELAB_LOCALOPERATOR_INTROSPECTION_HH

#include <type_traits>
#include <string>
#include <iomanip>
#include <sstream>

#include <dune/pdelab/assembler/context.hh>
#include <dune/pdelab/localoperator/guardedcalls.hh>

namespace Dune::PDELab::Experimental {

  namespace Introspection {

    template<typename GridOperator>
    auto patternContext(const GridOperator& go)
      -> std::add_lvalue_reference_t<
        decltype(
          makeContext(
            go.matrixPatternEngine()->context(
              go.assembler())))
          >;

    template<typename GridOperator>
    auto residualContext(const GridOperator& go)
      -> std::add_lvalue_reference_t<
        decltype(
          makeContext(
            go.residualEngine()->context(
              go.assembler())))
        >;

    template<typename GridOperator>
    auto jacobianContext(const GridOperator& go)
      -> std::add_lvalue_reference_t<
        decltype(
          makeContext(
            go.jacobianEngine()->context(
              go.assembler())))
          >;

    template<typename GridOperator>
    auto applyJacobianContext(const GridOperator& go)
      -> std::add_lvalue_reference_t<
        decltype(
          makeContext(
            go.applyJacobianEngine()->context(
              go.assembler())))
          >;

  }

  template<typename GO>
  class Introspector
  {

  public:

    using GridOperator = GO;
    //using LocalOperator = typename GridOperator::Traits::LocalOperator;

    constexpr Introspector(const GridOperator& go)
      : _go(&go)
    {}

  private:

    static const GridOperator& makeGridOperator();
    static const typename GridOperator::Traits::LocalOperator& makeLocalOperator();

    const GridOperator* _go;

    template<typename F, typename... Args>
    constexpr static auto invocable(F&& f, Args&&... args)
      -> std::is_invocable<F,decltype(std::forward<Args>(args))...>;

  public:

    constexpr const GridOperator& gridOperator() const
    {
      return *_go;
    }

    constexpr const typename GridOperator::Traits::LocalOperator& localOperator() const
    {
      return gridOperator().localOperator();
    }

    constexpr static bool intersectionsTwoSided()
    {
      // put this in unevaluated context to avoid requiring a definition of makeLocalOperator()
      return decltype(LocalOperator::intersectionsTwoSided(makeLocalOperator())){};
    }

    constexpr static bool functionSpaceFlavors()
    {
      // put this in unevaluated context to avoid requiring a definition of makeLocalOperator()
      return not decltype(LocalOperator::disableFunctionSpaceFlavors(makeLocalOperator())){};
    }

    constexpr static bool possiblyNonLinear()
    {
      return models<Concept::PossiblyNonLinear,typename GridOperator::Traits::LocalOperator>();
    }

    constexpr bool nonLinear() const
    {
      return isNonlinear(localOperator());
    }

    constexpr static bool start()
    {
      return decltype(
        invocable(
          LocalOperator::start(),
          makeLocalOperator(),
          Introspection::residualContext(makeGridOperator())
          )){};
    }

    constexpr static bool finish()
    {
      return decltype(
        invocable(
          LocalOperator::finish(),
          makeLocalOperator(),
          Introspection::residualContext(makeGridOperator())
          )){};
    }

    constexpr static bool skipCell()
    {
      return decltype(
        invocable(
          LocalOperator::skipCell(),
          makeLocalOperator(),
          Introspection::residualContext(makeGridOperator()
            ))){};
    }

    constexpr static bool startCell()
    {
      return decltype(
        invocable(
          LocalOperator::startCell(),
          makeLocalOperator(),
          Introspection::residualContext(makeGridOperator()
            ))){};
    }

    constexpr static bool finishCell()
    {
      return decltype(
        invocable(
          LocalOperator::finishCell(),
          makeLocalOperator(),
          Introspection::residualContext(makeGridOperator()
            ))){};
    }

    constexpr static bool startIntersections()
    {
      return decltype(
        invocable(
          LocalOperator::startIntersections(),
          makeLocalOperator(),
          Introspection::residualContext(makeGridOperator()
            ))){};
    }

    constexpr static bool finishIntersections()
    {
      return decltype(
        invocable(
          LocalOperator::finishIntersections(),
          makeLocalOperator(),
          Introspection::residualContext(makeGridOperator()
            ))){};
    }


    constexpr static bool volumePattern()
    {
      return decltype(
        invocable(
          LocalOperator::volumePattern(),
          makeLocalOperator(),
          Introspection::patternContext(makeGridOperator()
            ))){};
    }

    constexpr static bool skeletonPattern()
    {
      return decltype(
        invocable(
          LocalOperator::skeletonPattern(),
          makeLocalOperator(),
          Introspection::patternContext(makeGridOperator()
            ))){};
    }

    constexpr static bool boundaryPattern()
    {
      return decltype(
        invocable(
          LocalOperator::boundaryPattern(),
          makeLocalOperator(),
          Introspection::patternContext(makeGridOperator()
            ))){};
    }

    constexpr static bool volumePatternPostIntersections()
    {
      return decltype(
        invocable(
          LocalOperator::volumePatternPostIntersections(),
          makeLocalOperator(),
          Introspection::patternContext(makeGridOperator()
            ))){};
    }


    constexpr static bool volumeIntegral()
    {
      return decltype(
        invocable(
          LocalOperator::volumeIntegral(),
          makeLocalOperator(),
          Introspection::residualContext(makeGridOperator()
            ))){};
    }

    constexpr static bool skeletonIntegral()
    {
      return decltype(
        invocable(
          LocalOperator::skeletonIntegral(),
          makeLocalOperator(),
          Introspection::residualContext(makeGridOperator()
            ))){};
    }

    constexpr static bool boundaryIntegral()
    {
      return decltype(
        invocable(
          LocalOperator::boundaryIntegral(),
          makeLocalOperator(),
          Introspection::residualContext(makeGridOperator()
            ))){};
    }

    constexpr static bool volumeIntegralPostIntersections()
    {
      return decltype(
        invocable(
          LocalOperator::volumeIntegralPostIntersections(),
          makeLocalOperator(),
          Introspection::residualContext(makeGridOperator()
            ))){};
    }


    constexpr static bool volumeJacobian()
    {
      return decltype(
        invocable(
          LocalOperator::volumeJacobian(),
          makeLocalOperator(),
          Introspection::jacobianContext(makeGridOperator()
            ))){};
    }

    constexpr static bool skeletonJacobian()
    {
      return decltype(
        invocable(
          LocalOperator::skeletonJacobian(),
          makeLocalOperator(),
          Introspection::jacobianContext(makeGridOperator()
            ))){};
    }

    constexpr static bool boundaryJacobian()
    {
      return decltype(
        invocable(
          LocalOperator::boundaryJacobian(),
          makeLocalOperator(),
          Introspection::jacobianContext(makeGridOperator()
            ))){};
    }

    constexpr static bool volumeJacobianPostIntersections()
    {
      return decltype(
        invocable(
          LocalOperator::volumeJacobianPostIntersections(),
          makeLocalOperator(),
          Introspection::jacobianContext(makeGridOperator()
            ))){};
    }


    constexpr static bool volumeApplyJacobian()
    {
      return decltype(
        invocable(
          LocalOperator::volumeApplyJacobian(),
          makeLocalOperator(),
          Introspection::applyJacobianContext(makeGridOperator()
            ))){};
    }

    constexpr static bool skeletonApplyJacobian()
    {
      return decltype(
        invocable(
          LocalOperator::skeletonApplyJacobian(),
          makeLocalOperator(),
          Introspection::applyJacobianContext(makeGridOperator()
            ))){};
    }

    constexpr static bool boundaryApplyJacobian()
    {
      return decltype(
        invocable(
          LocalOperator::boundaryApplyJacobian(),
          makeLocalOperator(),
          Introspection::applyJacobianContext(makeGridOperator()
            ))){};
    }

    constexpr static bool volumeApplyJacobianPostIntersections()
    {
      return decltype(
        invocable(
          LocalOperator::volumeApplyJacobianPostIntersections(),
          makeLocalOperator(),
          Introspection::applyJacobianContext(makeGridOperator()
            ))){};
    }



    friend std::ostream& operator<<(std::ostream& os, const Introspector& inspect)
    {
      std::vector<std::pair<std::string_view,std::string>> info;
      auto stringize = [](const auto& data) {
                         std::stringstream stream;
                         stream << std::boolalpha << data;
                         return stream.str();
                       };

      info.emplace_back("intersectionsTwoSided",stringize(inspect.intersectionsTwoSided()));
      info.emplace_back("functionSpaceFlavors",stringize(inspect.functionSpaceFlavors()));
      info.emplace_back("possiblyNonLinear",stringize(inspect.possiblyNonLinear()));
      info.emplace_back("nonLinear",stringize(inspect.nonLinear()));

      info.emplace_back("","");
      info.emplace_back("start()",stringize(inspect.start()));
      info.emplace_back("finish()",stringize(inspect.finish()));
      info.emplace_back("skipCell()",stringize(inspect.skipCell()));
      info.emplace_back("startCell()",stringize(inspect.startCell()));
      info.emplace_back("finishCell()",stringize(inspect.finishCell()));
      info.emplace_back("startIntersections()",stringize(inspect.startIntersections()));
      info.emplace_back("finishIntersections()",stringize(inspect.finishIntersections()));

      info.emplace_back("","");
      info.emplace_back("volumePattern()",stringize(inspect.volumePattern()));
      info.emplace_back("skeletonPattern()",stringize(inspect.skeletonPattern()));
      info.emplace_back("boundaryPattern()",stringize(inspect.boundaryPattern()));
      info.emplace_back("volumePatternPostIntersections()",stringize(inspect.volumePatternPostIntersections()));

      info.emplace_back("","");
      info.emplace_back("volumeIntegral()",stringize(inspect.volumeIntegral()));
      info.emplace_back("skeletonIntegral()",stringize(inspect.skeletonIntegral()));
      info.emplace_back("boundaryIntegral()",stringize(inspect.boundaryIntegral()));
      info.emplace_back("volumeIntegralPostIntersections()",stringize(inspect.volumeIntegralPostIntersections()));

      info.emplace_back("","");
      info.emplace_back("volumeJacobian()",stringize(inspect.volumeJacobian()));
      info.emplace_back("skeletonJacobian()",stringize(inspect.skeletonJacobian()));
      info.emplace_back("boundaryJacobian()",stringize(inspect.boundaryJacobian()));
      info.emplace_back("volumeJacobianPostIntersections()",stringize(inspect.volumeJacobianPostIntersections()));

      info.emplace_back("","");
      info.emplace_back("volumeApplyJacobian()",stringize(inspect.volumeApplyJacobian()));
      info.emplace_back("skeletonApplyJacobian()",stringize(inspect.skeletonApplyJacobian()));
      info.emplace_back("boundaryApplyJacobian()",stringize(inspect.boundaryApplyJacobian()));
      info.emplace_back("volumeApplyJacobianPostIntersections()",stringize(inspect.volumeApplyJacobianPostIntersections()));

      auto width = std::max_element(
        begin(info),end(info),
        [](auto&& a, auto&& b)
        {
          return a.first.size() < b.first.size();
        })->first.size();

      for (auto&& [name, value] : info)
        if (not name.empty())
          os << std::left << std::setw(width) << name << " : " << value << std::endl;
        else
          os << std::endl;

      return os;
    }

  };


} // namespace Dune::PDELab::Experimental

#endif // DUNE_PDELAB_LOCALOPERATOR_INTROSPECTION_HH
