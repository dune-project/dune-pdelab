// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_ELECTRODYNAMIC_HH
#define DUNE_PDELAB_LOCALOPERATOR_ELECTRODYNAMIC_HH

#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/deprecated.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/pdelab/localoperator/numericalresidual.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    namespace ElectrodynamicImpl {

      //////////////////////////////////////////////////////////////////////
      // Deprecated interface handling

      // This is an adaptor that turns the old-style policy-object based
      // parameters into a function-object based one.
      //
      // This is here for backward-compatibility only.
      template<class Params>
      class Functor
      {
        const Params *params_;

      public:
        DUNE_DEPRECATED_MSG("You are using Dune::PDELab::Electrodynamic_T or "
                            "Dune::PDELab::Electrodynamic_S with an old-style "
                            "parameter class for eps/mu.  Please use a "
                            "callable (function object, (generic) lambda, or "
                            "a function pointer) instead.")
        Functor(const Params &params) : params_(&params) {}

        template<class Element, class Pos>
        typename Params::Traits::RangeType
        operator()(const Element &e, const Pos &xl) const
        {
          typename Params::Traits::RangeType result;
          params_->evaluate(e, xl, result);
          return result;
        };
      };

      template<class Params, class = void>
      struct IsOldstyleParams : std::false_type {};

      template<class Params>
      struct IsOldstyleParams
          <Params,
           std::enable_if_t<sizeof(typename Params::Traits::RangeType)> > :
        std::true_type
      {};

      template<class Params>
      using ConstRefOrFunctor =
        std::conditional_t<IsOldstyleParams<Params>::value,
                           Functor<Params>, Params>;

      //////////////////////////////////////////////////////////////////////
      // Curl manipulation

      // dimension of the curl for a given space dimension
      constexpr std::size_t dimOfCurl(std::size_t dimOfSpace)
      {
        return
          dimOfSpace == 1 ? 2 :
          dimOfSpace == 2 ? 1 :
          dimOfSpace == 3 ? 3 :
          // Dune exceptions are difficult to construct in constexpr functions
          throw std::invalid_argument("Only applicable for dimensions 1-3");
      }

      template<typename RF>
      void jacobianToCurl(FieldVector<RF, 1> &curl,
                          const FieldMatrix<RF, 2, 2> &jacobian)
      {
        curl[0] = jacobian[1][0] - jacobian[0][1];
      }
      template<typename RF>
      void jacobianToCurl(FieldVector<RF, 3> &curl,
                          const FieldMatrix<RF, 3, 3> &jacobian)
      {
        for(unsigned i = 0; i < 3; ++i)
          curl[i] = jacobian[(i+2)%3][(i+1)%3] - jacobian[(i+1)%3][(i+2)%3];
      }

    } // namespace ElectrodynamicImpl

    //! Construct matrix T for the Electrodynamic operator
    /**
     * Construct the matrix
     * \f[
     *    T_{ij}=\int_\Omega\epsilon\mathbf N_i\cdot\mathbf N_jdV
     * \f]
     * which appears inside the Electrodynamic operator.
     *
     * \tparam Eps    Type of function to evaluate \f$\epsilon\f$
     */
    template<typename Eps>
    class Electrodynamic_T
      : public FullVolumePattern
      , public LocalOperatorDefaultFlags
      , public JacobianBasedAlphaVolume<Electrodynamic_T<Eps> >
    {
      static constexpr bool oldstyle =
        ElectrodynamicImpl::IsOldstyleParams<Eps>::value;

    public:

      // pattern assembly flags
      static constexpr bool doPatternVolume = true;
      static constexpr bool doAlphaVolume = true;

      //! Construct an Electrodynamic_T localoperator
      /**
       * \param eps    (Reference to) Function object to evaluate.
       * \param qorder Quadrature order to use.
       *
       * \note When using the deprecated old style parameter interface, `eps`
       *       needed to be a reference to a parameter object, and the local
       *       operator would store a reference to that object.  With the new
       *       style interface, `eps` is a callable (function object, lambda,
       *       or even a plain function pointer), and the local operator
       *       stores a copy.  To get the old behaviour in a non-deprecated
       *       way, wrap the callable in a `std::reference_wrapper`.
       */
      Electrodynamic_T(const Eps &eps, int qorder = 2) :
        eps_(eps), qorder_(qorder)
      {}

      Electrodynamic_T(Eps &&eps, int qorder = 2) :
        eps_(std::move(eps)), qorder_(qorder)
      {}

      /**
       * \note We support only Galerkin method lfsu==lfsv
       */
      template<typename EG, typename LFS, typename X, typename M>
      void jacobian_volume (const EG& eg, const LFS& lfsu, const X& x,
                            const LFS& lfsv, M& mat) const
      {
        using BasisTraits =
          typename LFS::Traits::FiniteElementType::Traits::Basis::Traits;

        // static checks
        static constexpr unsigned dimR = BasisTraits::dimRange;
        static_assert(dimR == 3 || dimR == 2, "Works only in 2D or 3D");

        using ctype = typename EG::Geometry::ctype;
        using DF = typename BasisTraits::DomainField;
        static_assert(std::is_same<ctype, DF>::value, "Grids ctype and "
                      "Finite Elements DomainFieldType must match");

        using Range = typename BasisTraits::Range;
        std::vector<Range> phi(lfsu.size());

        // loop over quadrature points
        for(const auto &qp : quadratureRule(eg.geometry(), qorder_))
        {
          // values of basefunctions
          lfsu.finiteElement().basis().evaluateFunction(qp.position(),phi);

          // calculate T
          auto factor = qp.weight()
            * eg.geometry().integrationElement(qp.position())
            * eps_(eg.entity(), qp.position());

          for(unsigned i = 0; i < lfsu.size(); ++i)
            for(unsigned j = 0; j < lfsu.size(); ++j)
              mat.accumulate(lfsv,i,lfsu,j,factor * (phi[i] * phi[j]));
        }
      }

    private:
      ElectrodynamicImpl::ConstRefOrFunctor<Eps> eps_;
      int qorder_;
    };

    //! construct an Electrodynamic_T operator
    /**
     * This relieves the user from the need to construct the type of `eps`.
     */
    template<class Eps>
    Electrodynamic_T<std::decay_t<Eps> >
    makeLocalOperatorEdynT(Eps &&eps, int qorder = 2)
    {
      return { std::forward<Eps>(eps), qorder };
    }

    //! Contruct matrix S for the Electrodynamic operator
    /**
     * Construct the matrix
     * \f[
     *    S_{ij}=\int_\Omega\mu^{-1}(\nabla\times\mathbf N_i)\cdot
     *                              (\nabla\times\mathbf N_j)dV
     * \f]
     * which appears inside the Electrodynamic operator.
     *
     * \tparam Mu    Type of function to evaluate \f$\mu\f$
     */
    template<typename Mu>
    class Electrodynamic_S
      : public FullVolumePattern
      , public LocalOperatorDefaultFlags
      , public JacobianBasedAlphaVolume<Electrodynamic_S<Mu> >
    {
      static constexpr bool oldstyle =
        ElectrodynamicImpl::IsOldstyleParams<Mu>::value;

    public:

      // pattern assembly flags
      static constexpr bool doPatternVolume = true;
      static constexpr bool doAlphaVolume = true;

      //! Construct an Electrodynamic_S localoperator
      /**
       * \param mu     (Reference to) Function object to evaluate
       * \param qorder Quadrature order to use.
       *
       * \note When using the deprecated old style parameter interface, `mu`
       *       needed to be a reference to a parameter object, and the local
       *       operator would store a reference to that object.  With the new
       *       style interface, `mu` is a callable (function object, lambda,
       *       or even a plain function pointer), and the local operator
       *       stores a copy.  To get the old behaviour in a non-deprecated
       *       way, wrap the callable in a `std::reference_wrapper`.
       */
      Electrodynamic_S(const Mu &mu, int qorder = 2) :
        mu_(mu), qorder_(qorder)
      {}

      Electrodynamic_S(Mu &&mu, int qorder = 2) :
        mu_(std::move(mu)), qorder_(qorder)
      {}

      /**
       * \note We support only Galerkin method lfsu==lfsv
       */
      template<typename EG, typename LFS, typename X, typename M>
      void jacobian_volume (const EG& eg, const LFS& lfsu, const X& x,
                            const LFS& lfsv, M& mat) const
      {
        using ElectrodynamicImpl::dimOfCurl;
        using ElectrodynamicImpl::jacobianToCurl;

        using BasisTraits =
          typename LFS::Traits::FiniteElementType::Traits::Basis::Traits;

        // static checks
        static constexpr unsigned dimR = BasisTraits::dimRange;
        static_assert(dimR == 3 || dimR == 2, "Works only in 2D or 3D");

        using ctype = typename EG::Geometry::ctype;
        using DF = typename BasisTraits::DomainField;
        static_assert(std::is_same<ctype, DF>::value, "Grids ctype and "
                      "Finite Elements DomainFieldType must match");

        using Jacobian = typename BasisTraits::Jacobian;
        std::vector<Jacobian> J(lfsu.size());

        using RF = typename BasisTraits::RangeField;
        using Curl = FieldVector<RF, dimOfCurl(dimR)>;
        std::vector<Curl> rotphi(lfsu.size());

        // loop over quadrature points
        for(const auto &qp : quadratureRule(eg.geometry(), qorder_))
        {
          // curl of the basefunctions
          lfsu.finiteElement().basis().evaluateJacobian(qp.position(),J);
          for(unsigned i = 0; i < lfsu.size(); ++i)
            jacobianToCurl(rotphi[i], J[i]);

          // calculate S
          auto factor = qp.weight()
            * eg.geometry().integrationElement(qp.position())
            / mu_(eg.entity(), qp.position());

          for(unsigned i = 0; i < lfsu.size(); ++i)
            for(unsigned j = 0; j < lfsu.size(); ++j)
              mat.accumulate(lfsv,i,lfsu,j,factor * (rotphi[i] * rotphi[j]));

        }
      }

    private:
      ElectrodynamicImpl::ConstRefOrFunctor<Mu> mu_;
      int qorder_;
    };

    //! construct an Electrodynamic_S operator
    /**
     * This relieves the user from the need to construct the type of `mu`.
     */
    template<class Mu>
    Electrodynamic_S<std::decay_t<Mu> >
    makeLocalOperatorEdynS(Mu &&mu, int qorder = 2)
    {
      return { std::forward<Mu>(mu), qorder };
    }

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_ELECTRODYNAMIC_HH
