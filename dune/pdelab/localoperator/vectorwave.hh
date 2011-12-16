// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_LOCALOPERATOR_VECTORWAVE_HH
#define DUNE_PDELAB_LOCALOPERATOR_VECTORWAVE_HH

#include <algorithm>
#include <cstddef>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/static_assert.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/pdelab/common/jacobiantocurl.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/multistep/cache.hh>

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    namespace VectorWave {

      //! Parameter class interface for the vector wave local operators
      /**
       * This can also be used as a CRTP base class to implement evaluation by
       * element and local coordinates when the actual function is given in
       * global coordinates.
       *
       * \tparam Imp   The class derived from this class.  That class should
       *               either implement the global-coordinate functions or
       *               overwrite the local-coordinate-plus-element functions.
       * \tparam GV    Type of GridView to operate on.  Used to extract ctype
       *               and dimension.
       * \tparam RF    Field type of the values.
       * \tparam Time_ Type of temporal values.
       */
      template<class Imp, class GV, class RF = double, class Time_ = double>
      struct Parameters {
        //! export type of temporal values
        typedef Time_ Time;

        //! export dimension (both domain and range)
        static const std::size_t dimension = GV::dimension;

        //! field type of domain
        typedef typename GV::ctype DomainField;
        //! vector type of domain
        typedef FieldVector<DomainField, dimension> Domain;

        //! field type of range
        typedef RF RangeField;
        //! vector type of range
        typedef FieldVector<RangeField, dimension> Range;

        //! element type of GridView
        typedef typename GV::template Codim<0>::Entity Element;

        //! \brief evaluate dielectric permittivity
        //!        \f$\varepsilon=\varepsilon_0\varepsilon_r\f$
        /**
         * This implementation forwards calls to the derived class, which
         * should have an accessible member function
         * \code
RangeField epsilon(const Domain& xg) const;
         * \endcode
         */
        RangeField epsilon(const Element& e, const Domain& xl) const
        { return asImp().epsilonGlobal(e.geometry().global(xl)); }
        //! Whether epsilon changes between the given time steps
        /**
         * This is used by the policy class to determine when evaluation
         * values can be reused.  This default implemention always returns
         * true (change occured) to be on the safe side.
         */
        bool epsilonChanged(Time t1, Time t2) const { return true; }
        //! evaluate magnetic permeability \f$\mu=\mu_0\mu_r\f$
        /**
         * This implementation forwards calls to the derived class, which
         * should have an accessible member function
         * \code
RangeField mu(const Domain& xg) const;
         * \endcode
         */
        RangeField mu(const Element& e, const Domain& xl) const
        { return asImp().muGlobal(e.geometry().global(xl)); }
        //! Whether mu changes between the given time steps
        /**
         * This is used by the policy class to determine when evaluation
         * values can be reused.  This default implemention always returns
         * true (change occured) to be on the safe side.
         */
        bool muChanged(Time t1, Time t2) const { return true; }

        //! set the time for subsequent evaluation
        /**
         * This is a no-op in the default implementation.
         */
        void setTime(const Time &time) { }

      private:
        const Imp &asImp() const { return *static_cast<const Imp*>(this); }
      };

      //! Homogenous parameter class for the vector wave local operators
      /**
       * Parameter class with spatially constant \f$\varepsilon\f$ and
       * \f$\mu\f$.
       *
       * \tparam GV    Type of GridView to operate on.  Used to extract ctype
       *               and dimension.
       * \tparam RF    Field type of the values.
       * \tparam Time_ Type of temporal values.
       */
      template<class GV, class RF = double, class Time = double>
      class ConstantParameters :
        public Parameters<ConstantParameters<GV, RF, Time>, GV, RF, Time>
      {
        RF epsilon_;
        RF mu_;

      public:
        ConstantParameters(RF epsilon, RF mu) : epsilon_(epsilon), mu_(mu) { }

        template<typename Domain>
        RF epsilonGlobal(const Domain &) const { return epsilon_; }
        bool epsilonChanged(Time t1, Time t2) const { return false; }
        template<typename Domain>
        RF muGlobal(const Domain &) const { return mu_; }
        bool muChanged(Time t1, Time t2) const { return false; }
      };

      //! \brief Local operator for the vector wave problem,
      //!        no-temporal-derivatives part
      /**
       * The vector-wave equation in its simplest form:
       * \f[
       *    \partial_t^2(\varepsilon\mathbf E)
       *    +\nabla\times\mu^{-1}\nabla\times\mathbf E = 0
       * \f]
       *
       * This local operator implements the part without temporal derivatives
       * \f[
       *    \nabla\times\mu^{-1}\nabla\times\mathbf E
       * \f]
       * which boils down to the matrix
       * \f[
       *    S_{ij}=\int_\Omega\mu^{-1}
       *        (\nabla\times\varphi_j)\cdot(\nabla\times\psi_i)\,dV
       * \f]
       * where \f$\cdot\f$ denotes the scalar product and \f$\nabla\times\f$
       * the curl operator.
       *
       * \tparam Params Type of parameter class providing the values for
       *                \f$\mu\f$.  Should conform to the interface of
       *                Parameters.
       */
      template<class Params>
      class R0 :
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename Params::Time>,
        public JacobianBasedAlphaVolume<R0<Params> >
      {
        typedef InstationaryLocalOperatorDefaultMethods<typename Params::Time>
          IBase;

        Params &params;
        std::size_t qorder;

      public:
        // pattern assembly flags
        enum { doPatternVolume = true };
        enum { doAlphaVolume = true };

        //! Construct a local operator object
        /**
         * \param params_ Parameter object providing values for \f$\mu\f$.
         * \param qorder_ Quadrature order to use.
         *
         * \note The reference to the parameter objects should be valid for as
         *       long as this localoperators is evaluated.
         */
        R0(Params &params_, std::size_t qorder_ = 0) :
          params(params_), qorder(qorder_)
        {}

        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename M>
        void jacobian_volume(const EG& eg, const LFSU& lfsu, const X& x,
                             const LFSV& lfsv, M& mat) const
        {
          // domain and range field type
          typedef typename LFSU::Traits::FiniteElementType FEU;
          typedef typename LFSV::Traits::FiniteElementType FEV;

          const FEU &feU = lfsu.finiteElement();
          const FEV &feV = lfsv.finiteElement();

          typedef typename FEU::Traits::Basis BasisU;
          typedef typename FEV::Traits::Basis BasisV;

          const BasisU &basisU = feU.basis();
          const BasisV &basisV = feV.basis();

          typedef typename Params::RangeField RF;

          typedef typename BasisU::Traits::DomainField DF;
          static const std::size_t dimDLocal = BasisU::Traits::dimDomainLocal;

          typedef typename BasisU::Traits::Jacobian JacobianU;
          typedef typename BasisV::Traits::Jacobian JacobianV;

          typedef JacobianToCurl<JacobianU> J2CU;
          typedef JacobianToCurl<JacobianV> J2CV;

          static const J2CU &j2CU = J2CU();
          static const J2CV &j2CV = J2CV();

          typedef typename J2CU::Curl CurlU;
          typedef typename J2CV::Curl CurlV;

          dune_static_assert(J2CU::dimCurl == J2CV::dimCurl, "Curl dimension "
                             "of ansatz and test functions must match in "
                             "VectorWave::R0");

          // select quadrature rule
          typedef QuadratureRules<DF,dimDLocal> QRs;
          typedef QuadratureRule<DF,dimDLocal> QR;
          typedef typename QR::const_iterator QRIterator;
          GeometryType gt = eg.geometry().type();
          const QR& rule = QRs::rule(gt,qorder);

          std::vector<JacobianU> jacobianU(lfsu.size());
          std::vector<JacobianV> jacobianV(lfsv.size());

          std::vector<CurlU> curlU(lfsu.size());
          std::vector<CurlV> curlV(lfsv.size());

          // loop over quadrature points
          for(QRIterator it=rule.begin(); it != rule.end(); ++it) {
            // curls of basefunctions
            basisU.evaluateJacobian(it->position(), jacobianU);
            basisV.evaluateJacobian(it->position(), jacobianV);

            for(std::size_t j = 0; j < lfsu.size(); ++j)
              j2CU(jacobianU[j], curlU[j]);
            for(std::size_t i = 0; i < lfsv.size(); ++i)
              j2CV(jacobianV[i], curlV[i]);

            RF factor = it->weight()
              * eg.geometry().integrationElement(it->position())
              / params.mu(eg.entity(), it->position());

            for(std::size_t j = 0; j < lfsu.size(); ++j)
              for(std::size_t i = 0; i < lfsv.size(); ++i)
                mat(i,j) += factor * (curlU[j] * curlV[i]);
          }
        }

        //! set time on the parameter object
        void setTime(typename Params::Time time) {
          params.setTime(time);
          IBase::setTime(time);
        }
      };

      //! \brief Local operator for the vector wave problem,
      //!        one-temporal-derivative part
      /**
       * The vector-wave equation in its simplest form:
       * \f[
       *    \partial_t^2(\varepsilon\mathbf E)
       *    +\nabla\times\mu^{-1}\nabla\times\mathbf E = 0
       * \f]
       *
       * This local operator implements the part with one temporal derivatives
       * \f[
       *    0
       * \f]
       * Yes there is no part with only one temporal derivate in the above
       * equation, so this local operator is a dummy which implements nothing.
       *
       * \tparam Params Type of parameter class as for the other local
       *                operators.  For consistency mostly.
       */
      template<class Params>
      class R1 :
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename Params::Time>
      {
      public:
        //! Construct a local operator object
        /**
         * \param params_ Parameter object providing values for \f$\mu\f$.
         * \param qorder_ Quadrature order to use.
         *
         * \note This constructor does nothing and is present for consistency
         *       only.
         */
        R1(Params &params_, std::size_t qorder_ = 0) {}
      };

      //! \brief Local operator for the vector wave problem,
      //!        second-temporal-derivatives part
      /**
       * The vector-wave equation in its simplest form:
       * \f[
       *    \partial_t^2(\varepsilon\mathbf E)
       *    +\nabla\times\mu^{-1}\nabla\times\mathbf E = 0
       * \f]
       *
       * This local operator implements the part with two temporal derivatives
       * \f[
       *    \partial_t^2(\varepsilon\mathbf E)
       * \f]
       * which boils down to the matrix
       * \f[
       *    T_{ij}=\int_\Omega\varepsilon\varphi_j\cdot\psi_i\,dV
       * \f]
       * where \f$\cdot\f$ denotes the scalar product.
       *
       * \tparam Params Type of parameter class providing the values for
       *                \f$\varepsilon\f$.  Should conform to the interface of
       *                Parameters.
       */
      template<class Params>
      class R2 :
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename Params::Time>,
        public JacobianBasedAlphaVolume<R2<Params> >
      {
        typedef InstationaryLocalOperatorDefaultMethods<typename Params::Time>
          IBase;

        Params &params;
        std::size_t qorder;

      public:
        // pattern assembly flags
        enum { doPatternVolume = true };
        enum { doAlphaVolume = true };

        //! Construct a local operator object
        /**
         * \param params_ Parameter object providing values for \f$\mu\f$.
         * \param qorder_ Quadrature order to use.
         *
         * \note The reference to the parameter objects should be valid for as
         *       long as this localoperators is evaluated.
         */
        R2(Params &params_, std::size_t qorder_ = 2)
          : params(params_), qorder(qorder_)
        {}

        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename M>
        void jacobian_volume(const EG& eg, const LFSU& lfsu, const X& x,
                             const LFSV& lfsv, M& mat) const
        {
          // domain and range field type
          typedef typename LFSU::Traits::FiniteElementType FEU;
          typedef typename LFSV::Traits::FiniteElementType FEV;

          const FEU &feU = lfsu.finiteElement();
          const FEV &feV = lfsv.finiteElement();

          typedef typename FEU::Traits::Basis BasisU;
          typedef typename FEV::Traits::Basis BasisV;

          const BasisU &basisU = feU.basis();
          const BasisV &basisV = feV.basis();

          typedef typename Params::RangeField RF;

          typedef typename BasisU::Traits::DomainField DF;
          static const std::size_t dimDLocal = BasisU::Traits::dimDomainLocal;

          typedef typename BasisU::Traits::Range RangeU;
          typedef typename BasisV::Traits::Range RangeV;

          // select quadrature rule
          typedef QuadratureRules<DF,dimDLocal> QRs;
          typedef QuadratureRule<DF,dimDLocal> QR;
          typedef typename QR::const_iterator QRIterator;
          GeometryType gt = eg.geometry().type();
          const QR& rule = QRs::rule(gt,qorder);

          std::vector<RangeU> phiU(lfsu.size());
          std::vector<RangeV> phiV(lfsv.size());

          // loop over quadrature points
          for(QRIterator it=rule.begin(); it != rule.end(); ++it) {
            // curls of basefunctions
            basisU.evaluateFunction(it->position(), phiU);
            basisV.evaluateFunction(it->position(), phiV);

            RF factor = it->weight()
              * eg.geometry().integrationElement(it->position())
              * params.epsilon(eg.entity(), it->position());

            for(std::size_t j = 0; j < lfsu.size(); ++j)
              for(std::size_t i = 0; i < lfsv.size(); ++i)
                mat(i,j) += factor * (phiU[j] * phiV[i]);
          }
        }

        //! set time on the function object
        void setTime(typename Params::Time time) {
          params.setTime(time);
          IBase::setTime(time);
        }
      };

      //! MultiStepCachePolicy for VectorWave operators
      template<class Params, class Step = int, class Time = double>
      class CachePolicy :
        public MultiStepCachePolicy<Step, Time>
      {
        typedef MultiStepCachePolicy<Step, Time> Base;
        const Params &params;

      protected:
        using Base::currentStep;
        using Base::endTime;
        using Base::dt;

      public:
        //! construct a CachePolicy object
        CachePolicy(const Params &params_) : params(params_) { }

        //! All component operators are affine
        virtual bool isAffine(std::size_t order, Step step) const
        { return true; }
        //! The composed operator is affine
        virtual bool isComposedAffine(Step step) const
        { return true; }
        //! All local operators have pure linear alpha_*()
        virtual bool hasPureLinearAlpha(std::size_t order, Step step) const
        { return true; }
        virtual bool canReuseJacobian(std::size_t order,
                                      Step requested, Step available) const
        {
          Time t1 = endTime - dt*(currentStep-requested);
          Time t2 = endTime - dt*(currentStep-available);
          if(t1 > t2) std::swap(t1, t2);

          switch(order) {
          case 0: return !params.muChanged(t1, t2);
          case 1: return true;
          case 2: return !params.epsilonChanged(t1, t2);
          }

          DUNE_THROW(InvalidStateException,
                     "VectorWave::CachePolicy::canReuseJacobian(): Invalid "
                     "temporal derivative order=" << order);
        }
        virtual bool canReuseZeroResidual(std::size_t order,
                                          Step requested, Step available) const
        { return true; }
        virtual bool canReuseComposedJacobian(Step requested,
                                              Step available) const
        {
          Time t1 = endTime - dt*(currentStep-requested);
          Time t2 = endTime - dt*(currentStep-available);
          if(t1 > t2) std::swap(t1, t2);

          return !params.muChanged(t1, t2) && !params.epsilonChanged(t1, t2);
        }
      };

    } // namespace VectorWave

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_VECTORWAVE_HH
