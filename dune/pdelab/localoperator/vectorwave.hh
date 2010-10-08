// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_LOCALOPERATOR_VECTORWAVE_HH
#define DUNE_PDELAB_LOCALOPERATOR_VECTORWAVE_HH

#include <cstddef>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>
#include <dune/common/static_assert.hh>

#include <dune/grid/common/quadraturerules.hh>

#include <dune/pdelab/common/jacobiantocurl.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/gridoperatorspace/localmatrix.hh>

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{
    namespace VectorWave {

      template<class Imp, class GV, class RF = double, class Time_ = double>
      struct Parameters {
        typedef Time_ Time;

        static const std::size_t dimension = GV::dimension;

        typedef typename GV::ctype DomainField;
        typedef FieldVector<DomainField, dimension> Domain;

        typedef RF RangeField;
        typedef FieldVector<RangeField, dimension> Range;

        typedef typename GV::template Codim<0>::Entity Element;

        RangeField epsilon(const Element& e, const Domain& xl) const
        { return asImp().epsilonGlobal(e.geometry().global(xl)); }
        RangeField mu(const Element& e, const Domain& xl) const
        { return asImp().muGlobal(e.geometry().global(xl)); }

        void setTime(const Time &time) { }

      private:
        const Imp &asImp() const { return *static_cast<const Imp*>(this); }
      };

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
        template<typename Domain>
        RF muGlobal(const Domain &) const { return mu_; }
      };

      //! Construct matrix for a (possibly scaled) curl-of-base product
      /**
       * Construct the matrix
       * \f[
       *    S_{ij}=\int_\Omega\alpha
       *        (\nabla\times\psi_i)\cdot(\nabla\times\psi_j)\,dV
       * \f]
       * where "\f$\dcot\f$" denotes the scalar product and \f$\nabla\times\f$
       * the curl operator.
       *
       * \tparam Time  Type used for time values.
       * \tparam Alpha Type of function to evaluate \f$\alpha\f$.  Must support
       *               setTime().
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

        //! Construct a GlobalVolumeSFunctionRotBaseRotBaseLOP
        /**
         * \param alpha_  Reference to function object to evaluate
         * \param qorder_ Quadrature order to use.
         *
         * \note The references the the function objects should be valid for as
         *       long as this localoperators residual() method is used.
         */
        R0(Params &params_, std::size_t qorder_ = 0) :
          params(params_), qorder(qorder_)
        {}

        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename R>
        void jacobian_volume(const EG& eg, const LFSU& lfsu, const X& x,
                             const LFSV& lfsv, LocalMatrix<R>& mat) const
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

        //! set time on the function object
        void setTime(typename Params::Time time) {
          params.setTime(time);
          IBase::setTime(time);
        }
      };

      //! Construct matrix for a (possibly scaled) curl-of-base product
      /**
       * Construct the matrix
       * \f[
       *    S_{ij}=\int_\Omega\alpha
       *        (\nabla\times\psi_i)\cdot(\nabla\times\psi_j)\,dV
       * \f]
       * where "\f$\dcot\f$" denotes the scalar product and \f$\nabla\times\f$
       * the curl operator.
       *
       * \tparam Time  Type used for time values.
       * \tparam Alpha Type of function to evaluate \f$\alpha\f$.  Must support
       *               setTime().
       */
      template<class Params>
      class R1 :
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename Params::Time>
      {
      public:
        //! Construct a GlobalVolumeSFunctionRotBaseRotBaseLOP
        /**
         * \param alpha_  Reference to function object to evaluate
         * \param qorder_ Quadrature order to use.
         *
         * \note The references the the function objects should be valid for as
         *       long as this localoperators residual() method is used.
         */
        R1(Params &params_, std::size_t qorder_ = 0) {}
      };

      //! Construct matrix for a (possibly scaled) curl-of-base product
      /**
       * Construct the matrix
       * \f[
       *    S_{ij}=\int_\Omega\alpha
       *        (\nabla\times\psi_i)\cdot(\nabla\times\psi_j)\,dV
       * \f]
       * where "\f$\dcot\f$" denotes the scalar product and \f$\nabla\times\f$
       * the curl operator.
       *
       * \tparam Time  Type used for time values.
       * \tparam Alpha Type of function to evaluate \f$\alpha\f$.  Must support
       *               setTime().
       */
      template<class Params>
      class R2 :
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

        //! Construct a GlobalVolumeSFunctionRotBaseRotBaseLOP
        /**
         * \param alpha_  Reference to function object to evaluate
         * \param qorder_ Quadrature order to use.
         *
         * \note The references the the function objects should be valid for as
         *       long as this localoperators residual() method is used.
         */
        R2(Params &params_, std::size_t qorder_ = 2)
          : params(params_), qorder(qorder_)
        {}

        template<typename EG, typename LFSU, typename X, typename LFSV,
                 typename R>
        void jacobian_volume(const EG& eg, const LFSU& lfsu, const X& x,
                             const LFSV& lfsv, LocalMatrix<R>& mat) const
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

    } // namespace VectorWave

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_VECTORWAVE_HH
