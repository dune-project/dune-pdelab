// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_ELECTRODYNAMIC_HH
#define DUNE_PDELAB_LOCALOPERATOR_ELECTRODYNAMIC_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/fvector.hh>
#include<dune/common/typetraits.hh>

#include<dune/geometry/type.hh>
#include<dune/geometry/quadraturerules.hh>

#include"defaultimp.hh"
#include"pattern.hh"
#include"flags.hh"

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    //! Contruct matrix T for the Electrodynamic operator
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
    public:

      // pattern assembly flags
      enum { doPatternVolume = true };

      enum { doAlphaVolume = true };

      //! Construct an Electrodynamic_T localoperator
      /**
       * \param eps_    Reference to function object to evaluate
       * \param qorder_ Quadrature order to use.
       *
       * \note The references the the function objects should be valid for as
       *       long as this localoperators residual() method is used.
       */
      Electrodynamic_T(const Eps &eps_, int qorder_ = 2)
        : eps(eps_)
        , qorder(qorder_)
      {}

      /**
       * \note We support only Galerkin method lfsu==lfsv
       */
      template<typename EG, typename LFS, typename X, typename M>
      void jacobian_volume (const EG& eg, const LFS& lfsu, const X& x,
                            const LFS& lfsv, M& mat) const
      {
        // domain and range field type
        typedef typename LFS::Traits::FiniteElementType::
          Traits::Basis::Traits BasisTraits;

        typedef typename BasisTraits::DomainField DF;
        static const unsigned dimDLocal = BasisTraits::dimDomainLocal;

        typedef typename BasisTraits::RangeField RF;
        typedef typename BasisTraits::Range Range;
        static const unsigned dimR = BasisTraits::dimRange;

        // static checks
        static_assert(dimR == 3 || dimR == 2,
                      "Works only in 2D or 3D");
        static_assert
          ((std::is_same<typename EG::Geometry::ctype, DF>::value),
           "Grids ctype and Finite Elements DomainFieldType must match");

        // select quadrature rule
        typedef Dune::QuadratureRule<DF,dimDLocal> QR;
        typedef Dune::QuadratureRules<DF,dimDLocal> QRs;
        Dune::GeometryType gt = eg.geometry().type();
        const QR& rule = QRs::rule(gt,qorder);

        // loop over quadrature points
        for(typename QR::const_iterator it=rule.begin();
            it!=rule.end(); ++it) {
          // values of basefunctions
          std::vector<Range> phi(lfsu.size());
          lfsu.finiteElement().basis().evaluateFunction(it->position(),phi);

          // calculate T
          typename Eps::Traits::RangeType epsval;
          eps.evaluate(eg.entity(), it->position(), epsval);

          RF factor = it->weight()
            * eg.geometry().integrationElement(it->position()) * epsval;

          for(unsigned i = 0; i < lfsu.size(); ++i)
            for(unsigned j = 0; j < lfsu.size(); ++j)
              mat.accumulate(lfsv,i,lfsu,j,factor * (phi[i] * phi[j]));
        }
      }

    private:
      const Eps &eps;
      const int qorder;
    };

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
      //! size of FieldVector for holding the curl
      template <unsigned d>
      struct CurlTraits {
        static const unsigned dim =
          d == 1 ? 2 :
          d == 2 ? 1 :
          /*else*/ 3;
      };

      template<typename RF>
      static void
      jacobianToCurl(FieldVector<RF, 1> &curl,
                     const FieldMatrix<RF, 2, 2> &jacobian)
      {
        curl[0] = jacobian[1][0] - jacobian[0][1];
      }
      template<typename RF>
      static void
      jacobianToCurl(FieldVector<RF, 3> &curl,
                     const FieldMatrix<RF, 3, 3> &jacobian)
      {
        for(unsigned i = 0; i < 3; ++i)
          curl[i] = jacobian[(i+2)%3][(i+1)%3] - jacobian[(i+1)%3][(i+2)%3];
      }

    public:

      // pattern assembly flags
      enum { doPatternVolume = true };

      enum { doAlphaVolume = true };

      //! Construct an Electrodynamic_S localoperator
      /**
       * \param mu_     Reference to function object to evaluate
       * \param qorder_ Quadrature order to use.
       *
       * \note The references the the function objects should be valid for as
       *       long as this localoperators residual() method is used.
       */
      Electrodynamic_S(const Mu &mu_, int qorder_ = 2)
        : mu(mu_)
        , qorder(qorder_)
      {}

      /**
       * \note We support only Galerkin method lfsu==lfsv
       */
      template<typename EG, typename LFS, typename X, typename M>
      void jacobian_volume (const EG& eg, const LFS& lfsu, const X& x,
                            const LFS& lfsv, M& mat) const
      {
        // domain and range field type
        typedef typename LFS::Traits::FiniteElementType::Traits::Basis::Traits
          BasisTraits;

        typedef typename BasisTraits::DomainField DF;
        static const unsigned dimDLocal = BasisTraits::dimDomainLocal;

        typedef typename BasisTraits::RangeField RF;
        static const unsigned dimR = BasisTraits::dimRange;

        typedef typename BasisTraits::Jacobian Jacobian;
        typedef FieldVector<RF, CurlTraits<dimR>::dim> Curl;

        // static checks
        static_assert(dimR == 3 || dimR == 2,
                           "Works only in 2D or 3D");
        static_assert
          ((std::is_same<typename EG::Geometry::ctype, DF>::value),
           "Grids ctype and Finite Elements DomainFieldType must match");

        // select quadrature rule
        typedef Dune::QuadratureRule<DF,dimDLocal> QR;
        typedef Dune::QuadratureRules<DF,dimDLocal> QRs;
        Dune::GeometryType gt = eg.geometry().type();
        const QR& rule = QRs::rule(gt,qorder);

        // loop over quadrature points
        for(typename QR::const_iterator it=rule.begin();
            it!=rule.end(); ++it) {
          // curl of the basefunctions
          std::vector<Jacobian> J(lfsu.size());
          lfsu.finiteElement().basis().evaluateJacobian(it->position(),J);

          std::vector<Curl> rotphi(lfsu.size());
          for(unsigned i = 0; i < lfsu.size(); ++i)
            jacobianToCurl(rotphi[i], J[i]);

          // calculate S
          typename Mu::Traits::RangeType muval;
          mu.evaluate(eg.entity(), it->position(), muval);

          RF factor = it->weight()
            * eg.geometry().integrationElement(it->position()) / muval;

          for(unsigned i = 0; i < lfsu.size(); ++i)
            for(unsigned j = 0; j < lfsu.size(); ++j)
              mat.accumulate(lfsv,i,lfsu,j,factor * (rotphi[i] * rotphi[j]));

        }
      }

    private:
      const Mu &mu;
      const int qorder;
    };

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_ELECTRODYNAMIC_HH
