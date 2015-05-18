// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMMON_JACOBIANTOCURL_HH
#define DUNE_PDELAB_COMMON_JACOBIANTOCURL_HH

#include <cstddef>

#include <dune/common/fvector.hh>

namespace Dune {
  namespace PDELab {

    //! extract the curl of a function from the jacobian of that function
    /**
     * In 3D the curl \f$A=\nabla\times B\f$ is defined as
     * \f{align*}{
     *   a_x = \partial_yb_z-\partial_zb_y \\
     *   a_y = \partial_zb_x-\partial_xb_z \\
     *   a_z = \partial_xb_y-\partial_yb_x
     * \f}
     * In lower dimensions, some of the coordinates may be missing (because
     * the quantity does not vary in that direction), and some of the
     * quantity's components vanish and are thus missing as well.
     */
    template<typename Jacobian, std::size_t dimR = Jacobian::rows,
             std::size_t dimD = Jacobian::cols>
    class JacobianToCurl;

    //! \brief extract the curl of a 1D-valued function in 2D from the
    //!        jacobian of that function
    /**
     * The two coordinates are the \f$x\f$- and \f$y\f$ coordinates and the
     * one value component is the \f$z\f$-component of the quantity.  It is
     * assumed that the quantity shows no variation in the \f$z\f$-direction
     * (thus \f$\partial_z=0\f$) and that its \f$x\f$- and \f$y\f$-components
     * vanish.  From the general 3D formula for the curl
     * \f{align*}{
     *   A &=\nabla\times B \\
     *     & \Downarrow \\
     *   a_x &= \partial_yb_z-\partial_zb_y \\
     *   a_y &= \partial_zb_x-\partial_xb_z \\
     *   a_z &= \partial_xb_y-\partial_yb_x
     * \f}
     * only the first two survive:
     * \f{align*}{
     *   a_x &=  \partial_yb_z \\
     *   a_y &= -\partial_xb_z
     * \f}
     * Replacing \f$x\f$, \f$y\f$ and \f$z\f$ by the apropriate indices yields
     * \f{align*}{
     *   a_0 &=  \partial_1b_0 \\
     *   a_1 &= -\partial_0b_0
     * \f}
     */
    template<typename Jacobian>
    class JacobianToCurl<Jacobian, 1, 2> {
      static_assert
      ( Jacobian::rows == 1 && Jacobian::cols == 2, "This specialization "
        "works only for dimRange == 1 and dimDomain == 2");

    public:
      typedef typename Jacobian::block_type CurlField;
      static const std::size_t dimCurl = 2;
      typedef FieldVector<CurlField, dimCurl> Curl;

      void operator()(const Jacobian& jacobian, Curl& curl) const {
        curl[0] =  jacobian[0][1];
        curl[1] = -jacobian[0][0];
      }
    };

    //! \brief extract the curl of a 2D-valued function in 2D from the
    //!        jacobian of that function
    /**
     * The two coordinates are the \f$x\f$- and \f$y\f$ coordinates and the
     * two value components the \f$x\f$- and \f$y\f$-components of the
     * quantity.  It is assumed that the quantity shows no variation in the
     * \f$z\f$-direction (thus \f$\partial_z=0\f$) and that its
     * \f$z\f$-component vanishes.  From the general 3D formula for the curl
     * \f{align*}{
     *   A &=\nabla\times B \\
     *     & \Downarrow \\
     *   a_x &= \partial_yb_z-\partial_zb_y \\
     *   a_y &= \partial_zb_x-\partial_xb_z \\
     *   a_z &= \partial_xb_y-\partial_yb_x
     * \f}
     * only the last survives:
     * \f{align*}{
     *   a_z &= \partial_xb_y-\partial_yb_x
     * \f}
     * Replacing \f$x\f$, \f$y\f$ and \f$z\f$ by the apropriate indices yields
     * \f{align*}{
     *   a_0 &= \partial_0b_1-\partial_1b_0
     * \f}
     */
    template<typename Jacobian>
    class JacobianToCurl<Jacobian, 2, 2> {
      static_assert
      ( Jacobian::rows == 2 && Jacobian::cols == 2, "This specialization "
        "works only for dimRange == 2 and dimDomain == 2");

    public:
      typedef typename Jacobian::block_type CurlField;
      static const std::size_t dimCurl = 1;
      typedef FieldVector<CurlField, dimCurl> Curl;

      void operator()(const Jacobian& jacobian, Curl& curl) const {
        curl[0] = jacobian[1][0]-jacobian[0][1];
      }
    };

    //! \brief extract the curl of a 3D-valued function in 3D from the
    //!        jacobian of that function
    /**
     * In the general 3D formula for the curl
     * \f{align*}{
     *   A &=\nabla\times B \\
     *     & \Downarrow \\
     *   a_x &= \partial_yb_z-\partial_zb_y \\
     *   a_y &= \partial_zb_x-\partial_xb_z \\
     *   a_z &= \partial_xb_y-\partial_yb_x
     * \f}
     * we just need to replace \f$x\f$, \f$y\f$ and \f$z\f$ by the apropriate
     * indices
     * \f{align*}{
     *   a_0 &= \partial_1b_2-\partial_2b_1 \\
     *   a_1 &= \partial_2b_0-\partial_0b_2 \\
     *   a_2 &= \partial_0b_1-\partial_1b_0
     * \f}
     */
    template<typename Jacobian>
    class JacobianToCurl<Jacobian, 3, 3> {
      static_assert
      ( Jacobian::rows == 3 && Jacobian::cols == 3, "This specialization "
        "works only for dimRange == 3 and dimDomain == 3");

    public:
      typedef typename Jacobian::block_type CurlField;
      static const std::size_t dimCurl = 3;
      typedef FieldVector<CurlField, dimCurl> Curl;

      void operator()(const Jacobian& jacobian, Curl& curl) const {
        for(std::size_t alpha = 0; alpha < 3; ++alpha) {
          std::size_t beta  = (alpha+1)%3;
          std::size_t gamma = (alpha+2)%3;
          curl[alpha] = jacobian[gamma][beta]-jacobian[beta][gamma];
        }
      }
    };

  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_JACOBIANTOCURL_HH
