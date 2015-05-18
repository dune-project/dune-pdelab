// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMMON_CROSSPRODUCT_HH
#define DUNE_PDELAB_COMMON_CROSSPRODUCT_HH

#include <dune/common/fvector.hh>

namespace Dune {
  namespace PDELab {

    //////////////////////////////////////////////////////////////////////
    //
    //  Cross product
    //

    //! Calculate cross product
    /**
     * \f[
     *   A=B\times C
     * \f]
     *
     * \tparam dimB_ Dimension of B.
     * \tparam dimC_ Dimension of C.
     */
    template<unsigned dimB_, unsigned dimC_>
    class CrossProduct {
      static_assert(AlwaysFalse<CrossProduct>::value,
                    "CrossProduct cannot be used unspecialized");
    public:
      //! dimension of \f$A\f$
      static const unsigned dimA;
      //! dimension of \f$B\f$
      static const unsigned dimB;
      //! dimension of \f$C\f$
      static const unsigned dimC;

      //! calculate cross product of B and C
      template<typename AType, typename BType, typename CType>
      CrossProduct(AType& A, const BType& B, const CType& C);
    };

    //! Calculate cross product for dimB=3 and dimC=3
    /**
     * \begin{align*}{
     *   A&=B\times C \\
     *    &\Downarrow \\
     *   a_x&=b_yc_z-b_zc_y \\
     *   a_y&=b_zc_x-b_xc_z \\
     *   a_z&=b_xc_y-b_yc_x \\
     *    &\Downarrow \\
     *   a_0&=b_1c_2-b_2c_1 \\
     *   a_1&=b_2c_0-b_0c_2 \\
     *   a_2&=b_0c_1-b_1c_0
     * \}
     */
    template<>
    struct CrossProduct<3, 3> {
      //! dimension of \f$A\f$
      static const unsigned dimA = 3;
      //! dimension of \f$B\f$
      static const unsigned dimB = 3;
      //! dimension of \f$C\f$
      static const unsigned dimC = 3;

      //! calculate cross product of B and C
      template<typename AType, typename BType, typename CType>
      CrossProduct(AType& A, const BType& B, const CType& C) {
        for(unsigned i = 0; i < 3; ++i) {
          unsigned j = (i+1)%3;
          unsigned k = (i+2)%3;
          A[i] = B[j]*C[k] - B[k]*C[j];
        }
      }
    };

    //! Calculate cross product for dimB=2 and dimC=2
    /**
     * \begin{align*}{
     *   A&=B\times C \\
     *    &\Downarrow \\
     *   a_x&=b_yc_z-b_zc_y \\
     *   a_y&=b_zc_x-b_xc_z \\
     *   a_z&=b_xc_y-b_yc_x
     * \}
     * with
     * \begin{align*}{
     *   b_z&=c_z=0
     *    &\Downarrow \\
     *   a_x&=0 \\
     *   a_y&=0 \\
     *   a_z&=b_xc_y-b_yc_x \\
     *    &\Downarrow \\
     *   a_0&=b_0c_1-b_1c_0
     * \}
     */
    template<>
    struct CrossProduct<2, 2> {
      //! dimension of \f$A\f$
      static const unsigned dimA = 1;
      //! dimension of \f$B\f$
      static const unsigned dimB = 2;
      //! dimension of \f$C\f$
      static const unsigned dimC = 2;

      //! calculate cross product of B and C
      template<typename AType, typename BType, typename CType>
      CrossProduct(AType& A, const BType& B, const CType& C) {
        A[0] = B[0]*C[1] - B[1]*C[0];
      }
    };

    //! Calculate cross product for dimB=2 and dimC=1
    /**
     * \begin{align*}{
     *   A&=B\times C \\
     *    &\Downarrow \\
     *   a_x&=b_yc_z-b_zc_y \\
     *   a_y&=b_zc_x-b_xc_z \\
     *   a_z&=b_xc_y-b_yc_x
     * \}
     * with
     * \begin{align*}{
     *   b_z&=c_x=c_y=0
     *    &\Downarrow \\
     *   a_x&= b_yc_z \\
     *   a_y&=-b_xc_z \\
     *   a_z&=0 \\
     *    &\Downarrow \\
     *   a_0&= b_1c_0 \\
     *   a_1&=-b_0c_0
     * \}
     */
    template<>
    struct CrossProduct<2, 1> {
      //! dimension of \f$A\f$
      static const unsigned dimA = 2;
      //! dimension of \f$B\f$
      static const unsigned dimB = 2;
      //! dimension of \f$C\f$
      static const unsigned dimC = 1;

      //! calculate cross product of B and C
      template<typename AType, typename BType, typename CType>
      CrossProduct(AType& A, const BType& B, const CType& C) {
        A[0] =   B[1]*C[0];
        A[1] = - B[0]*C[0];
      }
    };

    //! Calculate cross product for dimB=1 and dimC=2
    /**
     * \begin{align*}{
     *   A&=B\times C \\
     *    &\Downarrow \\
     *   a_x&=b_yc_z-b_zc_y \\
     *   a_y&=b_zc_x-b_xc_z \\
     *   a_z&=b_xc_y-b_yc_x
     * \}
     * with
     * \begin{align*}{
     *   b_x&=b_y=c_z=0 \\
     *    &\Downarrow \\
     *   a_x&=-b_zc_y \\
     *   a_y&= b_zc_x \\
     *   a_z&= 0 \\
     *    &\Downarrow \\
     *   a_0&=-b_0c_1 \\
     *   a_1&= b_0c_0
     * \}
     */
    template<>
    struct CrossProduct<1, 2> {
      //! dimension of \f$A\f$
      static const unsigned dimA = 2;
      //! dimension of \f$B\f$
      static const unsigned dimB = 1;
      //! dimension of \f$C\f$
      static const unsigned dimC = 2;

      //! calculate cross product of B and C
      template<typename AType, typename BType, typename CType>
      CrossProduct(AType& A, const BType& B, const CType& C) {
        A[0] = - B[0]*C[1];
        A[1] =   B[0]*C[0];
      }
    };

    //////////////////////////////////////////////////////////////////////
    //
    //  Free-standing cross product functions for FieldVectors
    //

    //! calculate crossproduct
    /**
     * Compute the result of \f$b\times c\f$ in \c a.  This does the usual 3D
     * crossproduct:
     * \f[
     *   a_\alpha = \sum_{\beta\gamma}\varepsilon_{\alpha\beta\gamma}
     *                     b_\beta c_\gamma
     * \f]
     */
    template<typename T, int dimB, int dimC>
    void crossproduct(const FieldVector<T, dimB>& b,
                      const FieldVector<T, dimC>& c,
                      FieldVector<T, (CrossProduct<dimB, dimC>::dimA)>& a) {
      CrossProduct<dimB,dimC>(a,b,c);
    }

    //! calculate crossproduct
    /**
     * return the result of \f$b\times c\f$.
     */
    template<typename T, int dimB, int dimC>
    FieldVector<T, CrossProduct<dimB, dimC>::dimA>
    crossproduct(const FieldVector<T, dimB>& b,
                 const FieldVector<T, dimC>& c) {
      FieldVector<T, CrossProduct<dimB, dimC>::dimA> a;
      crossproduct(b,c,a);
      return a;
    }

  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_CROSSPRODUCT_HH
