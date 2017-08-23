// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENTMAP_QKDG_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_QKDG_HH

#include <dune/pdelab/finiteelement/qkdglagrange.hh>
#include <dune/pdelab/finiteelement/qkdglegendre.hh>
#include <dune/pdelab/finiteelement/qkdglobatto.hh>

#include <dune/pdelab/finiteelementmap/opbfem.hh>
#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>

namespace Dune {
  namespace PDELab {

    //! Switch between different basis for QkDGLocalFiniteElementMap
    enum class QkDGBasisPolynomial { lagrange, legendre, lobatto, l2orthonormal};

#ifndef DOXYGEN
    // Class declaration. Use template specialization below.
    template<class D, class R, int k, int d, QkDGBasisPolynomial p = QkDGBasisPolynomial::lagrange>
    class QkDGLocalFiniteElementMap;
#endif

    /** \brief Qk discontinuous Galerkin FiniteElementMap based on
     * Lagrange polynomials
     *
     * \tparam D Type used for coordinates
     * \tparam R Type used for shape function values
     * \tparam k Order of polynomial basis
     * \tparam d Grid dimension
     *
     * \ingroup FiniteElementMap
     */
    template<class D, class R, int k, int d>
    class QkDGLocalFiniteElementMap<D,R,k,d,QkDGBasisPolynomial::lagrange>
      : public Dune::PDELab::SimpleLocalFiniteElementMap< Dune::QkDGLagrangeLocalFiniteElement<D,R,k,d>,d>
    {
    public:

      static constexpr bool fixedSize()
      {
        return true;
      }

      static constexpr bool hasDOFs(int codim)
      {
        return codim == 0;
      }

      static constexpr std::size_t size(GeometryType gt)
      {
        if (gt == GeometryTypes::cube(d))
          return Dune::QkStuff::QkSize<k,d>::value;
        else
          return 0;
      }

      static constexpr std::size_t maxLocalSize()
      {
        return Dune::QkStuff::QkSize<k,d>::value;
      }

      //! return type of polynomial basis
      static constexpr QkDGBasisPolynomial polynomial()
      {
        return QkDGBasisPolynomial::lagrange;
      }

      //! return order of polynomial basis
      static constexpr std::size_t order()
      {
        return k;
      }
    };

    /** \brief Qk discontinuous Galerkin FiniteElementMap based on
     * Legendre polynomials
     *
     * \tparam D Type used for coordinates
     * \tparam R Type used for shape function values
     * \tparam k Order of polynomial basis
     * \tparam d Grid dimension
     *
     * \ingroup FiniteElementMap
     */
    template<class D, class R, int k, int d>
    class QkDGLocalFiniteElementMap<D,R,k,d,QkDGBasisPolynomial::legendre>
      : public Dune::PDELab::SimpleLocalFiniteElementMap< Dune::QkDGLegendreLocalFiniteElement<D,R,k,d>,d>
    {
    public:

      static constexpr bool fixedSize()
      {
        return true;
      }

      static constexpr bool hasDOFs(int codim)
      {
        return codim == 0;
      }

      static constexpr std::size_t size(GeometryType gt)
      {
        if (gt == GeometryTypes::cube(d))
          return Dune::QkStuff::QkSize<k,d>::value;
        else
          return 0;
      }

      static constexpr std::size_t maxLocalSize()
      {
        return Dune::LegendreStuff::LegendreSize<k,d>::value;
      }

      //! return type of polynomial basis
      static constexpr QkDGBasisPolynomial polynomial()
      {
        return QkDGBasisPolynomial::legendre;
      }

      //! return order of polynomial basis
      static constexpr std::size_t order()
      {
        return k;
      }
    };

     /** \brief Qk discontinuous Galerkin FiniteElementMap based on
     * Legendre polynomials at Gauss-Lobatto points
     *
     * \tparam D Type used for coordinates
     * \tparam R Type used for shape function values
     * \tparam k Order of polynomial basis
     * \tparam d Grid dimension
     *
     * \ingroup FiniteElementMap
     */
    template<class D, class R, int k, int d>
    class QkDGLocalFiniteElementMap<D,R,k,d,QkDGBasisPolynomial::lobatto>
      : public Dune::PDELab::SimpleLocalFiniteElementMap< Dune::QkDGGLLocalFiniteElement<D,R,k,d>,d>
    {
    public:

      static constexpr bool fixedSize()
      {
        return true;
      }

      static constexpr bool hasDOFs(int codim)
      {
        return codim == 0;
      }

      static constexpr std::size_t size(GeometryType gt)
      {
        if (gt == GeometryTypes::cube(d))
          return Dune::QkStuff::QkSize<k,d>::value;
        else
          return 0;
      }

      static constexpr std::size_t maxLocalSize()
      {
        return Dune::QkStuff::QkSize<k,d>::value;
      }

      //! return type of polynomial basis
      static constexpr QkDGBasisPolynomial polynomial()
      {
        return QkDGBasisPolynomial::lobatto;
      }

      //! return order of polynomial basis
      static constexpr std::size_t order()
      {
        return k;
      }
    };


    /** \brief Qk discontinuous Galerkin FiniteElementMap based on
     * an L2 orthonormal polynomials
     *
     * If you have gmp the computation field for the l2 orthonormal
     * basis polynomials is Dune::GMPField<512>.
     *
     * \tparam D Type used for coordinates
     * \tparam R Type used for shape function values
     * \tparam k Order of polynomial basis
     * \tparam d Grid dimension
     *
     * \note Use OPBLocalFiniteElementMap directly if you need more
     * customization.
     *
     * \ingroup FiniteElementMap
     */
   template<class D, class R, int k, int d>
    class QkDGLocalFiniteElementMap<D,R,k,d,QkDGBasisPolynomial::l2orthonormal>
      : public OPBLocalFiniteElementMap<D,R,k,d,Dune::GeometryType::cube,
#if HAVE_GMP
                                        Dune::GMPField<512>,
#else
                                        R,
#endif
                                        Dune::PB::BasisType::Qk>
    {
    public:

      //! return type of polynomial basis
      static constexpr QkDGBasisPolynomial polynomial()
      {
        return QkDGBasisPolynomial::l2orthonormal;
      }

      //! return order of polynomial basis
      static constexpr std::size_t order()
      {
        return k;
      }

    };

  }
}

#endif // DUNE_PDELAB_FINITEELEMENTMAP_QKDG_HH
