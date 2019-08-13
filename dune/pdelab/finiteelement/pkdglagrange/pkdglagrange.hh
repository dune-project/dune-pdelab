// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENT_PKDGLAGRANGE_HH
#define DUNE_PDELAB_FINITEELEMENT_PKDGLAGRANGE_HH

#include <numeric>
#include "pkdg1d.hh"
#include "pkdg2d.hh"
#include "pkdg3d.hh"




namespace Dune
{


  /** \brief General Lagrange finite element with arbitrary dimension and polynomial order
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam d dimension of the reference element
   * \tparam k polynomial order
   */
  template<class D, class R, int d, int k>
  class PkDGLocalFiniteElement
  {
  public:
    PkDGLocalFiniteElement()
    {}

    /** Constructor for variants with permuted vertices.

        \param vertexmap The permutation of the vertices.  This
        can for instance be generated from the global indices of
        the vertices by reducing those to the integers 0...k+1
     */
    PkDGLocalFiniteElement(const unsigned int vertexmap[k+1])
    {}
  };
  /** \brief General Lagrange finite element -- specialization for a 2d reference element
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam k polynomial order
   */
  template<class D, class R, int k>
  class PkDGLocalFiniteElement<D, R, 1, k>
    : public PkDG1DLocalFiniteElement<D, R, k>
  {
  public:
    PkDGLocalFiniteElement()
    {}

    PkDGLocalFiniteElement(const unsigned int vertexmap[2]) :
      Pk1DLocalFiniteElement<D, R, k>(vertexmap)
    {}
  };

  /** \brief General Lagrange finite element -- specialization for a 2d reference element
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam k polynomial order
   */
  template<class D, class R, int k>
  class PkDGLocalFiniteElement<D, R, 2, k>
    : public PkDG2DLocalFiniteElement<D, R, k>
  {
  public:
    PkDGLocalFiniteElement()
    {}

    PkDGLocalFiniteElement(const unsigned int vertexmap[3]) :
      Pk2DLocalFiniteElement<D, R, k>(vertexmap)
    {}
  };

  /** \brief General Lagrange finite element -- specialization for a 3d reference element
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam k polynomial order
   */
  template<class D, class R, int k>
  class PkDGLocalFiniteElement<D, R, 3, k>
    : public PkDG3DLocalFiniteElement<D, R, k>
  {
  public:
    PkDGLocalFiniteElement()
    {}

    PkDGLocalFiniteElement(const unsigned int vertexmap[4]) :
      PkDG3DLocalFiniteElement<D, R, k>(vertexmap)
    {}
  };




};


#endif // DUNE_PDELAB_FINITEELEMENT_QKDGLAGRANGE_HH
