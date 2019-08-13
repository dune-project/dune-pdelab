// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PKNDDGLOCALCOEFFICIENTS_HH
#define DUNE_PKNDDGLOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>
// #include<dune/pdelab/finiteelement/pkdglagrange/pkdglagrange.hh>

namespace Dune
{


  namespace PkStuff
  {
    // This is the main class
    // usage: QkSize<2,3>::value
    // k is the polynomial degree,
    // n is the space dimension
    // compile time evaluation of number of DoF for simplical DG element

    template<unsigned int k, unsigned int dim>
    struct PkDGSize
    {
      enum{
        value = ((k+dim)*PkDGSize<k,dim-1>::value)/dim
      };
    };

    template<unsigned int k>
    struct PkDGSize<k,0>
    {
      enum{
        value = 1
      };
    };
  }

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for P0 elements

         \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
   */
  template<unsigned int k, unsigned int d>
  class PkDGNDLocalCoefficients
  {
    enum {N = Dune::PkStuff::PkDGSize<k, d>::value};

  public:
    //! \brief Standard constructor
    PkDGNDLocalCoefficients () : li(N)
    {
      for (std::size_t i = 0 ; i < N; i++)
        li[i] = LocalKey(0,0, i);
      }

    //! number of coefficients
    std::size_t size () const
    {
      return N;
    }

    //! get i'th index
    const LocalKey& localKey (std::size_t i) const
    {
      return li[i];
    }

  private:
    std::vector<LocalKey> li;
  };

}

#endif
