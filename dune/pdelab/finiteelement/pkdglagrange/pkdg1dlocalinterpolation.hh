// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P1DG_LOCALINTERPOLATION_HH
#define DUNE_P1DG_LOCALINTERPOLATION_HH

#include <vector>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{
  template<int dim, class LB>
  class PkDG1DLocalInterpolation
  {
  public:
    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;

      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);

      out.resize(dim+1);

      // vertex 0
      for (int i=0; i<dim; i++)
        x[i] = 0;
      out[0] = f(x);

      // remaining vertices
      for (int i=0; i<dim; i++) {
        for (int j=0; j<dim; j++)
          x[j] = (i==j);

        out[i+1] = f(x);

      }

    }

  };
}

#endif
