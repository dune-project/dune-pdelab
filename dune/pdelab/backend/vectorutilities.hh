// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMMON_VECTORUTILITIES_HH
#define DUNE_PDELAB_COMMON_VECTORUTILITIES_HH

#warning "dune/pdelab/backend/vectorutilities.hh" is deprecated.
#warning It has not been ported to the new DOF infrastructure and will be removed shortly.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>

namespace Dune {
  namespace PDELab {

    //! make sure all vector entries are below a certain limit
    template<class Vector>
    DUNE_DEPRECATED bool checkVectorLimit(const Vector &v, typename Vector::ElementType limit)
    {
      typedef typename Vector::Backend BE;
      for(std::size_t i = 0; i < v.flatsize(); ++i)
        // This condition is carefully crafted such that 'return false' will
        // always happen for infinite or NaN entries.
        if(!(std::abs(BE::access(v, i)) < limit))
          return false;
      return true;
    }

    //! make sure all vector entries are below a certain limit
    template<class Communication, class Vector>
    DUNE_DEPRECATED bool checkVectorLimit(const Communication &comm, const Vector &v,
                                          typename Vector::ElementType limit)
    {
      bool good = checkVectorLimit(v, limit);
      return comm.min(good);
    }

    //! Calculate infinity norm of all vector entries.
    template<class Vector>
    DUNE_DEPRECATED typename Vector::ElementType vectorInfNorm(const Vector &v) {
      typedef typename Vector::Backend BE;
      typename Vector::ElementType norm = 0;
      for(std::size_t i = 0; i < v.flatsize(); ++i)
        norm = std::max(norm, std::abs(BE::access(v, i)));
      return norm;
    }

    //! Calculate infinity norm of all vector entries.
    template<class Communication, class Vector>
    DUNE_DEPRECATED typename Vector::ElementType vectorInfNorm(const Communication &comm,
                                                               const Vector &v)
    {
      typename Vector::ElementType norm = vectorInfNorm(v);
      return comm.max(norm);
    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_COMMON_VECTORUTILITIES_HH
