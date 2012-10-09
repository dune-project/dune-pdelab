// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMMON_SIMPLEDOFINDEX_HH
#define DUNE_PDELAB_COMMON_SIMPLEDOFINDEX_HH

#include <dune/common/fvector.hh>

namespace Dune {
  namespace PDELab {


    template<typename F>
    struct SimpleDOFIndex
      : public FieldVector<F,1>
    {

      SimpleDOFIndex()
      {}

      SimpleDOFIndex(const F& v)
        : FieldVector<F,1>(v)
      {}

      F& back()
      {
        return (*this)[0];
      }

      const F& back() const
      {
        return (*this)[0];
      }

    };


    template<typename F>
    struct SimpleContainerIndex
      : public FieldVector<F,1>
    {

      SimpleContainerIndex()
      {}

      SimpleContainerIndex(const F& v)
        : FieldVector<F,1>(v)
      {}

      F& back()
      {
        return (*this)[0];
      }

      const F& back() const
      {
        return (*this)[0];
      }

    };

  } // namespace PDELab
} // namespace Dune



#endif // DUNE_PDELAB_COMMON_SIMPLEDOFINDEX_HH
