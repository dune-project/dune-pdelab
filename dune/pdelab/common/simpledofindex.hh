// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMMON_SIMPLEDOFINDEX_HH
#define DUNE_PDELAB_COMMON_SIMPLEDOFINDEX_HH

#include <dune/common/fvector.hh>
#include <dune/common/hash.hh>

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

    template<typename F>
    inline std::size_t hash_value(const SimpleDOFIndex<F>& di)
    {
      return di.back();
    }

  } // namespace PDELab
} // namespace Dune

DUNE_DEFINE_HASH(DUNE_HASH_TEMPLATE_ARGS(typename F),DUNE_HASH_TYPE(Dune::PDELab::SimpleDOFIndex<F>))


#endif // DUNE_PDELAB_COMMON_SIMPLEDOFINDEX_HH
