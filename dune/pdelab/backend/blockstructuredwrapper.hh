//
// Created by marckoch on 11/05/18.
//

#ifndef DUNE_PDELAB_BLOCKSTRUCTUREDWRAPPER_HH
#define DUNE_PDELAB_BLOCKSTRUCTUREDWRAPPER_HH

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/common/blockstructureduncachedvectorview.hh>
#include <dune/pdelab/backend/istl/vectorhelpers.hh>

namespace Dune{
  namespace Blockstructured{


    template<typename ISTL_BACKEND>
    struct VectorBackendWrapper
        : public ISTL_BACKEND
    {

    };

    template<typename ISTL_CONTAINER>
    class VectorWrapper
        : public ISTL_CONTAINER
    {
    public:

      template<typename LFSCache>
      using LocalView = BlockstructuredUncachedVectorView<ISTL_CONTAINER,LFSCache>;

      template<typename LFSCache>
      using ConstLocalView = BlockstructuredConstUncachedVectorView<const ISTL_CONTAINER,LFSCache>;

      using E = typename ISTL_CONTAINER::E;

      using ISTL_CONTAINER::ISTL_CONTAINER;

      ISTL_CONTAINER& operator= (const VectorWrapper& r)
      {
        return ISTL_CONTAINER::operator=(r);
      }

      ISTL_CONTAINER& operator= (const E& e)
      {
        return ISTL_CONTAINER::operator=(e);
      }

      ISTL_CONTAINER& operator*= (const E& e)
      {
        return ISTL_CONTAINER::operator*=(e);
      }

      ISTL_CONTAINER& operator+= (const E& e)
      {
        return ISTL_CONTAINER::operator+=(e);
      }

      ISTL_CONTAINER& operator+= (const VectorWrapper& e)
      {
        return ISTL_CONTAINER::operator+=(e);
      }

      ISTL_CONTAINER& operator-= (const VectorWrapper& e)
      {
        return ISTL_CONTAINER::operator-=(e);
      }

      E operator*(const VectorWrapper& y) const
      {
        return ISTL_CONTAINER::operator*(y);
      }

      E dot(const VectorWrapper& y) const
      {
        return ISTL_CONTAINER::dot(y);
      }

      ISTL_CONTAINER& axpy(const E& a, const VectorWrapper& y)
      {
        return ISTL_CONTAINER::axpy(a, y);
      }
    };
  }

  namespace PDELab{
    namespace Backend {
      namespace impl {

        template<typename B, typename GFS, typename E>
        struct BackendVectorSelectorHelper<
            Dune::Blockstructured::VectorBackendWrapper<B>,
            GFS,
            E>
        {
          using NativeType = typename Dune::PDELab::Backend::impl::BackendVectorSelectorHelper<B, GFS, E>::Type;
          using Type = Dune::Blockstructured::VectorWrapper<NativeType>;
        };
      }
    }
  }
}

#endif //DUNE_PDELAB_BLOCKSTRUCTUREDWRAPPER_HH
