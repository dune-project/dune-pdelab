// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_INDEXWRAPPER_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_INDEXWRAPPER_HH

#include <vector>
#include <array>
#include <dune/geometry/referenceelements.hh>

namespace Dune{
  namespace PDELab{
    namespace Blockstructured{

      template <typename ContainerOrDOFIndex, int d>
      class SubentityWiseIndexWrapper{
      public:

        SubentityWiseIndexWrapper()
            : storage()
        {
          auto refEl = Dune::ReferenceElements<double,d>::general(Dune::GeometryTypes::cube(d));

          for (int c = 0; c < d + 1; ++c) {
            storage[c].resize(refEl.size(c));
          }
        }

        const ContainerOrDOFIndex& index(const int s, const int c) const
        {
          return storage[c][s];
        }

        ContainerOrDOFIndex& index(const int s, const int c)
        {
          return storage[c][s];
        }

        typename ContainerOrDOFIndex::View indexView(const int s, const int c) const
        {
          return storage[c][s].view();
        }

        void clear()
        {
          for(auto& codim: storage)
            for(auto& subentity: codim)
              subentity.clear();
        }

        auto begin()
        {
          return storage.begin();
        }
        auto cbegin() const
        {
          return storage.cbegin();
        }

        auto end()
        {
          return storage.end();
        }
        auto cend() const
        {
          return storage.cend();
        }

      private:
        std::array<std::vector<ContainerOrDOFIndex>, d + 1> storage;
      };
    }
  }
}

#endif //DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_INDEXWRAPPER_HH
