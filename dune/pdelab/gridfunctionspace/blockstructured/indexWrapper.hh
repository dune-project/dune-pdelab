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
      class SubentityIndexWrapper{
      public:

        SubentityIndexWrapper()
            : storage(), codimOffset(), numberOfUsedSubentities(0)
        {
          auto refEl = Dune::ReferenceElements<double,d>::general(Dune::GeometryTypes::cube(d));

          std::size_t offset = 0;
          for (int c = 0; c < d + 1; ++c) {
            codimOffset[c] = offset;
            offset += refEl.size(c);
          }
          numberOfUsedSubentities = offset;
        }

        const ContainerOrDOFIndex& index(const int s, const int c) const
        {
          return storage[codimOffset[c] + s];
        }

        ContainerOrDOFIndex& index(const int s, const int c)
        {
          return storage[codimOffset[c] + s];
        }

        typename ContainerOrDOFIndex::View indexView(const int s, const int c) const
        {
          return storage[codimOffset[c] + s].view();
        }

        void clear()
        {
          for(auto& index: *this)
              index.clear();
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
          return begin() + numberOfUsedSubentities;
        }
        auto cend() const
        {
          return cbegin() + numberOfUsedSubentities;
        }

      private:
        constexpr static std::size_t maxNumberOfSubentities = 27; // Number of subentities for a hexahedron

        std::array<ContainerOrDOFIndex, maxNumberOfSubentities> storage;
        std::array<std::size_t, d + 1> codimOffset;
        std::size_t numberOfUsedSubentities;
      };
    }
  }
}

#endif //DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_INDEXWRAPPER_HH
