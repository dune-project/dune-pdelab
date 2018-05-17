//
// Created by marckoch on 15/05/18.
//

#ifndef DUNE_PDELAB_DOFINDEX_HH
#define DUNE_PDELAB_DOFINDEX_HH

#include <vector>
#include <array>
#include <dune/geometry/referenceelements.hh>

namespace Dune{
  namespace Blockstructured{

    template <typename Index>
    class SubentityWiseIndexWrapper{
    public:
      constexpr static std::size_t dim = 2;

      SubentityWiseIndexWrapper()
          : storage()
      {
        // 2D case
        auto refEl = Dune::ReferenceElements<double,dim>::general(Dune::GeometryTypes::cube(dim));

        for (int c = 0; c < dim + 1; ++c) {
          storage[c].resize(refEl.size(c));
        }
      }

      const Index& index(const int s, const int c) const
      {
        return storage[c][s];
      }

      Index& index(const int s, const int c)
      {
        return storage[c][s];
      }

      typename Index::View indexView(const int s, const int c) const
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
      std::array<std::vector<Index>,3> storage;
    };
  }
}

#endif //DUNE_PDELAB_DOFINDEX_HH
