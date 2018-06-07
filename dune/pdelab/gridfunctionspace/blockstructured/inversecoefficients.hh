// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_INVERSECOEFFICIENTS_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_INVERSECOEFFICIENTS_HH

#include <vector>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/common/power.hh>
#include <dune/geometry/referenceelements.hh>

namespace Dune{
  namespace PDELab {
    namespace Blockstructured {

      template<int d>
      struct InverseQkLocalCoefficients {
        template<typename FE>
        explicit InverseQkLocalCoefficients(const FE &fe) {
          setupCodims(fe.type());
          setupSubentities(fe.localCoefficients());

          const auto &coeffs = fe.localCoefficients();

          for (std::size_t i = 0; i < coeffs.size(); ++i) {
            localDOF(coeffs.localKey(i)) = i;
          }
        }

        std::size_t localDOF(const Dune::LocalKey& l) const
        {
          return container[subentityOffset[l.codim()][l.subEntity()] + l.index()];
        }

        std::size_t& localDOF(const Dune::LocalKey& l)
        {
          return container[subentityOffset[l.codim()][l.subEntity()] + l.index()];
        }


        std::size_t size(std::size_t codim) const
        {
          return subentitySizes[codim].size();
        }

        std::size_t size(std::size_t subentity, std::size_t codim) const
        {
          return subentitySizes[codim][subentity];
        }

      private:
        void setupCodims(const Dune::GeometryType gt) {
          const auto refEl = Dune::ReferenceElements<double, d>::general(gt);

          for (int c = 0; c < d + 1; ++c) {
            subentitySizes[c].resize(refEl.size(c));
            subentityOffset[c].resize(refEl.size(c));
          }
        }

        template<typename Coeffs>
        void setupSubentities(const Coeffs& coeffs) {
          for (std::size_t i = 0; i < coeffs.size(); ++i) {
            const auto &l = coeffs.localKey(i);
            subentitySizes[l.codim()][l.subEntity()]++;
          }

          std::size_t offset = 0;
          for (int c = 0; c < d + 1; ++c) {
            for (int s = 0; s < subentitySizes[c].size(); ++s) {
              subentityOffset[c][s] = offset;
              offset += subentitySizes[c][s];
            }
          }

          container.resize(offset);
        }

        std::vector<std::size_t> container;
        std::array<std::vector<std::size_t>, d + 1> subentitySizes;
        std::array<std::vector<std::size_t>, d + 1> subentityOffset;
      };
    }
  }
}

#endif //DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_INVERSECOEFFICIENTS_HH
