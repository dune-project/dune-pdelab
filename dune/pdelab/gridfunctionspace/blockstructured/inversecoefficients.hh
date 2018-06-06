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
        std::array<std::vector<std::vector<unsigned long>>, d + 1> container;

        template<typename FE>
        explicit InverseQkLocalCoefficients(const FE &fe) {
          setupCodims(fe.type());
          setupSubentities(fe.localCoefficients());

          const auto &coeffs = fe.localCoefficients();

          for (std::size_t i = 0; i < coeffs.size(); ++i) {
            const auto &l = coeffs.localKey(i);
            container[l.codim()][l.subEntity()][l.index()] = i;
          }
        }

      private:
        void setupCodims(const Dune::GeometryType gt) {
          const auto refEl = Dune::ReferenceElements<double, d>::general(gt);

          for (int c = 0; c < d + 1; ++c) {
            container[c].resize(refEl.size(c));
            subentitySizes[c].resize(refEl.size(c));
          }
        }

        template<typename Coeffs>
        void setupSubentities(const Coeffs& coeffs) {
          for (std::size_t i = 0; i < coeffs.size(); ++i) {
            const auto &l = coeffs.localKey(i);
            subentitySizes[l.codim()][l.subEntity()]++;
          }

          for (int c = 0; c < d + 1; ++c) {
            for (int s = 0; s < container[c].size(); ++s) {
              container[c][s].resize(subentitySizes[c][s]);
            }
          }
        }

        std::array<std::vector<std::size_t>, d + 1> subentitySizes;
      };
    }
  }
}

#endif //DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_INVERSECOEFFICIENTS_HH
