// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_INVERSECOEFFICIENTS_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_INVERSECOEFFICIENTS_HH

#include <dune/localfunctions/common/localkey.hh>
#include <dune/common/power.hh>
#include <vector>

namespace Dune{
  namespace PDELab {
    namespace Blockstructured {

      template<int d>
      struct InverseQkLocalCoefficients {
        std::array<std::vector<std::vector<unsigned long>>, d + 1> container;
        int k;
        int blocks;

        template<typename FE>
        InverseQkLocalCoefficients(const FE &fe)
            : k(FE::k), blocks(FE::blocks) {
          setupCodims(fe.type());
          setupSubentities();

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
          }
        }

        void setupSubentities() {
          if (k == 0) {
            container[0][0].resize(Dune::Power<d>::eval(blocks));
          } else {
            container[0][0].resize(Dune::Power<d>::eval(k * blocks - 1));
            setupFaces();
            setupEdges();
            setupVertices();
          }
        }

        void setupFaces() {
          if (d > 2) {
            for (std::size_t i = 0; i < container[1].size(); ++i) {
              container[1][i].resize(Dune::Power<d - 1>::eval(k * blocks - 1));
            }
          }
        }

        void setupEdges() {
          if (d > 1) {
            for (std::size_t i = 0; i < container[d - 1].size(); ++i) {
              container[d - 1][i].resize(k * blocks - 1);
            }
          }
        }

        void setupVertices() {
          for (std::size_t i = 0; i < container[d].size(); ++i) {
            container[d][i].resize(1);
          }
        }
      };
    }
  }
}

#endif //DUNE_PDELAB_GRIDFUNCTIONSPACE_BLOCKSTRUCTURED_INVERSECOEFFICIENTS_HH
