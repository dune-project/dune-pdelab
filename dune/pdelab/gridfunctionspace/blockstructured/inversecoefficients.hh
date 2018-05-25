// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_INVERSECOEFFICIENTS_HH
#define DUNE_PDELAB_INVERSECOEFFICIENTS_HH

#include <dune/localfunctions/common/localkey.hh>
#include <dune/common/power.hh>
#include <vector>

namespace Dune{
  namespace Blockstructured{

    template<int d>
    struct InverseQkLocalCoefficients
    {
      std::array<std::vector<std::vector<unsigned long>>, d + 1> container;
      const int k;
      const int blocks;

      template <typename FE>
      InverseQkLocalCoefficients(const FE& fe)
          : k(FE::k), blocks(FE::blocks)
      {
        setupCodims();
        setupSubentities();

        const auto& coeffs = fe.localCoefficients();

        for (int i = 0; i < coeffs.size(); ++i) {
          const auto& l = coeffs.localKey(i);
          container[l.codim()][l.subEntity()][l.index()] = i;
        }
      }

      unsigned long localIndex(const Dune::LocalKey& l) const
      {
        return container[l.codim()][l.subEntity()][l.index()];
      }


    private:
      void setupCodims() {
        switch (d){
          case 1:
            setupCodims1d();
            break;
          case 2:
            setupCodims2d();
            break;
          case 3:
            setupCodims3d();
            break;
        }
      }
      void setupSubentities(){
        if(k == 0) {
          container[0][0].resize(Dune::Power<d>::eval(blocks));
        }
        else {
          container[0][0].resize(Dune::Power<d>::eval(k * blocks - 1));
          setupFaces();
          setupEdges();
          setupVertices();
        }
      }

      void setupCodims1d(){
        container[0].resize(1);
        container[1].resize(2);
      }
      void setupCodims2d(){
        container[0].resize(1);
        container[1].resize(4);
        container[2].resize(4);
      }
      void setupCodims3d(){
        container[0].resize(1);
        container[1].resize(6);
        container[2].resize(12);
        container[3].resize(8);
      }

      void setupFaces() {
        if (d > 2) {
          for (int i = 0; i < container[1].size(); ++i) {
            container[d - 2][i].resize(Dune::Power<d - 1>::eval(k * blocks - 1));
          }
        }
      }
      void setupEdges() {
        if (d > 1){
          for (int i = 0; i < container[1].size(); ++i) {
            container[d - 1][i].resize(k * blocks - 1);
          }
        }
      }
      void setupVertices() {
        for (int i = 0; i < container[d].size(); ++i) {
          container[d][i].resize(1);
        }
      }
    };
  }
}

#endif //DUNE_PDELAB_INVERSECOEFFICIENTS_HH
