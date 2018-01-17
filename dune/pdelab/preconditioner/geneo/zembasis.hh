#ifndef DUNE_GENEO_ZEMBASIS_HH
#define DUNE_GENEO_ZEMBASIS_HH

#include "subdomainbasis.hh"

namespace Dune {
  namespace PDELab {

    /*!
     * \brief A basis implementing zero-energy eigenmodes for linear elasticity.
     */
    template<class GFS, class LFS, class X, int dim, const int d>
    class ZEMBasis : public SubdomainBasis<X>
    {

    public:
      ZEMBasis(const GFS& gfs, LFS& lfs, X& part_unity) {
        using Dune::PDELab::Backend::native;

        this->local_basis.resize(6);

        //First add displacements in each direction
        std::vector<Dune::FieldVector<double,dim>> a(3);
        a[0] = {1., 0. ,0.};
        a[1] = {0., 1., 0.};
        a[2] = {0., 0., 1.};
        for (int k = 0; k < 3; k++){
            this->local_basis[k] = std::make_shared<X>(part_unity);
            for (auto iter = native(*this->local_basis[k]).begin(); iter != native(*this->local_basis[k]).end(); iter++){
                for(int j = 0; j< dim; j++)
                    (*iter)[j] *= a[k][j];
            }
        }

        //Then add rotations
        for(int k = 0; k < 3; k++){
            auto zv = X(gfs,0.0);
            this->local_basis[k+3] = std::make_shared<X>(part_unity);

            for (auto it = gfs.gridView().template begin<0>(); it != gfs.gridView().template end<0>(); ++it) {
                lfs.bind(*it);

                auto geo  = it->geometry();
                const auto gt = geo.type();
                const auto& ref_el = Dune::ReferenceElements<double, d>::general(gt);

                auto& coeffs = lfs.finiteElement().localCoefficients();
                for (std::size_t i = 0; i < coeffs.size(); ++i) {

                    auto local_pos = ref_el.position (coeffs.localKey(i).subEntity(), coeffs.localKey(i).codim());
                    auto global_pos = geo.global(local_pos);

                    auto subindex = gfs.entitySet().indexSet().subIndex(*it, coeffs.localKey(i).subEntity(), coeffs.localKey(i).codim());

                    //ensure each subindex is only written once
                    if(!(native(zv)[subindex][0] > 0.0)){
                        auto val = global_pos;
                        double val_h = val[k];
                        val[k] = val[(k+1)%3];
                        val[(k+1)%2] = val_h;

                        for(int j = 0; j < dim; j++){
                            native(*this->local_basis[k+3])[subindex][j] *= (global_pos[j] - val[j]);
                        }
                        native(zv)[subindex][0] = 1;
                    }
                }
            }
        }
      }

    };

  }
}
#endif //DUNE_GENEO_ZEMBASIS_HH
