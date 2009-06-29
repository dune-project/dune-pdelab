#ifndef DUNE_PDELAB_Pk3DFEM_HH
#define DUNE_PDELAB_Pk3DFEM_HH

#include<dune/common/exceptions.hh>

#include<dune/finiteelements/pk3d.hh>
#include"finiteelementmap.hh"

namespace Dune
{
    namespace PDELab
    {

        //! wrap up element from local functions
        //! \ingroup FiniteElementMap

        template<typename GV, typename D, typename R, unsigned int k>
        class Pk3DLocalFiniteElementMap :
            public LocalFiniteElementMapInterface<LocalFiniteElementMapTraits< Dune::Pk3DLocalFiniteElement<D,R,k> >,
                                                  Pk3DLocalFiniteElementMap<GV,D,R,k> >,
            public Countable
        {
            typedef Dune::Pk3DLocalFiniteElement<D,R,k> FE;
            typedef typename GV::IndexSet IndexSet;

        public:
            //! \brief export type of the signature
            typedef LocalFiniteElementMapTraits<FE> Traits;

            //! \brief Use when Imp has a standard constructor
            Pk3DLocalFiniteElementMap (const GV& gv_) : is(gv_.indexSet())
            {
                for (int i = 0; i < 256; ++i)
                    perm_index[i] = 0;

                // create all variants by iterating over all permutations
                unsigned int n = 0;
                unsigned int vertexmap[4];
                for(vertexmap[0] = 0; vertexmap[0] < 4; ++vertexmap[0])
                {
                    for(vertexmap[1] = 0; vertexmap[1] < 4; ++vertexmap[1])
                    {
                        if (vertexmap[0] == vertexmap[1])
                            continue;
                        for(vertexmap[2] = 0; vertexmap[2] < 4; ++vertexmap[2])
                        {
                            if (vertexmap[0] == vertexmap[2] ||
                                vertexmap[1] == vertexmap[2])
                                continue;
                            vertexmap[3] = 6 - vertexmap[0] - vertexmap[1] - vertexmap[2];
                            variant[n] = FE(vertexmap);
                            perm_index[compressPerm(vertexmap)] = n++;
                        }
                    }
                }
            }

            //! \brief get local basis functions for entity
            template<class EntityType>
            const typename Traits::LocalFiniteElementType& find (const EntityType& e) const
            {
                // get the vertex indices
                unsigned int vertexmap[4];
                for (unsigned int i = 0; i < 4; ++i)
                    vertexmap[i] = is.subIndex(e,i,3);

                // reduce the indices to the interval 0..3
                for (unsigned int i = 0; i < 4; ++i)
                {
                    int min_index = -1;
                    for (unsigned int j = 0; j < 4; ++j)
                        if ((min_index < 0 || vertexmap[j] < vertexmap[min_index]) && vertexmap[j] >= i)
                            min_index = j;
                    assert(min_index >= 0);
                    vertexmap[min_index] = i;
                }
                return variant[perm_index[compressPerm(vertexmap)]];
            }

        private:
            FE variant[24];
            unsigned int perm_index[256];
            const IndexSet& is;

            unsigned int compressPerm(const unsigned int vertexmap[4]) const
            {
                return vertexmap[0] + (vertexmap[1]<<2) + (vertexmap[2]<<4) + (vertexmap[3]<<6);
            }
        };
    }
}

#endif
