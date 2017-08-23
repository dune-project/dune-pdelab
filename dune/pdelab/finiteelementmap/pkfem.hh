// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_PKFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_PKFEM_HH

#include <algorithm>

#include <dune/common/deprecated.hh>
#include <dune/common/array.hh>
#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/lagrange/pk.hh>

#include "finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    namespace fem {

      template<typename GV, typename D, typename R, unsigned int k, unsigned int d>
      class PkLocalFiniteElementMapBase;


      // ********************************************************************************
      // 1D version
      // ********************************************************************************

      template<typename GV, typename D, typename R, unsigned int k>
      class PkLocalFiniteElementMapBase<GV,D,R,k,1>
        : public SimpleLocalFiniteElementMap<Dune::PkLocalFiniteElement<D,R,1,k>,1>
      {

      public:

        PkLocalFiniteElementMapBase(const GV& gv)
        {}

        static constexpr bool fixedSize()
        {
          return true;
        }

        static constexpr bool hasDOFs(int codim)
        {
          switch (codim)
            {
            case 1: // vertex
              return k > 0;
            case 0: // line
              return k != 1;
            default:
              assert(k >= 0 and k <= 1);
            }
          return false;
        }

        static constexpr std::size_t size(GeometryType gt)
        {
          if (gt == GeometryTypes::vertex)
            return k > 0 ? 1 : 0;
          if (gt == GeometryTypes::line)
            return k > 0 ? k - 1 : 1;
          return 0;
        }

        static constexpr std::size_t maxLocalSize()
        {
          return k + 1;
        }

      };


      // ********************************************************************************
      // 2D version
      // ********************************************************************************

      template<typename GV, typename D, typename R, unsigned int k>
      class PkLocalFiniteElementMapBase<GV,D,R,k,2>
        : public LocalFiniteElementMapInterface<LocalFiniteElementMapTraits<
                                                  Dune::PkLocalFiniteElement<D,R,2,k>
                                                  >,
                                                PkLocalFiniteElementMapBase<GV,D,R,k,2>
                                                >
      {

        typedef Dune::PkLocalFiniteElement<D,R,2,k> FE;

      public:

        //! \brief export type of the signature
        typedef LocalFiniteElementMapTraits<FE> Traits;

        PkLocalFiniteElementMapBase(const GV& gv)
          : _gv(gv)
        {
          // generate permutations
          unsigned int p[3] = {0,1,2};
          for (int i = 0; i < 6; ++i)
            {
              _variant[i] = FE(p);
              std::next_permutation(p,p+3);
            }
        }

        //! \brief get local basis functions for entity
        template<typename Entity>
        const typename Traits::FiniteElementType& find (const Entity& e) const
        {

          if (!e.type().isSimplex())
            DUNE_THROW(InvalidGeometryType,"PkLocalFiniteElementMap only works for simplices!");

          const typename GV::IndexSet& is = _gv.indexSet();
          unsigned int n0 = is.subIndex(e,0,2);
          unsigned int n1 = is.subIndex(e,1,2);
          unsigned int n2 = is.subIndex(e,2,2);
          // compress first number to [0,2]
          unsigned int n0_compressed = (n0 > n1) + (n0 > n2);
          // translate to permutation index
          return _variant[2 * n0_compressed + (n1 > n2)];
        }

        static constexpr bool fixedSize()
        {
          return true;
        }

        static constexpr bool hasDOFs(int codim)
        {
          switch (codim)
            {
            case 2: // vertex
              return k > 0;
            case 1: // line
              return k > 1;
            case 0: // triangle
              return k > 2 || k == 0;
            default:
              assert(false && "Invalid codim specified!");
            }
          return false;
        }

        static constexpr std::size_t size(GeometryType gt)
        {
          if (gt == GeometryTypes::vertex)
            return k > 0 ? 1 : 0;
          if (gt == GeometryTypes::line)
            return k > 1 ? k - 1 : 0;
          if (gt == GeometryTypes::triangle)
            return k > 2 ? (k-2)*(k-1)/2 : (k == 0);
          return 0;
        }

        static constexpr std::size_t maxLocalSize()
        {
          return (k+1)*(k+2)/2;
        }

      private:
        std::array<FE,6> _variant;
        GV _gv;

      };


      // ********************************************************************************
      // 3D version
      // ********************************************************************************

      template<typename GV, typename D, typename R, unsigned int k>
      class PkLocalFiniteElementMapBase<GV,D,R,k,3>
        : public LocalFiniteElementMapInterface<LocalFiniteElementMapTraits<
                                                  Dune::PkLocalFiniteElement<D,R,3,k>
                                                  >,
                                                PkLocalFiniteElementMapBase<GV,D,R,k,3>
                                                >
      {

        typedef Dune::PkLocalFiniteElement<D,R,3,k> FE;

      public:

        //! \brief export type of the signature
        typedef LocalFiniteElementMapTraits<FE> Traits;

        PkLocalFiniteElementMapBase(const GV& gv)
          : _gv(gv)
        {
          std::fill(_perm_index.begin(),_perm_index.end(),0);

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
                      _variant[n] = FE(vertexmap);
                      _perm_index[compressPerm(vertexmap)] = n++;
                    }
                }
            }
        }

        //! \brief get local basis functions for entity
        template<typename Entity>
        const typename Traits::FiniteElementType& find (const Entity& e) const
        {

          if (!e.type().isSimplex())
            DUNE_THROW(InvalidGeometryType,"PkLocalFiniteElementMap only works for simplices!");

          // get the vertex indices
          const typename GV::IndexSet& is = _gv.indexSet();
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
          return _variant[_perm_index[compressPerm(vertexmap)]];
        }

        static constexpr bool fixedSize()
        {
          return true;
        }

        static constexpr bool hasDOFs(int codim)
        {
          switch (codim)
            {
            case 3: // vertex
              return k > 0;
            case 2: // line
              return k > 1;
            case 1: // triangle
              return k > 2;
            case 0: // tetrahedron
              return k == 0 || k > 3;
            default:
              assert(false && "Invalid codim specified!");
            }
          return false;
        }

        static constexpr std::size_t size(GeometryType gt)
        {
          if (gt == GeometryTypes::vertex)
            return k > 0 ? 1 : 0;
          if (gt == GeometryTypes::line)
            return k > 1 ? k - 1 : 0;
          if (gt == GeometryTypes::triangle)
            return k > 2 ? (k-2)*(k-1)/2 : 0;
          if (gt == GeometryTypes::tetrahedron)
            return k == 0 ? 1 : (k-3)*(k-2)*(k-1)/6;
          return 0;
        }

        static constexpr std::size_t maxLocalSize()
        {
          return (k+1)*(k+2)*(k+3)/6;
        }

      private:

        unsigned int compressPerm(const unsigned int vertexmap[4]) const
        {
          return vertexmap[0] + (vertexmap[1]<<2) + (vertexmap[2]<<4) + (vertexmap[3]<<6);
        }

        std::array<FE,24> _variant;
        std::array<unsigned int,256> _perm_index;
        GV _gv;

      };

    } // namespace fem


    template<typename GV, typename D, typename R, unsigned int k>
    class PkLocalFiniteElementMap
      : public fem::PkLocalFiniteElementMapBase<GV,D,R,k,GV::dimension>
    {

    public:

      //! The dimension of the finite elements returned by this map.
      static constexpr int dimension = GV::dimension;

      PkLocalFiniteElementMap(const GV& gv)
        : fem::PkLocalFiniteElementMapBase<GV,D,R,k,GV::dimension>(gv)
      {}

    };


  }
}

#endif // DUNE_PDELAB_FINITEELEMENTMAP_PKFEM_HH
