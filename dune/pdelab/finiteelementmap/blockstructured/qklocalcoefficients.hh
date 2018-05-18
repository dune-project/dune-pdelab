// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_BLOCKSTRUCTURED_QK_LOCALCOEFFICIENTS_HH
#define DUNE_BLOCKSTRUCTURED_QK_LOCALCOEFFICIENTS_HH

#include <array>
#include <cassert>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/power.hh>

#include <dune/localfunctions/common/localkey.hh>
#include <dune/pdelab/gridfunctionspace/blockstructured/lfsindexcache.hh>

namespace Dune
{
  namespace Blockstructured {
    template<int d>
    using SubentityWiseLocalIndexContainer = std::array<std::vector<std::vector<std::size_t>>, d + 1>;


    /** \brief Attaches a shape function to an entity
     *
     * \tparam k Polynomial order
     * \tparam d Dimension of the reference cube
     */
    template<std::size_t k, std::size_t d, std::size_t blocks>
    class QkLocalCoefficients {

      static constexpr std::size_t DOFs1d = k * blocks + 1;

      // Return i as a d-digit number in the (k+1)-nary system
      static std::array<unsigned int, d> multiindex(unsigned int i) {
        std::array<unsigned int, d> alpha;
        for (int j = 0; j < d; j++) {
          alpha[j] = i % DOFs1d;
          i = i / DOFs1d;
        }
        return alpha;
      }

      /** \brief Set the 'subentity' field for each dof for a 1d element */
      void setup1d(std::vector<unsigned int> &subEntity) {
        // Special-handling for piecewise constant elements
        if (k == 0) {
          subEntity[0] = 0;
          return;
        }

        unsigned lastIndex = 0;

        /* edge and vertex numbering
           0----0----1
         */

        // edge (0)
        subEntity[lastIndex++] = 0;                 // corner 0
        for (unsigned i = 0; i < DOFs1d - 1; ++i)
          subEntity[lastIndex++] = 0;               // inner dofs of element (0)

        subEntity[lastIndex++] = 1;                 // corner 1

        assert((StaticPower<DOFs1d, d>::power == lastIndex));
      }

      void setup2d(std::vector<unsigned int> &subEntity) {
        // Special-handling for piecewise constant elements
        if (k == 0) {
          subEntity[0] = 0;
          return;
        }

        unsigned lastIndex = 0;

        // LocalKey: entity number, entity codim, dof indices within each entity
        /* edge and vertex numbering
         2----3----3
         |         |
         |         |
         0         1
         |         |
         |         |
         0----2----1
         */

        // lower edge (2)
        subEntity[lastIndex++] = 0;                 // corner 0
        for (unsigned i = 0; i < DOFs1d - 2; ++i)
          subEntity[lastIndex++] = 2;           // inner dofs of lower edge (2)

        subEntity[lastIndex++] = 1;                 // corner 1

        // iterate from bottom to top over inner edge dofs
        for (unsigned e = 0; e < DOFs1d - 2; ++e) {
          subEntity[lastIndex++] = 0;                   // left edge (0)
          for (unsigned i = 0; i < DOFs1d - 2; ++i)
            subEntity[lastIndex++] = 0;                     // face dofs
          subEntity[lastIndex++] = 1;                   // right edge (1)
        }

        // upper edge (3)
        subEntity[lastIndex++] = 2;                 // corner 2
        for (unsigned i = 0; i < DOFs1d - 2; ++i)
          subEntity[lastIndex++] = 3;                   // inner dofs of upper edge (3)

        subEntity[lastIndex++] = 3;                 // corner 3

        assert((StaticPower<DOFs1d, d>::power == lastIndex));
      }


      void setup3d(std::vector<unsigned int> &subEntity) {
        // Special-handling for piecewise constant elements
        if (k == 0) {
          subEntity[0] = 0;
          return;
        }

        unsigned lastIndex = 0;
#ifndef NDEBUG
        const unsigned numIndices = StaticPower<DOFs1d, d>::power;
        const unsigned numFaceIndices = StaticPower<DOFs1d, d - 1>::power;
#endif
        const unsigned numInnerEdgeDofs = DOFs1d - 2;

        // LocalKey: entity number, entity codim, dof indices within each entity
        /* edge and vertex numbering

                6---(11)--7              6---------7
               /|        /|             /|  (5)   /|
             (8)|      (9)|            / | top   / |
             / (2)     / (3)          /  |(3)bac/k |
            4---(10)--5   |          4---------5   |
            |   |     |   |      left|(0)|     |(1)|right
            |   2--(7)|---3          |   2-----|---3
           (0) /     (1) /           |(2)front |  /
            |(4)      |(5)           | /  (4)  | /
            |/        |/             |/ bottom |/
            0---(6)---1              0---------1
         */

        // bottom face (4)
        lastIndex = 0;
        // lower edge (6)
        subEntity[lastIndex++] = 0;              // corner 0
        for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
          subEntity[lastIndex++] = 6;                // inner dofs of lower edge (6)

        subEntity[lastIndex++] = 1;              // corner 1

        // iterate from bottom to top over inner edge dofs
        for (unsigned e = 0; e < numInnerEdgeDofs; ++e) {
          subEntity[lastIndex++] = 4;                // left edge (4)
          for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
            subEntity[lastIndex++] = 4;                       // inner face dofs
          subEntity[lastIndex++] = 5;                 // right edge (5)
        }

        // upper edge (7)
        subEntity[lastIndex++] = 2;              // corner 2
        for (unsigned i = 0; i < DOFs1d - 2; ++i)
          subEntity[lastIndex++] = 7;                // inner dofs of upper edge (7)
        subEntity[lastIndex++] = 3;                // corner 3

        assert(numFaceIndices == lastIndex);       // added 1 face so far
        /////////////////////////////////////////// end bottom face (4)

        ///////////////////// inner faces
        for (unsigned f = 0; f < numInnerEdgeDofs; ++f) {

          // lower edge (connecting  edges 0 and 1)
          subEntity[lastIndex++] = 0;                // dof on edge 0
          for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
            subEntity[lastIndex++] = 2;                            // dof in front face
          subEntity[lastIndex++] = 1;                // dof on edge 1

          // iterate from bottom to top over inner edge dofs
          for (unsigned e = 0; e < numInnerEdgeDofs; ++e) {
            subEntity[lastIndex++] = 0;                  // on left face (0)
            for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
              subEntity[lastIndex++] = 0;                    // volume dofs
            subEntity[lastIndex++] = 1;                  // right face (1)
          }

          // upper edge (connecting  edges 0 and 1)
          subEntity[lastIndex++] = 2;                // dof on edge 2
          for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
            subEntity[lastIndex++] = 3;                  // dof on rear face (3)
          subEntity[lastIndex++] = 3;                // dof on edge 3

          assert(lastIndex == (f + 1 + 1) * numFaceIndices);
        }

        ////////////////////////////////////////// top face (5)
        // lower edge (10)
        subEntity[lastIndex++] = 4;              // corner 4
        for (unsigned i = 0; i < DOFs1d - 2; ++i)
          subEntity[lastIndex++] = 10;                // inner dofs on lower edge (10)
        subEntity[lastIndex++] = 5;              // corner 5

        // iterate from bottom to top over inner edge dofs
        for (unsigned e = 0; e < DOFs1d - 2; ++e) {
          subEntity[lastIndex++] = 8;                // left edge (8)
          for (unsigned i = 0; i < DOFs1d - 2; ++i)
            subEntity[lastIndex++] = 5;                  // face dofs
          subEntity[lastIndex++] = 9;                // right edge (9)
        }

        // upper edge (11)
        subEntity[lastIndex++] = 6;              // corner 6
        for (unsigned i = 0; i < DOFs1d - 2; ++i)
          subEntity[lastIndex++] = 11;                // inner dofs of upper edge (11)
        subEntity[lastIndex++] = 7;              // corner 7

        assert(numIndices == lastIndex);
      }

    public:

      //! \brief Default constructor
      QkLocalCoefficients() : li(StaticPower<DOFs1d, d>::power) {
        // Set up array of codimension-per-dof-number
        std::vector<unsigned int> codim(li.size());

        for (std::size_t i = 0; i < codim.size(); i++) {
          codim[i] = 0;
          if (k == 0)
            continue;
          // Codimension gets increased by 1 for each coordinate direction
          // where dof is on boundary
          std::array<unsigned int, d> mIdx = multiindex(i);
          for (int j = 0; j < d; j++)
            if (mIdx[j] == 0 or mIdx[j] == (DOFs1d - 1))
              codim[i]++;
        }

        // Set up index vector (the index of the dof in the set of dofs of a given subentity)
        // Algorithm: the 'index' has the same ordering as the dof number 'i'.
        // To make it consecutive we interpret 'i' in the (k+1)-adic system, omit all digits
        // that correspond to axes where the dof is on the element boundary, and transform the
        // rest to the (k-1)-adic system.
        std::vector<unsigned int> index(size());

        for (std::size_t i = 0; i < size(); i++) {

          index[i] = 0;

          std::array<unsigned int, d> mIdx = multiindex(i);

          for (int j = d - 1; j >= 0; j--)
            if (mIdx[j] > 0 and mIdx[j] < (DOFs1d - 1))
              index[i] = (DOFs1d - 2) * index[i] + (mIdx[j] - 1);

        }

        // Set up entity and dof numbers for each (supported) dimension separately
        std::vector<unsigned int> subEntity(li.size());

        if (DOFs1d == 2) {  // We can handle the first-order case in any dimension

          for (std::size_t i = 0; i < size(); i++)
            subEntity[i] = i;

        } else if (d == 1) {

          setup1d(subEntity);

        } else if (d == 2) {

          setup2d(subEntity);

        } else if (d == 3) {

          setup3d(subEntity);

        } else
          DUNE_THROW(Dune::NotImplemented, "QkLocalCoefficients for k==" << k << " and d==" << d);

        // set the number of subentities for each codim
        if (d == 1) {

          il[0].resize(1);
          il[1].resize(2);

        } else if (d == 2) {

          il[0].resize(1);
          il[1].resize(4);
          il[2].resize(4);

        } else if (d == 3) {

          il[0].resize(1);
          il[1].resize(6);
          il[2].resize(12);
          il[3].resize(8);

        }

        // set the number of dof for each (codim,subentity)
        if (k == 0) {

          il[0][0].resize(1);

        } else if (DOFs1d == 2) {

          for (int i = 0; i < il[d].size(); ++i) {
            il[d][i].resize(1);
          }

        } else {

          il[0][0].resize(Dune::StaticPower<DOFs1d - 2, d>::power);
          for (int i = 0; i < il[1].size(); ++i) {
            il[1][i].resize(Dune::StaticPower<DOFs1d - 2, d - 1>::power);
          }
          if (d == 2) {
            for (int i = 0; i < il[2].size(); ++i) {
              il[2][i].resize(1);
            }
          } else if (d == 3) {
            for (int i = 0; i < il[2].size(); ++i) {
              il[2][i].resize(DOFs1d - 2);
            }
            for (int i = 0; i < il[3].size(); ++i) {
              il[3][i].resize(1);
            }
          }
        }

        for (size_t i = 0; i < li.size(); i++) {
          const LocalKey l(subEntity[i], codim[i], index[i]);
          li[i] = l;
          il[l.codim()][l.subEntity()][l.index()] = i;
        }
      }

      //! number of coefficients
      std::size_t size() const {
        return StaticPower<DOFs1d, d>::power;
      }

      //! get i'th index
      const LocalKey &localKey(std::size_t i) const {
        return li[i];
      }

      std::size_t size_index(std::size_t subentity, std::size_t codim) const {
        return il[codim][subentity].size();
      }

      std::size_t localDOF(const LocalKey &l) const {
        return il[l.codim()][l.subEntity()][l.index()];
      }

      const SubentityWiseLocalIndexContainer<d>& getLocalIndexContainer() const
      {
        return il;
      }

    private:
      std::vector<LocalKey> li;

      SubentityWiseLocalIndexContainer<d> il;
    };
  }

}

#endif
