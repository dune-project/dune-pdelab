// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_MATRIXHELPERS_HH
#define DUNE_PDELAB_BACKEND_ISTL_MATRIXHELPERS_HH

#include<utility>
#include<vector>
#include <unordered_map>
#include <unordered_set>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dune/pdelab/ordering/orderingbase.hh>
#include <dune/pdelab/backend/istl/tags.hh>

namespace Dune {
  namespace PDELab {

#ifndef DOXYGEN

    namespace istl {

      template<typename RV, typename CV, typename block_type>
      struct matrix_for_vectors;

      template<typename B1, typename A1, typename B2, typename A2, typename block_type>
      struct matrix_for_vectors<Dune::BlockVector<B1,A1>,Dune::BlockVector<B2,A2>,block_type>
      {
        typedef Dune::BCRSMatrix<block_type> type;
      };

      template<typename B1, int n1, typename B2, int n2, typename block_type>
      struct matrix_for_vectors<Dune::FieldVector<B1,n1>,Dune::FieldVector<B2,n2>,block_type>
      {
        typedef Dune::FieldMatrix<block_type,n1,n2> type;
      };

      template<typename E, typename RV, typename CV, std::size_t blocklevel>
      struct recursive_build_matrix_type
      {
        typedef typename matrix_for_vectors<RV,CV,typename recursive_build_matrix_type<E,typename RV::block_type,typename CV::block_type,blocklevel-1>::type>::type type;
      };

      template<typename E, typename RV, typename CV>
      struct recursive_build_matrix_type<E,RV,CV,1>
      {
        typedef Dune::FieldMatrix<E,RV::dimension,CV::dimension> type;
      };


      template<typename E, typename RV, typename CV>
      struct build_matrix_type
      {

        static_assert(static_cast<int>(RV::blocklevel) == static_cast<int>(CV::blocklevel),"Both vectors must have identical blocking depth");

        typedef typename recursive_build_matrix_type<E,RV,CV,RV::blocklevel>::type type;

      };


      template<typename RI, typename CI, typename Block>
      typename Block::field_type&
      access_matrix_element(tags::field_matrix_1_1, Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        assert(i == -1);
        assert(j == -1);
        return b[0][0];
      }

      template<typename RI, typename CI, typename Block>
      typename Block::field_type&
      access_matrix_element(tags::field_matrix_n_m, Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        assert(i == 0);
        assert(j == 0);
        return b[ri[0]][ci[0]];
      }

      template<typename RI, typename CI, typename Block>
      typename Block::field_type&
      access_matrix_element(tags::field_matrix_1_m, Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        assert(i == -1);
        assert(j == 0);
        return b[0][ci[0]];
      }

      template<typename RI, typename CI, typename Block>
      typename Block::field_type&
      access_matrix_element(tags::field_matrix_n_1, Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        assert(i == 0);
        assert(j == -1);
        return b[ri[0]][0];
      }

      template<typename RI, typename CI, typename Block>
      typename Block::field_type&
      access_matrix_element(tags::bcrs_matrix, Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        return access_matrix_element(container_tag(b[ri[i]][ci[j]]),b[ri[i]][ci[j]],ri,ci,i-1,j-1);
      }


      template<typename RI, typename CI, typename Block>
      const typename Block::field_type&
      access_matrix_element(tags::field_matrix_1_1, const Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        assert(i == -1);
        assert(j == -1);
        return b[0][0];
      }

      template<typename RI, typename CI, typename Block>
      const typename Block::field_type&
      access_matrix_element(tags::field_matrix_n_m, const Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        assert(i == 0);
        assert(j == 0);
        return b[ri[0]][ci[0]];
      }

      template<typename RI, typename CI, typename Block>
      const typename Block::field_type&
      access_matrix_element(tags::field_matrix_1_m, const Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        assert(i == -1);
        assert(j == 0);
        return b[0][ci[0]];
      }

      template<typename RI, typename CI, typename Block>
      const typename Block::field_type&
      access_matrix_element(tags::field_matrix_n_1, const Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        assert(i == 0);
        assert(j == -1);
        return b[ri[0]][0];
      }

      template<typename RI, typename CI, typename Block>
      const typename Block::field_type&
      access_matrix_element(tags::bcrs_matrix, const Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        return access_matrix_element(container_tag(b[ri[i]][ci[j]]),b[ri[i]][ci[j]],ri,ci,i-1,j-1);
      }



      template<typename RI, typename Block>
      void clear_matrix_row(tags::field_matrix_1_any, Block& b, const RI& ri, int i)
      {
        assert(i == -1);
        b[0] = 0;
      }

      template<typename RI, typename Block>
      void clear_matrix_row(tags::field_matrix_n_any, Block& b, const RI& ri, int i)
      {
        assert(i == 0);
        b[ri[0]] = 0;
      }

      template<typename RI, typename Block>
      void clear_matrix_row(tags::bcrs_matrix, Block& b, const RI& ri, int i)
      {
        typedef typename Block::ColIterator col_iterator_type;
        const col_iterator_type end = b[ri[i]].end();
        for(col_iterator_type cit = b[ri[i]].begin(); cit != end; ++cit)
          clear_matrix_row(container_tag(*cit),*cit,ri,i-1);
      }


      template<typename T, typename RI, typename CI, typename Block>
      void write_matrix_element_if_exists(const T& v, tags::field_matrix_1_1, Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        assert(i == -1);
        assert(j == -1);
        b[0][0] = v;
      }

      template<typename T, typename RI, typename CI, typename Block>
      void write_matrix_element_if_exists(const T& v, tags::field_matrix_n_m, Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        assert(i == 0);
        assert(j == 0);
        b[ri[0]][ci[0]] = v;
      }

      template<typename T, typename RI, typename CI, typename Block>
      void write_matrix_element_if_exists(const T& v, tags::field_matrix_1_m, Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        assert(i == -1);
        assert(j == 0);
        b[0][ci[0]] = v;
      }

      template<typename T, typename RI, typename CI, typename Block>
      void write_matrix_element_if_exists(const T& v, tags::field_matrix_n_1, Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        assert(i == 0);
        assert(j == -1);
        b[ri[0]][0] = v;
      }

      template<typename T, typename RI, typename CI, typename Block>
      void write_matrix_element_if_exists(const T& v, tags::bcrs_matrix, Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        if (b.exists(ri[i],ci[j]))
          write_matrix_element_if_exists(v,container_tag(b[ri[i]][ci[j]]),b[ri[i]][ci[j]],ri,ci,i-1,j-1);
      }


      template<typename OrderingV, typename OrderingU, typename Pattern, typename Container>
      typename std::enable_if<
        !std::is_same<typename Pattern::SubPattern,void>::value &&
      requires_pattern<Container>::value
      >::type
      allocate_matrix(const OrderingV& ordering_v,
                      const OrderingU& ordering_u,
                      const Pattern& p,
                      Container& c)
      {
        c.setSize(ordering_v.blockCount(),ordering_u.blockCount(),false);
        c.setBuildMode(Container::random);

        for (std::size_t i = 0; i < c.N(); ++i)
          c.setrowsize(i,p[i].size());
        c.endrowsizes();

        for (std::size_t i = 0; i < c.N(); ++i)
          for (typename Pattern::value_type::const_iterator cit = p[i].begin(); cit != p[i].end(); ++cit)
            c.addindex(i,cit->first);
        c.endindices();

        for (std::size_t i = 0; i < c.N(); ++i)
          for (typename Pattern::value_type::const_iterator cit = p[i].begin(); cit != p[i].end(); ++cit)
            {
              allocate_matrix(ordering_v.childOrdering(i),
                              ordering_u.childOrdering(cit->first),
                              cit->second,
                              c[i][cit->first]);
            }
      }

      template<typename OrderingV, typename OrderingU, typename Pattern, typename Container>
      typename std::enable_if<
        !std::is_same<typename Pattern::SubPattern,void>::value &&
      !requires_pattern<Container>::value
      >::type
      allocate_matrix(const OrderingV& ordering_v,
                      const OrderingU& ordering_u,
                      const Pattern& p,
                      Container& c)
      {
        for (std::size_t i = 0; i < c.N(); ++i)
          for (typename Pattern::value_type::iterator cit = p[i].begin(); cit != p[i].end(); ++cit)
            {
              allocate_matrix(ordering_v.childOrdering(i),
                              ordering_u.childOrdering(cit->first),
                              cit->second,
                              c[i][cit->first]);
            }
      }

      template<typename OrderingV, typename OrderingU, typename Pattern, typename Container>
      typename std::enable_if<
        std::is_same<typename Pattern::SubPattern,void>::value
        >::type
      allocate_matrix(const OrderingV& ordering_v,
                      const OrderingU& ordering_u,
                      const Pattern& p,
                      Container& c)
      {
        c.setSize(ordering_v.blockCount(),ordering_u.blockCount(),false);
        c.setBuildMode(Container::random);

        for (std::size_t i = 0; i < c.N(); ++i)
          c.setrowsize(i,p[i].size());
        c.endrowsizes();

        for (std::size_t i = 0; i < c.N(); ++i)
          for (typename Pattern::value_type::const_iterator cit = p[i].begin(); cit != p[i].end(); ++cit)
            c.addindex(i,*cit);
        c.endindices();
      }

    } // namespace istl

#endif // DOXYGEN

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_MATRIXHELPERS_HH
