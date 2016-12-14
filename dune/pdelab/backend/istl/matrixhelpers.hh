// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_MATRIXHELPERS_HH
#define DUNE_PDELAB_BACKEND_ISTL_MATRIXHELPERS_HH

// this is here for backwards compatibility and deprecation warnings, remove after 2.5.0
#include "ensureistlinclude.hh"

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

    namespace ISTL {

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

      template<typename RowOrdering, typename ColOrdering, typename SubPattern_ = void>
      class Pattern
        : public std::vector<std::unordered_map<std::size_t,SubPattern_> >
      {

      public:

        typedef SubPattern_ SubPattern;

        template<typename RI, typename CI>
        void add_link(const RI& ri, const CI& ci)
        {
          recursive_add_entry(ri.view(),ci.view());
        }

        template<typename RI, typename CI>
        void recursive_add_entry(const RI& ri, const CI& ci)
        {
          this->resize(_row_ordering.blockCount());
          std::pair<typename std::unordered_map<std::size_t,SubPattern>::iterator,bool> r = (*this)[ri.back()].insert(make_pair(ci.back(),SubPattern(_row_ordering.childOrdering(ri.back()),_col_ordering.childOrdering(ci.back()))));
          r.first->second.recursive_add_entry(ri.back_popped(),ci.back_popped());
        }

        Pattern(const RowOrdering& row_ordering, const ColOrdering& col_ordering)
          : _row_ordering(row_ordering)
          , _col_ordering(col_ordering)
        {}

      private:

        const RowOrdering& _row_ordering;
        const ColOrdering& _col_ordering;

      };

      template<typename RowOrdering, typename ColOrdering>
      class Pattern<RowOrdering,ColOrdering,void>
        : public std::vector<std::unordered_set<std::size_t> >
      {

      public:

        typedef void SubPattern;

        template<typename RI, typename CI>
        void add_link(const RI& ri, const CI& ci)
        {
          recursive_add_entry(ri,ci);
        }

        template<typename RI, typename CI>
        void recursive_add_entry(const RI& ri, const CI& ci)
        {
          this->resize(_row_ordering.blockCount());
          (*this)[ri.back()].insert(ci.back());
        }

        Pattern(const RowOrdering& row_ordering, const ColOrdering& col_ordering)
          : _row_ordering(row_ordering)
          , _col_ordering(col_ordering)
        {}

      private:

        const RowOrdering& _row_ordering;
        const ColOrdering& _col_ordering;

      };

      template<typename M, int blocklevel = M::blocklevel>
      struct requires_pattern
      {
        static const bool value = requires_pattern<typename M::block_type,blocklevel-1>::value;
      };

      template<typename M>
      struct requires_pattern<M,0>
      {
        static const bool value = false;
      };

      template<typename B, typename A, int blocklevel>
      struct requires_pattern<Dune::BCRSMatrix<B,A>,blocklevel>
      {
        static const bool value = true;
      };

      template<typename M, typename RowOrdering, typename ColOrdering, bool pattern>
      struct _build_pattern_type
      {
        typedef void type;
      };

      template<typename M, typename RowOrdering, typename ColOrdering>
      struct _build_pattern_type<M,RowOrdering,ColOrdering,true>
      {
        typedef Pattern<RowOrdering,ColOrdering,typename _build_pattern_type<typename M::block_type,RowOrdering,ColOrdering,requires_pattern<typename M::block_type>::value>::type> type;
      };

      template<typename M, typename GFSV, typename GFSU, typename Tag>
      struct build_pattern_type
      {

        typedef OrderingBase<
          typename GFSV::Ordering::Traits::DOFIndex,
          typename GFSV::Ordering::Traits::ContainerIndex
          > RowOrdering;

        typedef OrderingBase<
          typename GFSU::Ordering::Traits::DOFIndex,
          typename GFSU::Ordering::Traits::ContainerIndex
          > ColOrdering;

        typedef typename _build_pattern_type<M,RowOrdering,ColOrdering,requires_pattern<M>::value>::type type;
      };

      template<typename M, typename GFSV, typename GFSU>
      struct build_pattern_type<M,GFSV,GFSU,FlatContainerAllocationTag>
      {
        typedef Pattern<typename GFSV::Ordering, typename GFSU::Ordering> type;
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


      template<typename RI, typename Block>
      void clear_matrix_row_block(tags::field_matrix_1_1, Block& b, const RI& ri, int i)
      {
        assert(i == -1);
        b = 0;
      }

      template<typename RI, typename Block>
      void clear_matrix_row_block(tags::field_matrix_1_any, Block& b, const RI& ri, int i)
      {
        DUNE_THROW(Dune::Exception,"Should never get here!");
      }

      template<typename RI, typename Block>
      void clear_matrix_row_block(tags::field_matrix_n_any, Block& b, const RI& ri, int i)
      {
        assert(i == 0);
        b = 0;
      }

      template<typename RI, typename Block>
      void clear_matrix_row_block(tags::bcrs_matrix, Block& b, const RI& ri, int i)
      {
        typedef typename Block::ColIterator col_iterator_type;
        const col_iterator_type end = b[ri[i]].end();
        for(col_iterator_type cit = b[ri[i]].begin(); cit != end; ++cit)
          clear_matrix_row_block(container_tag(*cit),*cit,ri,i-1);
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




      template<typename T, typename RI, typename CI, typename Block>
      void write_matrix_element_if_exists_to_block(const T& v, tags::field_matrix_1_1, Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        assert(i == -1);
        assert(j == -1);
        b[0][0] = v;
      }

      template<typename T, typename RI, typename CI, typename Block>
      void write_matrix_element_if_exists_to_block(const T& v, tags::field_matrix_n_m, Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        assert(i == 0);
        assert(j == 0);
        for (std::size_t i = 0; i < b.rows; ++i)
          b[i][i] = v;
      }

      template<typename T, typename RI, typename CI, typename Block>
      void write_matrix_element_if_exists_to_block(const T& v, tags::field_matrix_1_m, Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        DUNE_THROW(Dune::Exception,"Should never get here!");
      }

      template<typename T, typename RI, typename CI, typename Block>
      void write_matrix_element_if_exists_to_block(const T& v, tags::field_matrix_n_1, Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        DUNE_THROW(Dune::Exception,"Should never get here!");
      }

      template<typename T, typename RI, typename CI, typename Block>
      void write_matrix_element_if_exists_to_block(const T& v, tags::bcrs_matrix, Block& b, const RI& ri, const CI& ci, int i, int j)
      {
        if (b.exists(ri[i],ci[j]))
          write_matrix_element_if_exists_to_block(v,container_tag(b[ri[i]][ci[j]]),b[ri[i]][ci[j]],ri,ci,i-1,j-1);
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

    } // namespace ISTL

#endif // DOXYGEN

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_MATRIXHELPERS_HH
