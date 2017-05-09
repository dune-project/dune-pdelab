// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDOPERATOR_COMMON_BORDERDOFEXCHANGER_HH
#define DUNE_PDELAB_GRIDOPERATOR_COMMON_BORDERDOFEXCHANGER_HH

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <dune/common/deprecated.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/geometry/typeindex.hh>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/pdelab/common/borderindexidcache.hh>
#include <dune/pdelab/common/globaldofindex.hh>
#include <dune/pdelab/gridfunctionspace/entityindexcache.hh>

namespace Dune {
  namespace PDELab {


    //! \addtogroup GridOperator
    //! \ingroup PDELab
    //! \{

    //! Helper class for adding up matrix entries on border.
    /**
     *  Utility class for accumulating matrix entries for border-border
     *  couplings in parallel computations with nonoverlapping domain
     *  decomposition.
     *
     *  For nonoverlapping grids, there is an additional problem related to
     *  matrix entries that should exist on a given node, but which are only
     *  present on a remote node. For these entries, we not only need to
     *  accumulate the matrix values, we also have to create those matrix
     *  entries in the sparsity pattern, as they cannot be discovered during
     *  the node-local matrix creation:
     *
     *  We look at a 2D example with a nonoverlapping grid,
     *  two processes and no ghosts with Q1 discretization.
     *  Process 0 has the left part of the domain
     *  with three cells and eight vertices (1-8),
     *  Process 1 the right part with three cells
     *  and eight vertices (2,4,7-12).
     *  <pre>
     *  1 _ 2        2 _ 9 _ 10
     *  |   |        |   |   |
     *  3 _ 4 _ 7    4 _ 7 _ 11
     *  |   |   |        |   |
     *  5 _ 6 _ 8        8 _ 12
     *  </pre>
     *  If we look at vertex 7 and the corresponding entries in the matrix for P0,
     *  there will be entries for (7,4) and (7,8), but not for (7,2).
     *  Unfortunately, local pattern creation will not create the link (7,2) on process
     *  0. This class will find this kind of entry and extend the sparsity pattern appropriately.
     *
     * @tparam GridOperator The grid operator to work on.
     */
    template<typename GridOperator>
    class NonOverlappingBorderDOFExchanger
    {
      typedef typename GridOperator::Traits::Jacobian M;
      typedef M Matrix;
      typedef typename GridOperator::Traits GridOperatorTraits;
      typedef typename GridOperatorTraits::JacobianField Scalar;
      typedef typename GridOperatorTraits::TrialGridFunctionSpace GFSU;
      typedef typename GridOperatorTraits::TestGridFunctionSpace GFSV;
      using EntitySet = typename GFSV::Traits::EntitySet;
      static const int dim = EntitySet::dimension;
      using Grid = typename EntitySet::Traits::GridView::Traits::Grid;
      typedef typename Matrix::block_type BlockType;
      typedef typename Grid::Traits::GlobalIdSet IdSet;
      typedef typename IdSet::IdType IdType;

      //! Extended DOF index, which globally unique
      typedef Dune::PDELab::GlobalDOFIndex<
        typename GFSV::Ordering::Traits::DOFIndex::value_type,
        GFSV::Ordering::Traits::DOFIndex::max_depth,
        typename IdSet::IdType
        > GlobalDOFIndex;

    public:
      //! Data structure for storing border-border matrix pattern entries in a communication-optimized form
      typedef std::unordered_map<
        typename GFSV::Ordering::Traits::DOFIndex,
        std::unordered_set<GlobalDOFIndex>
        > BorderPattern;

    private:
      typedef typename GFSV::Ordering::Traits::DOFIndex RowDOFIndex;
      typedef typename GFSU::Ordering::Traits::DOFIndex ColDOFIndex;

      typedef std::pair<
        typename RowDOFIndex::TreeIndex,
        typename BorderPattern::mapped_type::value_type
        > PatternMPIData;

      typedef std::tuple<
        typename RowDOFIndex::TreeIndex,
        typename BorderPattern::mapped_type::value_type,
        typename M::field_type
        > ValueMPIData;

    public:
      /*!
       * \brief Constructor. Sets up the local to global relations.
       *
       * \param[in] grid_operator The grid operator to access grid view
       *                          and communication cache.
       */
      NonOverlappingBorderDOFExchanger(const GridOperator& grid_operator)
        : _communication_cache(std::make_shared<CommunicationCache>(grid_operator))
        , _entity_set(grid_operator.testGridFunctionSpace().entitySet())
      {}

      void update(const GridOperator& grid_operator)
      {
        _communication_cache = std::make_shared<CommunicationCache>(grid_operator);
      }

      class CommunicationCache
        : public BorderIndexIdCache<GFSV>
      {

        friend class NonOverlappingBorderDOFExchanger;
        typedef BorderIndexIdCache<GFSV> BaseT;

      public:

        CommunicationCache(const GridOperator& go)
          : BaseT(go.testGridFunctionSpace())
          , _gfsu(go.trialGridFunctionSpace())
          , _initialized(false)
          , _entity_cache(go.testGridFunctionSpace())
        {}

        typedef IdType EntityID;
        typedef typename GFSU::Ordering::Traits::DOFIndex::TreeIndex ColumnTreeIndex;
        typedef std::size_t size_type;

        bool initialized() const
        {
          return _initialized;
        }

        void finishInitialization()
        {
          _initialized = true;
        }

        void update()
        {
          BaseT::update();
          _border_pattern.clear();
          _initialized = false;
        }


        const BorderPattern& pattern() const
        {
          assert(initialized());
          return _border_pattern;
        }

        template<typename LFSVCache, typename LFSUCache, typename LocalPattern>
        void addEntries(const LFSVCache& lfsv_cache, const LFSUCache& lfsu_cache, const LocalPattern& pattern)
        {
          assert(!initialized());

          for (typename LocalPattern::const_iterator it = pattern.begin(),
                 end_it = pattern.end();
               it != end_it;
               ++it)
            {
              // skip constrained entries for now. TODO: Is this correct??
              if (lfsv_cache.isConstrained(it->i()) || lfsu_cache.isConstrained(it->j()))
                continue;

              const typename LFSVCache::DOFIndex& di = lfsv_cache.dofIndex(it->i());
              const typename LFSUCache::DOFIndex& dj = lfsu_cache.dofIndex(it->j());

              size_type row_gt_index = GFSV::Ordering::Traits::DOFIndexAccessor::geometryType(di);
              size_type row_entity_index = GFSV::Ordering::Traits::DOFIndexAccessor::entityIndex(di);

              size_type col_gt_index = GFSU::Ordering::Traits::DOFIndexAccessor::geometryType(dj);
              size_type col_entity_index = GFSU::Ordering::Traits::DOFIndexAccessor::entityIndex(dj);

              // We are only interested in connections between two border entities.
              if (!this->isBorderEntity(row_gt_index,row_entity_index) ||
                  !this->isBorderEntity(col_gt_index,col_entity_index))
                continue;

              _border_pattern[di].insert(GlobalDOFIndex(this->id(col_gt_index,col_entity_index),dj.treeIndex()));
            }
        }


        template<typename Entity>
        size_type size(const Entity& e) const
        {
          if (!_gfsu.entitySet().contains(e))
            return 0;
          _entity_cache.update(e);
          size_type n = 0;
          for (size_type i = 0; i < _entity_cache.size(); ++i)
            {
              typename BorderPattern::const_iterator it = _border_pattern.find(_entity_cache.dofIndex(i));
              if (!transfer_dof(i,it))
                continue;
              n += it->second.size();
            }

          return n;
        }

        template<typename Buffer, typename Entity>
        void gather_pattern(Buffer& buf, const Entity& e) const
        {
          if (!_gfsu.entitySet().contains(e))
            return;
          _entity_cache.update(e);
          for (size_type i = 0; i < _entity_cache.size(); ++i)
            {
              typename BorderPattern::const_iterator it = _border_pattern.find(_entity_cache.dofIndex(i));
              if (!transfer_dof(i,it))
                continue;
              for (typename BorderPattern::mapped_type::const_iterator col_it = it->second.begin(),
                     col_end = it->second.end();
                   col_it != col_end;
                   ++col_it)
                buf.write(std::make_pair(_entity_cache.dofIndex(i).treeIndex(),*col_it));
            }
        }

        template<typename Buffer, typename Entity>
        void gather_data(Buffer& buf, const Entity& e, const M& matrix) const
        {
          if (!_gfsu.entitySet().contains(e))
            return;
          _entity_cache.update(e);
          for (size_type i = 0; i < _entity_cache.size(); ++i)
            {
              typename BorderPattern::const_iterator it = _border_pattern.find(_entity_cache.dofIndex(i));
              if (!transfer_dof(i,it))
                continue;
              for (typename BorderPattern::mapped_type::const_iterator col_it = it->second.begin(),
                     col_end = it->second.end();
                   col_it != col_end;
                   ++col_it)
                {
                  typename BaseT::EntityIndex col_entity = this->index(col_it->entityID());

                  ColDOFIndex dj;
                  GFSU::Ordering::Traits::DOFIndexAccessor::store(dj,col_entity.geometryTypeIndex(),col_entity.entityIndex(),col_it->treeIndex());
                  buf.write(std::make_tuple(_entity_cache.dofIndex(i).treeIndex(),*col_it,matrix(_entity_cache.containerIndex(i),_gfsu.ordering().mapIndex(dj))));
                }
            }
        }

      private:

        bool transfer_dof(size_type i, typename BorderPattern::const_iterator it) const
        {
          // not a border DOF
          if (it == _border_pattern.end())
            return false;
          else
            return true;

          /* Constraints check moved to addEntry()
          // check for constraint
          typename GridOperator::Traits::TrialGridFunctionSpaceConstraints::const_iterator cit = _constraints->find(_entity_cache.containerIndex(i));

          // transfer if DOF is not constrained
          // TODO: What about non-Dirichlet constraints??
          return cit == _constraints->end();
          */
        }

        const GFSU& _gfsu;
        BorderPattern _border_pattern;
        bool _initialized;
        mutable EntityIndexCache<GFSV,true> _entity_cache;

      };

      //! A DataHandle class to exchange matrix sparsity patterns
      template<typename Pattern>
      class PatternExtender
        : public CommDataHandleIF<PatternExtender<Pattern>, PatternMPIData>
      {

        typedef std::size_t size_type;

      public:
        //! Export type of data for message buffer
        typedef PatternMPIData DataType;

        bool contains (int dim, int codim) const
        {
          return
            codim > 0 &&
            (_gfsu.dataHandleContains(codim) ||
             _gfsv.dataHandleContains(codim));
        }

        bool fixedsize (int dim, int codim) const
        {
          return false;
        }

        /** @brief How many objects of type DataType have to be sent for a given entity
        */
        template<typename Entity>
        size_type size (Entity& e) const
        {
          if (Entity::codimension == 0)
            return 0;

          return _communication_cache.size(e);
        }

        /** @brief Pack data from user to message buffer
        */
        template<typename MessageBuffer, typename Entity>
        void gather (MessageBuffer& buff, const Entity& e) const
        {
          if (Entity::codimension == 0)
            return;

          _communication_cache.gather_pattern(buff,e);
        }

        /** @brief Unpack data from message buffer to user
        */
        template<typename MessageBuffer, typename Entity>
        void scatter (MessageBuffer& buff, const Entity& e, size_t n)
        {
          if (Entity::codimension == 0)
            return;

          for (size_type i = 0; i < n; ++i)
            {
              DataType data;
              buff.read(data);

              std::pair<bool,typename CommunicationCache::EntityIndex> col_index = _communication_cache.findIndex(data.second.entityID());
              if (!col_index.first)
                continue;

              RowDOFIndex di;
              GFSV::Ordering::Traits::DOFIndexAccessor::store(di,
                                                              e.type(),
                                                              _entity_set.indexSet().index(e),
                                                              data.first);

              ColDOFIndex dj;
              GFSU::Ordering::Traits::DOFIndexAccessor::store(dj,
                                                              col_index.second.geometryTypeIndex(),
                                                              col_index.second.entityIndex(),
                                                              data.second.treeIndex());

              _pattern.add_link(_gfsv.ordering().mapIndex(di),_gfsu.ordering().mapIndex(dj));
            }
        }

        PatternExtender(const NonOverlappingBorderDOFExchanger& dof_exchanger,
                        const GFSU& gfsu,
                        const GFSV& gfsv,
                        Pattern& pattern)
          : _communication_cache(dof_exchanger.communicationCache())
          , _entity_set(dof_exchanger.entitySet())
          , _gfsu(gfsu)
          , _gfsv(gfsv)
          , _pattern(pattern)
        {}

      private:

        const CommunicationCache& _communication_cache;
        const EntitySet _entity_set;
        const GFSU& _gfsu;
        const GFSV& _gfsv;
        Pattern& _pattern;

      };

      //! A DataHandle class to exchange matrix entries
      class EntryAccumulator
        : public CommDataHandleIF<EntryAccumulator,ValueMPIData>
      {

        typedef std::size_t size_type;

      public:
        //! Export type of data for message buffer
        typedef ValueMPIData DataType;

        bool contains(int dim, int codim) const
        {
          return
            codim > 0 &&
            (_gfsu.dataHandleContains(codim) ||
             _gfsv.dataHandleContains(codim));
        }

        bool fixedsize(int dim, int codim) const
        {
          return false;
        }

        template<typename Entity>
        size_type size(Entity& e) const
        {
          if (Entity::codimension == 0)
            return 0;

          return _communication_cache.size(e);
        }

        template<typename MessageBuffer, typename Entity>
        void gather(MessageBuffer& buff, const Entity& e) const
        {
          if (Entity::codimension == 0)
            return;

          _communication_cache.gather_data(buff,e,_matrix);
        }

        /** @brief Unpack data from message buffer to user
        */
        template<typename MessageBuffer, typename Entity>
        void scatter(MessageBuffer& buff, const Entity& e, size_type n)
        {
          if (Entity::codimension == 0)
            return;

          for (size_type i = 0; i < n; ++i)
            {
              DataType data;
              buff.read(data);

              std::pair<bool,typename CommunicationCache::EntityIndex> col_index = _communication_cache.findIndex(std::get<1>(data).entityID());
              if (!col_index.first)
                continue;

              RowDOFIndex di;
              GFSV::Ordering::Traits::DOFIndexAccessor::store(di,
                                                              e.type(),
                                                              _entity_set.indexSet().index(e),
                                                              std::get<0>(data));

              ColDOFIndex dj;
              GFSU::Ordering::Traits::DOFIndexAccessor::store(dj,
                                                              col_index.second.geometryTypeIndex(),
                                                              col_index.second.entityIndex(),
                                                              std::get<1>(data).treeIndex());

              _matrix(_gfsv.ordering().mapIndex(di),_gfsu.ordering().mapIndex(dj)) += std::get<2>(data);
            }
        }


        EntryAccumulator(const NonOverlappingBorderDOFExchanger& dof_exchanger,
                         const GFSU& gfsu,
                         const GFSV& gfsv,
                         Matrix& matrix)
          : _communication_cache(dof_exchanger.communicationCache())
          , _entity_set(dof_exchanger.entitySet())
          , _gfsu(gfsu)
          , _gfsv(gfsv)
          , _matrix(matrix)
        {}

      private:

        const CommunicationCache& _communication_cache;
        EntitySet _entity_set;
        const GFSU& _gfsu;
        const GFSV& _gfsv;
        Matrix& _matrix;

      };

      /**
       * \brief Sums up the entries corresponding to border vertices.
       *
       * \param grid_operator Grid operator.
       * \param matrix Matrix to operate on.
       */
      void accumulateBorderEntries(const GridOperator& grid_operator, Matrix& matrix)
      {
        if (_entity_set.gridView().comm().size() > 1)
          {
            EntryAccumulator data_handle(*this,
                                         grid_operator.testGridFunctionSpace(),
                                         grid_operator.trialGridFunctionSpace(),
                                         matrix);
            _entity_set.gridView().communicate(data_handle,
                                               InteriorBorder_InteriorBorder_Interface,
                                               ForwardCommunication);
          }
      }

      CommunicationCache& communicationCache()
      {
        return *_communication_cache;
      }

      const CommunicationCache& communicationCache() const
      {
        return *_communication_cache;
      }

      shared_ptr<CommunicationCache> communicationCacheStorage()
      {
        return _communication_cache;
      }

      const EntitySet& entitySet() const
      {
        return _entity_set;
      }

    private:

      shared_ptr<CommunicationCache> _communication_cache;
      EntitySet _entity_set;

    };


    template<typename GridOperator>
    class NoDataBorderDOFExchanger
    {

    public:

      typedef NoDataBorderDOFExchanger CommunicationCache;

      //! Data structure for storing border-border matrix pattern entries in a communication-optimized form
      typedef Empty BorderPattern;

      NoDataBorderDOFExchanger()
      {}

      NoDataBorderDOFExchanger(const GridOperator& grid_operator)
      {}

      void accumulateBorderEntries(const GridOperator& grid_operator, typename GridOperator::Traits::Jacobian& matrix)
      {}

      CommunicationCache& communicationCache()
      {
        return *this;
      }

      const CommunicationCache& communicationCache() const
      {
        return *this;
      }

      void update(const GridOperator& grid_operator)
      {}

    };


    template<typename GridOperator>
    class OverlappingBorderDOFExchanger :
      public NoDataBorderDOFExchanger<GridOperator>
    {

    public:

      OverlappingBorderDOFExchanger()
      {}

      OverlappingBorderDOFExchanger(const GridOperator& grid_operator)
      {}

    };


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDOPERATOR_COMMON_BORDERDOFEXCHANGER_HH
