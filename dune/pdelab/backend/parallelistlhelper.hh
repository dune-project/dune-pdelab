// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PARALLELISTLHELPER_HH
#define DUNE_PARALLELISTLHELPER_HH

#include <dune/common/deprecated.hh>
#include <dune/common/mpihelper.hh>

#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/io.hh>
#include <dune/istl/superlu.hh>

#include "../constraints/constraints.hh"
#include "../gridfunctionspace/genericdatahandle.hh"
#include "../newton/newton.hh"
#include "istlvectorbackend.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup Backend
    //! \ingroup PDELab
    //! \{

    //========================================================
    // A parallel helper class providing a nonoverlapping
    // decomposition of all degrees of freedom
    //========================================================

    // operator that resets result to zero at constrained DOFS
    template<typename GFS>
    class ParallelISTLHelper
    {
      /**
       * @brief Writes 1<<24 to each data item (of the container) that is gathered or scattered
       * and is neither interior nor border.
       *
       * Can be used to mark ghost cells.
       */
      class GhostGatherScatter
      {
      public:
        template<class MessageBuffer, class EntityType, class DataType>
        void gather (MessageBuffer& buff, const EntityType& e, DataType& data)
        {
          if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
            data = (1<<24);
          buff.write(data);
        }

        template<class MessageBuffer, class EntityType, class DataType>
        void scatter (MessageBuffer& buff, const EntityType& e, DataType& data)
        {
          DataType x;
          buff.read(x);
          if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
            data = (1<<24);
        }
      };

      /**
       * @brief GatherScatter handle that sets 1<<24 for data items neither associated to
       * the interior or border and take the minimum when scattering.
       *
       * Used to compute an owner rank for each unknown.
       */
      class InteriorBorderGatherScatter
      {
      public:
        template<class MessageBuffer, class EntityType, class DataType>
        void gather (MessageBuffer& buff, const EntityType& e, DataType& data)
        {
          if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
            data = (1<<24);
          buff.write(data);
        }

        template<class MessageBuffer, class EntityType, class DataType>
        void scatter (MessageBuffer& buff, const EntityType& e, DataType& data)
        {
          DataType x;
          buff.read(x);
          if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
            data = x;
          else
            data = std::min(data,x);
        }
      };

      /**
       * @brief GatherScatter handle for finding out about neighbouring processor ranks.
       *
       */
      template<typename T>
      struct NeighbourGatherScatter
      {
        NeighbourGatherScatter(int rank_, std::set<int>& neighbours_)
          : myrank(rank_), neighbours(neighbours_)
        {}

        template<class MessageBuffer, class DataType>
        void gather (MessageBuffer& buff, DataType& data)
        {
          buff.write(myrank);
        }

        template<class MessageBuffer, class DataType>
        void scatter (MessageBuffer& buff, DataType& data)
        {
          DataType x;
          buff.read(x);
          neighbours.insert((int)x);
        }

        T myrank;
        std::set<int>& neighbours;
      };


      /**
       * @brief GatherScatter handle for finding out about neighbouring processor ranks.
       *
       */
      struct SharedGatherScatter
      {
        template<class MessageBuffer, class DataType>
        void gather (MessageBuffer& buff, DataType& data)
        {
          data=true;
          buff.write(data);
        }

        template<class MessageBuffer, class DataType>
        void scatter (MessageBuffer& buff, DataType& data)
        {
          bool x;
          buff.read(x);
          data = data || x;
        }
      };

      /**
       * @brief GatherScatter handle for finding out about neighbouring processor ranks.
       *
       */
      template<typename B, typename V1>
      struct GlobalIndexGatherScatter
      {
        GlobalIndexGatherScatter(const V1& mask_)
          : mask(mask_)
        {}

        template<class MessageBuffer, class DataType>
        void gather (MessageBuffer& buff, typename B::size_type i, DataType& data)
        {
          //if(B::access(mask, i)>0)
          if(data < std::numeric_limits<DataType>::max())
            // We now the global index and therefore write it
            buff.write(data);
          else
            buff.write(std::numeric_limits<DataType>::max());
        }

        template<class MessageBuffer, class DataType>
        void scatter (MessageBuffer& buff, typename B::size_type i, DataType& data)
        {
          DataType x;
          buff.read(x);
          data = std::min(data, x);
        }
        V1 mask;
      };

      
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,double>::Type V;

    public:

      ParallelISTLHelper (const GFS& gfs_, int verbose_=1)
        : gfs(gfs_), v(gfs,(double)gfs.gridview().comm().rank()), g(gfs,0.0), verbose(verbose_)
      {
        // find out about ghosts
        Dune::PDELab::GenericDataHandle2<GFS,V,GhostGatherScatter> gdh(gfs,v,GhostGatherScatter());
        if (gfs.gridview().comm().size()>1)
          gfs.gridview().communicate(gdh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);

        g = v;

        // partition interior/border
        Dune::PDELab::GenericDataHandle2<GFS,V,InteriorBorderGatherScatter> dh(gfs,v,InteriorBorderGatherScatter());
        if (gfs.gridview().comm().size()>1)
          gfs.gridview().communicate(dh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);

        // convert vector into mask vector
        for (typename V::size_type i=0; i<v.base().N(); ++i)
          for (typename V::size_type j=0; j<v.base()[i].N(); ++j)
            if (v.base()[i][j]==gfs.gridview().comm().rank())
              v.base()[i][j] = 1.0;
            else
              v.base()[i][j] = 0.0;

      }

      // keep only DOFs assigned to this processor
      template<typename W>
      void mask (W& w) const
      {
        for (typename V::size_type i=0; i<v.N(); ++i)
          for (typename V::size_type j=0; j<v[i].N(); ++j)
            w.base()[i][j] *= v.base()[i][j];
      }

      // access to mask vector
      double mask (typename V::size_type i, typename V::size_type j) const
      {
        return v.base()[i][j];
      }

      // access to ghost vector
      double ghost (typename V::size_type i, typename V::size_type j) const
      {
        return g.base()[i][j];
      }

#if HAVE_MPI

      /**
       * @brief Creates a matrix suitable for parallel AMG and the parallel information
       *
       * It is silently assumed that the unknows are associated with vertices.
       *
       * @tparam MatrixType The type of the ISTL matrix used.
       * @tparam Comm The type of the OwnerOverlapCopyCommunication
       * @param m The local matrix.
       * @param c The parallel information object providing index set, interfaces and
       * communicators.
       */
      template<typename MatrixType, typename Comm>
      void createIndexSetAndProjectForAMG(MatrixType& m, Comm& c);
#endif
    private:
      const GFS& gfs;
      V v; // vector to identify unique decomposition
      V g; //vector to identify ghost dofs
      int verbose; //verbosity
    };


    namespace
    {
      template<typename GFS, bool Comp>
      struct BlockwiseIndicesHelper
      {
        enum{ value = false };
      };

      template<typename M, typename B, int k>
      struct BlockSizeIsEqual
      {
        enum{ value = false };
      };

      template<typename M, int k>
      struct BlockSizeIsEqual<M,ISTLVectorBackend<k>,k>
      {
        enum{ value = false };
        static_assert(AlwaysFalse<M>::value, "Unsupported GridFunctionSpace mapper");
      };

      template<int k>
      struct BlockSizeIsEqual<GridFunctionSpaceComponentBlockwiseMapper<1>,ISTLVectorBackend<k>,k>
      {
        enum{ value = true };
      };

      template<int k>
      struct BlockSizeIsEqual<ComponentBlockwiseOrderingTag<1>,ISTLVectorBackend<k>,k>
      {
        enum{ value = true };
      };

      template<typename GFS>
      struct BlockwiseIndicesHelper<GFS,true>
      {
        enum{ value =  BlockSizeIsEqual<typename GFS::Traits::MapperType,
              typename GFS::Traits::BackendType, GFS::CHILDREN>::value};
      };

      template<typename GFS>
      struct BlockwiseIndices
      {
        enum{
          value = BlockwiseIndicesHelper<GFS,GFS::Traits::isComposite>::value
        };
      };

      template<typename GFS, bool b, int k>
      struct BlockProcessorHelper
      {};

      template<typename GFS>
      struct BlockProcessorHelper<GFS,false,1>
      {

        template<typename T>
        struct AMGVectorTypeSelector
        {
          typedef typename T::BaseT  Type;
        };

        template<typename T>
        static typename AMGVectorTypeSelector<T>::Type& getVector(T& t)
        {
          return t.base();
        }

        template<typename G>
        static void postProcessCount(G& g)
        {}

        template<typename G>
        static void increment(G& g, std::size_t i)
        {
          ++g;
        }
        template<typename M, typename TI>
        static void addIndex(const typename TI::GlobalIndex& gi, std::size_t i,
                             typename TI::LocalIndex::Attribute attr, M& m, TI& idxset)
        {
          // Add index
          idxset.add(gi, typename TI::LocalIndex(i, attr));
        }
      };

      template<typename GFS>
      struct BlockProcessorHelper<GFS,true,1>
        : public BlockProcessorHelper<GFS,false,1>
      {};

      template<typename GFS, int k>
      struct BlockProcessorHelper<GFS, true, k>
      {
        template<typename T>
        struct AMGVectorTypeSelector
        {
          typedef typename T::BaseT Type;
        };

        template<typename T>
        static typename AMGVectorTypeSelector<T>::Type& getVector(T& t)
        {
          return t.base();
        }
        template<typename G>
        static void postProcessCount(G& g)
        {
          g=g;
        }

        template<typename G>
        static void increment(G& g, std::size_t i)
        {
          if((i+1)%GFS::Traits::noChilds==0)
            ++g;
        }
        template<typename M, typename TI>
        static void addIndex(const typename TI::GlobalIndex& gi, std::size_t i,
                             typename TI::LocalIndex::Attribute attr, M& m, TI& idxset)
        {
          if(i%GFS::Traits::noChilds==0)
            BlockProcessorHelper<GFS, false, 1>::addIndex(gi, i/GFS::Traits::noChilds,
                                                          attr, m, idxset);
        }

      };

      template<typename GFS>
      struct BlockProcessor
        : public BlockProcessorHelper<GFS, BlockwiseIndices<GFS>::value,
                                      GFS::Traits::BackendType::BlockSize>
      {
      };

    } // end anonymous namspace


#if HAVE_MPI
    template<typename GFS>
    template<typename M, typename C>
    void ParallelISTLHelper<GFS>::createIndexSetAndProjectForAMG(M& m, C& c)
    {
      typedef typename GFS::Traits::GridViewType GV;
      const GV& gv = gfs.gridview();
      //      static const std::size_t dim = GV::Grid::dimension;
      //      if(gv.comm().size()>1 && gv.grid().overlapSize(dim)<1)
      //        DUNE_THROW(Dune::InvalidStateException, "ParallelISTLHelper::createIndexSetAndProjectForAMG: "
      //                   <<"Only grids with at least one layer of overlap cells are supported");

      // First find out which dofs we share with other processors
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,bool>::Type BoolVector;
      BoolVector sharedDOF(gfs, false);
      Dune::PDELab::GenericDataHandle<GFS,BoolVector,SharedGatherScatter> gdh(gfs,sharedDOF,SharedGatherScatter());

      if (gfs.gridview().comm().size()>1)
        gfs.gridview().communicate(gdh,Dune::All_All_Interface,Dune::ForwardCommunication);

      // Count shared dofs that we own
      typedef typename C::ParallelIndexSet::GlobalIndex GlobalIndex;
      GlobalIndex count=0;
      std::size_t noScalars=0;

      for (typename V::size_type i=0; i<v.N(); ++i)
        for (typename V::size_type j=0; j<v[i].N(); ++j, ++noScalars)
          if(v[i][j]==1.0 && sharedDOF[i][j])
            ++count;

      if (verbose > 1)
        std::cout<<gv.comm().rank()<<": shared count is "<< count.touint()<<std::endl;

      // Maybe divide by block size?
      BlockProcessor<GFS>::postProcessCount(count);
      if (verbose > 1)
      std::cout<<gv.comm().rank()<<": shared block count is "<< count.touint()<<std::endl;

      std::vector<GlobalIndex> counts(gfs.gridview().comm().size());
      MPI_Allgather(&count, 1, MPITraits<GlobalIndex>::getType(), &(counts[0]),
                    1, MPITraits<GlobalIndex>::getType(),
                    gfs.gridview().comm());

      // Compute start index start_p = \sum_{i=0}^{i<p} counts_i
      GlobalIndex start=0;
      for(int i=0; i<gfs.gridview().comm().rank(); ++i)
        start=start+counts[i];
      //std::cout<<gv.comm().rank()<<": start index = "<<start.touint()<<std::endl;

      
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,GlobalIndex>::Type GIVector;
      GIVector scalarIndices(gfs, std::numeric_limits<GlobalIndex>::max());


      for (typename V::size_type i=0, ii=0; i<v.N(); ++i)
        for (typename V::size_type j=0; j<v[i].N(); ++j)
          if(v[i][j]==1.0 && sharedDOF[i][j]){
            scalarIndices[i][j]=start;
            BlockProcessor<GFS>::increment(start, ii++);
          }

      // publish global indices for the shared DOFS to other processors.
      typedef GlobalIndexGatherScatter<typename GFS::Traits::BackendType,V> GIGS;
      Dune::PDELab::GenericDataHandle3<GFS,GIVector,GIGS> gdhgi(gfs, scalarIndices, GIGS(v));
       if (gfs.gridview().comm().size()>1)
                gfs.gridview().communicate(gdhgi,Dune::All_All_Interface,Dune::ForwardCommunication);


      // Setup the index set
      c.indexSet().beginResize();
      for (typename V::size_type i=0, ii=0; i<v.N(); ++i)
        for (typename V::size_type j=0; j<v[i].N(); ++j, ++ii){
          Dune::OwnerOverlapCopyAttributeSet::AttributeSet attr;
          if(scalarIndices[i][j]!=std::numeric_limits<GlobalIndex>::max()){
            // global index exist in index set
            if(v[i][j]>0){
              // This dof is managed by us.
              attr = Dune::OwnerOverlapCopyAttributeSet::owner;
            }
            else if ( g[i][j]==(1<<24) && ( c.getSolverCategory() == 
                                            static_cast<int>(SolverCategory::nonoverlapping)) ){
              //use attribute overlap for ghosts in novlp grids
              attr = Dune::OwnerOverlapCopyAttributeSet::overlap;
            }
            else {
              attr = Dune::OwnerOverlapCopyAttributeSet::copy;                
            }
            BlockProcessor<GFS>::
              addIndex(scalarIndices[i][j], ii, attr, m, c.indexSet());
          }
        }
      c.indexSet().endResize();
      //std::cout<<gv.comm().rank()<<": index set size = "<<c.indexSet().size()<<std::endl;
      //std::cout<<gv.comm().rank()<<": "<<c.indexSet()<<std::endl;

      // Compute neighbours using communication
      typedef NeighbourGatherScatter<typename V::ElementType> NeighbourGS;
      std::set<int> neighbours;
      Dune::PDELab::GenericDataHandle<GFS,V,NeighbourGS> gdhn(gfs, v, NeighbourGS(gfs.gridview().comm().rank(),
                                                                                  neighbours));
      if (gfs.gridview().comm().size()>1)
        gfs.gridview().communicate(gdhn,Dune::All_All_Interface,Dune::ForwardCommunication);
      c.remoteIndices().setNeighbours(neighbours);
      //std::cout<<gv.comm().rank()<<": no neighbours="<<neighbours.size()<<std::endl;

      c.remoteIndices().template rebuild<false>();
      //std::cout<<c.remoteIndices()<<std::endl;
       
    }
#endif

    template<int s, bool isFakeMPIHelper>
    struct CommSelector
    {
      typedef Dune::Amg::SequentialInformation type;
    };

    // Need MPI for OwnerOverlapCopyCommunication
#if HAVE_MPI
    template<int s>
    struct CommSelector<s,false>
    {
      typedef OwnerOverlapCopyCommunication<bigunsignedint<s>,int> type;
    };
#endif

    //! \} group Backend

  } // namespace PDELab
} // namespace Dune

#endif
