// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_EXTEND_OVERLAP_TOOLS_HH
#define DUNE_PDELAB_EXTEND_OVERLAP_TOOLS_HH

// dune-common includes
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/enumset.hh>
#include<dune/common/parallel/indexset.hh>
#include<dune/common/parallel/plocalindex.hh>
#include<dune/common/parallel/interface.hh>
#include<dune/common/parallel/remoteindices.hh>
#include<dune/common/parallel/communicator.hh>
#include<dune/common/parallel/variablesizecommunicator.hh>

#include "communicator_with_rank.hh"
#include "variablesizecommunicator_with_rank.hh"

namespace Dune {



  // nested name space for classes local to this file
  namespace OverlapTools {

    template<typename V>
    struct AddGatherScatter
    {
      static typename V::value_type gather(const V& a, int i)
      {
        return a[i]; // I am sending my value
      }
      static void scatter (V& a, typename V::value_type v, int i)
      {
        a[i]+=v; // add what I receive to my value
      }
    };

    // data handle to construct owner partitioning on the degrees of freedom
    template<typename ScalarVector>
    class MakeOwnerDataHandle
    {
      const ScalarVector& distance;
      int rank;
      int p;
      ScalarVector& owner;

    public:
      // export type of items to be sent around
      using DataType = int;

      // constructor
      MakeOwnerDataHandle (const ScalarVector& distance_,  // gives distance from the boundary for each dof; owners must have distance >0
                           int rank_, // my rank
                           int p_, // total number of procs
                           ScalarVector& owner_ // value 1 indicates this index is owned by that rank, otherwise it is 0
                           )
        : distance(distance_), rank(rank_), p(p_), owner(owner_)
      {
        owner = 1.0;
      }

      // enable variable size
      bool fixedsize()
      {
        return true;
      }

      // return size for given (local) index
      std::size_t size (int i)
      {
        return 1; // we need a variable communicator to pass and out in our data
      }

      // gather to buffer
      template<class B>
      void gather (B& buffer, int i)
      {
        if (distance[i]>0.0)
          buffer.write(rank); // send my number if I could possibly own it
        else
          {
            buffer.write(p); // send a number that is larger than any rank
            owner[i] = 0.0; // we will never own it
          }
      }

      template<class B>
      void scatter(B& buffer, int i, int count)
      {
        DataType otherrank;
        buffer.read(otherrank);
        if (distance[i]>0.0 && otherrank<rank)
          owner[i] = 0.0; // the other owns it since its rank is smaller
      }
    };


    // data handle to find out degrees of freedom on the boundary with other subdomains
    template<typename ParallelIndexSet, typename Matrix, typename Vector>
    class CountNonlocalEdgesDataHandle
    {
      const Matrix& M;
      const ParallelIndexSet& pis;
      using GlobalId = typename ParallelIndexSet::GlobalIndex;
      std::vector<GlobalId> globalid;
      Vector& z;

    public:
      // export type of items to be sent around
      using DataType = GlobalId;

      // constructor
      CountNonlocalEdgesDataHandle (const Matrix& M_,  // the fully constructed matrix in the overlapping subdomains
                                    const ParallelIndexSet& pis_, // our parallel index set
                                    Vector& z_ // every entry will contain the number of graph edges that are NOT present in this rank
                                    )
        : M(M_), pis(pis_), globalid(pis_.size()), z(z_)
      {
        // construct map from local index to global id on all indices
        auto endit = pis.end();
        for (auto it=pis.begin(); it!=endit; ++it)
          globalid[it->local().local()] = it->global();

        // initialize z with zero
        z = 0.0;
      }

      // enable variable size
      bool fixedsize()
      {
        return false;
      }

      // return size for given (local) index
      std::size_t size (int i)
      {
        std::size_t count = 0; // we dont need the trick because there is at least one nonzero per row
        auto MEndIt = M[i].end();
        for (auto It = M[i].begin(); It!=MEndIt; ++It)
          count++;
        //std::cout << "sending " << count << " entries for index " << i << std::endl;
        return count;
      }

      // gather to buffer
      template<class B>
      void gather(B& buffer, int i)
      {
        auto MEndIt = M[i].end();
        for (auto It = M[i].begin(); It!=MEndIt; ++It)
          buffer.write(globalid[It.index()]);
      }

      template<class B>
      void scatter(B& buffer, int i, int count)
      {
        DataType gid;
        for (int k=0; k<count; k++)
          {
            buffer.read(gid);
            if (!pis.exists(gid))
              z[i] += 1.0;
          }
      }
    };


    // data handle to find out degrees of freedom on the boundary with other subdomains
    template<typename Matrix, typename ParallelIndexSet, typename Vector>
    class AddNonlocalEntriesToDiagonal
    {
      using GlobalId = typename ParallelIndexSet::GlobalIndex;
      using BlockType = typename Matrix::block_type;

      Matrix& M;
      const ParallelIndexSet& pis;
      const Vector& owner;
      std::vector<GlobalId> globalid;

    public:
      // export type of items to be sent around
      using DataType = std::pair<BlockType,GlobalId>;

      // constructor
      AddNonlocalEntriesToDiagonal (Matrix& M_,  // the fully constructed matrix in the overlapping subdomains
                                    const ParallelIndexSet& pis_, // our parallel index set
                                    const Vector& owner_) // who owns a dof has the complete row, only owner sends
        : M(M_), pis(pis_), owner(owner_), globalid(pis_.size())
      {
        // construct map from local index to global id on all indices
        auto endit = pis.end();
        for (auto it=pis.begin(); it!=endit; ++it)
          globalid[it->local().local()] = it->global();
      }

      // enable variable size
      bool fixedsize()
      {
        return false;
      }

      // return size for given (local) index
      std::size_t size (int i)
      {
        std::size_t count = 1; // the trick, send at least one entry
        if (owner[i]==0.0) return count;
        auto MEndIt = M[i].end();
        for (auto It = M[i].begin(); It!=MEndIt; ++It)
          count++;
        return count;
      }

      // gather to buffer
      template<class B>
      void gather(B& buffer, int i)
      {
        buffer.write(std::make_pair(BlockType(),GlobalId())); // write a default value
        if (owner[i]==0.0) return;
        auto MEndIt = M[i].end();
        for (auto It = M[i].begin(); It!=MEndIt; ++It)
          buffer.write(std::make_pair(*It,globalid[It.index()]));
      }

      template<class B>
      void scatter(B& buffer, int i, int count)
      {
        DataType entry;
        buffer.read(entry); // skip default entry
        for (int k=1; k<count; k++)
          {
            buffer.read(entry);
            if (!pis.exists(entry.second)) // I don't have this matrix entry;
              M[i][i] += entry.first;
          }
      }
    };


    // data handle to insert matrix entries in the overlap region
    template<typename Matrix, typename ParallelIndexSet>
    class MatrixConstructDataHandle
    {
      using LocalIndex = typename Matrix::size_type;
      using GlobalId = typename ParallelIndexSet::GlobalIndex;
      Matrix& M;       // the new matrix on the overlapping index set
      const Matrix& A; // the original input matrix
      std::shared_ptr<ParallelIndexSet> pis;
      const std::vector<LocalIndex>& new2old_localindex;
      const std::vector<LocalIndex>& old2new_localindex;
      const std::vector<GlobalId>& globalid;

      using BlockType = typename Matrix::block_type;

    public:
      // export type of items to be sent around
      using DataType = std::pair<BlockType,GlobalId>;

      // constructor
      MatrixConstructDataHandle (Matrix& M_, // the new matrix
                                 const Matrix& A_, // the input matrix in additive form
                                 std::shared_ptr<ParallelIndexSet> pis_,
                                 const std::vector<LocalIndex>& new2old_localindex_,
                                 const std::vector<LocalIndex>& old2new_localindex_,
                                 const std::vector<GlobalId>& globalid_)
        : M(M_), A(A_), pis(pis_), new2old_localindex(new2old_localindex_), old2new_localindex(old2new_localindex_), globalid(globalid_)
      {}

      // enable variable size
      bool fixedsize()
      {
        return false;
      }

      // return size for given (local) index
      std::size_t size (int i, int proc) // i is a new local index
      {
        std::size_t count = 1; // always send one
        if ((unsigned)i>=new2old_localindex.size()) return count; // not our business since this index is not in the original range
        auto I = new2old_localindex[i]; // I is the local index in the input index set
        auto MEndIt = A[I].end();
        for (auto It = A[I].begin(); It!=MEndIt; ++It) // iterate over row and count entries
          if (old2new_localindex[It.index()]<new2old_localindex.size()) // column index is also non ghost
            count++;
        //std::cout << "sending " << count << " entries for index " << i << std::endl;
        return count;
      }

      // gather to buffer
      template<class B>
      void gather(B& buffer, int i, int proc)
      {
        //std::cout << "proc: " << proc << std::endl;
        buffer.write(std::make_pair(BlockType(),GlobalId())); // write a default value
        if ((unsigned)i>=new2old_localindex.size()) return; // not our business since this index is not in the original range
        auto I = new2old_localindex[i];
        auto MEndIt = A[I].end();
        for (auto It = A[I].begin(); It!=MEndIt; ++It) // iterate over row and count entries
          if (old2new_localindex[It.index()]<new2old_localindex.size())
            {
              buffer.write(std::make_pair(*It,globalid[It.index()])); // yes globalid is the right thing since column index is in input index set
              //std::cout << "gather " << globalid[I] << " " << globalid[It.index()] << std::endl;
            }
      }

      template<class B>
      void scatter(B& buffer, int i, int count, int proc)
      {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        DataType pair;
        buffer.read(pair); // throw away the first one
        for (int k=1; k<count; k++)
          {
            buffer.read(pair);
            if (pis->exists(pair.second))
              {
                M.entry(i,pis->at(pair.second).local().local()) += pair.first;
                //std::cout << "scatter " << pair.second << std::endl;
              }
          }
      }
    };

    // data handle to insert matrix entries in the overlap region
    template<typename Matrix, typename ParallelIndexSet>
    class LambdaMultiMatrixConstructDataHandle
    {
      using LocalIndex = typename Matrix::size_type;
      using GlobalId = typename ParallelIndexSet::GlobalIndex;
      Matrix& M;       // the new matrix on the overlapping index set
      Matrix& M2;
      std::function<std::shared_ptr<Matrix>(int)> MatrixLambda; // the original input matrix
      std::shared_ptr<Matrix> A = nullptr;
      std::shared_ptr<ParallelIndexSet> pis;
      const std::vector<LocalIndex>& new2old_localindex;
      const std::vector<LocalIndex>& old2new_localindex;
      const std::vector<GlobalId>& globalid;

      using BlockType = typename Matrix::block_type;

      int previousRemoteProc = -1;
      int loadedMatrixProc = -1;

    public:
      // export type of items to be sent around
      using DataType = std::pair<BlockType,GlobalId>;

      // constructor
      LambdaMultiMatrixConstructDataHandle (Matrix& M_, Matrix& M2_, // the new matrix
                                      std::function<std::shared_ptr<Matrix>(int)> MatrixLambda_, // the input matrix in additive form
                                      std::shared_ptr<ParallelIndexSet> pis_,
                                      const std::vector<LocalIndex>& new2old_localindex_,
                                      const std::vector<LocalIndex>& old2new_localindex_,
                                      const std::vector<GlobalId>& globalid_)
      : M(M_), M2(M2_), MatrixLambda(MatrixLambda_), pis(pis_), new2old_localindex(new2old_localindex_), old2new_localindex(old2new_localindex_), globalid(globalid_)
      {}

      // enable variable size
      bool fixedsize()
      {
        return false;
      }

      // return size for given (local) index
      std::size_t size (int i, int proc) // i is a new local index
      {
        // Doesn't matter which matrix we send here (we assume identical pattern), so just make sure one is ready
        if (loadedMatrixProc == -1) {
            std::cout << "Assembling for " << proc << std::endl;
          A = MatrixLambda(proc);
          loadedMatrixProc = proc;
        }

        std::size_t count = 1; // always send one
        if ((unsigned)i>=new2old_localindex.size()) return count; // not our business since this index is not in the original range

        auto I = new2old_localindex[i]; // I is the local index in the input index set
        auto MEndIt = (*A)[I].end(); // Can simply use first matrix here, as we assume equal sparsity patterns
        for (auto It = (*A)[I].begin(); It!=MEndIt; ++It) // iterate over row and count entries
          if (old2new_localindex[It.index()]<new2old_localindex.size()) // column index is also non ghost
            count++;

        return count;
      }

      // gather to buffer
      template<class B>
      void gather(B& buffer, int i, int proc)
      {
        if (loadedMatrixProc != proc) {
          std::cout << "Assembling for " << proc << std::endl;
          A = MatrixLambda(proc);
          loadedMatrixProc = proc;
        }

        buffer.write(std::make_pair(BlockType(),GlobalId())); // write a default value
        if ((unsigned)i>=new2old_localindex.size()) return; // not our business since this index is not in the original range
        auto I = new2old_localindex[i];
        auto MEndIt = (*A)[I].end();
        for (auto It = (*A)[I].begin(); It!=MEndIt; ++It) // iterate over row and count entries
          if (old2new_localindex[It.index()]<new2old_localindex.size())
          {
            buffer.write(std::make_pair(*It,globalid[It.index()])); // yes globalid is the right thing since column index is in input index set
          }
      }

      template<class B>
      void scatter(B& buffer, int i, int count, int proc)
      {
        DataType pair;
        buffer.read(pair); // throw away the first one
        for (int k=1; k<count; k++)
        {
          buffer.read(pair);
          if (pis->exists(pair.second))
          {
            M.entry(i,pis->at(pair.second).local().local()) += pair.first;
            M2.entry(i,pis->at(pair.second).local().local()) += pair.first;
          }
        }
      }
    };

    // data handle to insert matrix entries in the overlap region
    template<typename Matrix, typename ParallelIndexSet>
    class MultiMatrixConstructDataHandle
    {
      using LocalIndex = typename Matrix::size_type;
      using GlobalId = typename ParallelIndexSet::GlobalIndex;
      Matrix& M;       // the new matrix on the overlapping index set
      Matrix& M2;
      std::map<int,std::shared_ptr<Matrix>> A; // the original input matrix
      std::shared_ptr<ParallelIndexSet> pis;
      const std::vector<LocalIndex>& new2old_localindex;
      const std::vector<LocalIndex>& old2new_localindex;
      const std::vector<GlobalId>& globalid;

      using BlockType = typename Matrix::block_type;

    public:
      // export type of items to be sent around
      using DataType = std::pair<BlockType,GlobalId>;

      // constructor
      MultiMatrixConstructDataHandle (Matrix& M_, Matrix& M2_, // the new matrix
                                 std::map<int,std::shared_ptr<Matrix>> A_, // the input matrix in additive form
                                 std::shared_ptr<ParallelIndexSet> pis_,
                                 const std::vector<LocalIndex>& new2old_localindex_,
                                 const std::vector<LocalIndex>& old2new_localindex_,
                                 const std::vector<GlobalId>& globalid_)
      : M(M_), M2(M2_), A(A_), pis(pis_), new2old_localindex(new2old_localindex_), old2new_localindex(old2new_localindex_), globalid(globalid_)
      {}

      // enable variable size
      bool fixedsize()
      {
        return false;
      }

      // return size for given (local) index
      std::size_t size (int i, int proc) // i is a new local index
      {
        std::size_t count = 1; // always send one
        if (i>=new2old_localindex.size()) return count; // not our business since this index is not in the original range
        auto I = new2old_localindex[i]; // I is the local index in the input index set
        auto MEndIt = (*A[proc])[I].end(); // Can simply use first matrix here, as we assume equal sparsity patterns
        for (auto It = (*A[proc])[I].begin(); It!=MEndIt; ++It) // iterate over row and count entries
          if (old2new_localindex[It.index()]<new2old_localindex.size()) // column index is also non ghost
            count++;
        //std::cout << "sending " << count << " entries for index " << i << " to rank " << proc << std::endl;
        return count;
      }

      // gather to buffer
      template<class B>
      void gather(B& buffer, int i, int proc)
      {
        //std::size_t count = 1; // always send one
        //std::cout << "proc: " << proc << std::endl;
        buffer.write(std::make_pair(BlockType(),GlobalId())); // write a default value
        if (i>=new2old_localindex.size()) return; // not our business since this index is not in the original range
        auto I = new2old_localindex[i];
        auto MEndIt = (*A[proc])[I].end();
        for (auto It = (*A[proc])[I].begin(); It!=MEndIt; ++It) // iterate over row and count entries
          if (old2new_localindex[It.index()]<new2old_localindex.size())
          {
            buffer.write(std::make_pair(*It,globalid[It.index()])); // yes globalid is the right thing since column index is in input index set
            //count++;
            //std::cout << "gather " << globalid[I] << " " << globalid[It.index()] << std::endl;
          }
        //std::cout << "sent " << count << " entries for index " << i << " to rank " << proc << std::endl;
      }

      template<class B>
      void scatter(B& buffer, int i, int count, int proc)
      {
        DataType pair;
        buffer.read(pair); // throw away the first one
        for (int k=1; k<count; k++)
        {
          buffer.read(pair);
          if (pis->exists(pair.second))
          {
            M.entry(i,pis->at(pair.second).local().local()) += pair.first;
            M2.entry(i,pis->at(pair.second).local().local()) += pair.first;
            //std::cout << "scatter " << pair.second << std::endl;
          }
        }
      }
    };

    /**
     * \brief communication data handle for finding out neighbors
     *
     * Simple version that only finds out ranks with overlap in the border entities
     *
     */
    template<typename GV>
    class UpdateNeighborsDataHandle
      : public Dune::CommDataHandleIF<UpdateNeighborsDataHandle<GV>,int>
    {
      const GV gv;
      const typename GV::IndexSet& indexset;
      std::vector<std::set<int>>& origins;

    public:
      typedef int DataType;

      UpdateNeighborsDataHandle (const GV& gv_, std::vector<std::set<int>>& origins_)
        : gv(gv_), indexset(gv.indexSet()), origins(origins_)
      {}

      bool contains (int dim, int codim) const
      {
        return (codim==dim);
      }

      bool fixedSize (int dim, int codim) const
      {
        return false;
      }

      template<class EntityType>
      size_t size (const EntityType& e) const
      {
        auto i = indexset.index(e);
        return origins[i].size();
      }

      template<class MessageBufferImp, class EntityType>
      void gather (MessageBufferImp& buff, const EntityType& e) const
      {
        auto& myset = origins[indexset.index(e)];
        for (auto it=myset.begin(); it!=myset.end(); ++it)
          buff.write(*it);
      }

      template<class MessageBufferImp, class EntityType>
      void scatter (MessageBufferImp& buff, const EntityType& e, size_t n)
      {
        auto& myset = origins[indexset.index(e)];
        DataType x;
        for (size_t i=0; i<n; ++i)
          {
            buff.read(x);
            myset.insert(x);
          }
      }
    };

  } // end namespace OverlapTools


    // Attribute type for communication and transfering partition type from the grid
  enum EPISAttribute {interior=0,border=1,overlap=2,front=3,ghost=4};


  /** A class extending a given index set by adding overlap with respect to a matrix graph
   *
   * \tparam CollectiveCommunication     Type collective communication
   * \tparam GlobalId                    Type representing globally unique ids
   * \tparam Matrix                      ISTL sparse matrix type
   */
  template<typename CollectiveCommunication, typename GlobalId, typename Matrix>
  class ExtendedParallelIndexSet
  {
  public:
    // make enum type for attributes public
    using AttributedLocalIndex = Dune::ParallelLocalIndex<EPISAttribute>;
    using ParallelIndexSet = Dune::ParallelIndexSet<GlobalId,AttributedLocalIndex,256>;
    using RemoteIndices = Dune::RemoteIndices<ParallelIndexSet>;
    using LocalIndex = typename Matrix::size_type;

  private:
    // some types we need
    using FieldType = typename Matrix::field_type;

    // local data
    const CollectiveCommunication& comm; // collective communication object
    int rank;     // my rank
    int p;        // number of subdomains in total
    int overlapsize;// overlap in terms of graph distance
    size_t N_orig; // number of degrees of freedom in the input
    size_t N_prec; // number of degrees of freedom of vectors used in preconditioner
    std::vector<LocalIndex> new2old_localindex; // maps new local index to a local index in the original input
    std::vector<LocalIndex> old2new_localindex; // maps local index from input set to its new value or the value N_orig if it is not contained
    std::shared_ptr<ParallelIndexSet> pis; // newly built up index set
    std::shared_ptr<RemoteIndices> si; // shared index information

    // data handle that constructs a set of new indices on each rank
    // - the parallel index set is not modified
    // - this works on the input matrix, so we must care about ghosts
    // - after completion the new indices can be read off
    class ExtenderDataHandle
    {
      const int rank;
      const ParallelIndexSet& pis;
      const Matrix& A;
      const std::vector<LocalIndex>& n2o;
      const std::vector<LocalIndex>& o2n;
      std::vector<GlobalId> globalid;
      std::set<GlobalId> newindices;

    public:
      // export type of items to be sent around
      typedef GlobalId DataType;

      // constructor
      ExtenderDataHandle (int rank_,
                          const ParallelIndexSet& pis_,
                          const Matrix& A_,
                          const std::vector<LocalIndex>& n2o_,
                          const std::vector<LocalIndex>& o2n_)
        : rank(rank_), pis(pis_), A(A_), n2o(n2o_), o2n(o2n_), globalid(pis.size())
      {
        // build a map from local index to global id
        auto endit = pis.end();
        for (auto it=pis.begin(); it!=endit; ++it)
          globalid[it->local().local()] = it->global();
      }

      // enable variable size
      bool fixedsize()
      {
        return false; // we have variable size
      }

      // return size for given (local) index
      std::size_t size (int i)
      {
        // since local indices are always appended at the end (i<n2o.size())
        // indicates that this local index was already in the initial set (including the ghosts)
        // then we have part of the graph and send information about our part of the row
        std::size_t count=1;
        if ((unsigned)i<n2o.size()) // i is a new local index since it comes from the parallel index set
          {
            auto I = n2o[i]; // I is the corresponding index in the original local index set of A
            auto cIt = A[I].begin();
            auto cEndIt = A[I].end();
            for (; cIt!=cEndIt; ++cIt)
              if (cIt.index()<o2n.size()) {
                if (o2n[cIt.index()]<n2o.size()) // this indicates a valid index
                  count++;
              } else std::cout << "BING THIS SHOULD NOT HAPPEN" << std::endl;
          }
        return count;
      }

      // gather to buffer
      template<class B>
      void gather (B& buffer, int i)
      {
        buffer.write(GlobalId()); // write dummy since we want message length >= 1
        if ((unsigned)i<n2o.size()) // i is a new local index since it comes from the parallel index set
          {
            auto I = n2o[i]; // I is the corresponding index in the original local index set of A
            auto cIt = A[I].begin();
            auto cEndIt = A[I].end();
            for (; cIt!=cEndIt; ++cIt)
              if (o2n[cIt.index()]<n2o.size()) // this indicates a valid index
                {
                  buffer.write(globalid[o2n[cIt.index()]]);
                  //std::cout << rank << ": gather " << globalid[i] << "," << globalid[o2n[cIt.index()]] << std::endl;
                }
          }
      }

      template<class B>
      void scatter(B& buffer, int i, int size)
      {
        DataType x;
        buffer.read(x); // read dummy
        // now read size items from buffer
        for (int k=1; k<size; k++)
          {
            buffer.read(x); // here we receive a global index
            if (!pis.exists(x))
              {
                newindices.insert(x);
                // std::cout << rank << ": scatter " << globalid[i] << " -> " << x << std::endl;
              }
          }
      }

      // give out the result
      std::vector<GlobalId> getNewIndices ()
      {
        // convert set to std::vector
        std::vector<GlobalId> temp;
        for (auto it=newindices.begin(); it!=newindices.end(); ++it)
          temp.push_back(*it);
        return temp;
      }
    };

  public:

    //! \brief Constructor.
    ExtendedParallelIndexSet (const CollectiveCommunication& comm_,    // collective communication object
                              const std::vector<int>& neighbors,       // ranks with whom we possibly ever share degrees of freedom
                              const Matrix& A,                         // matrix assembled on interior elements and with dirichlet rows replaced
                              const std::vector<EPISAttribute>& partitiontype,   // vector giving partitiontype for each degree of freedom
                              const std::vector<GlobalId>& globalid,   // vector mapping local index to globally unique id
                              int overlap_,                            // layers of overlap to add
                              bool verbose=false)                      // be verbose if true

      : comm(comm_), rank(comm.rank()), p(comm.size()), overlapsize(overlap_), N_orig(A.N()), old2new_localindex(A.N())
    {
      // compute dist of every dof to the boundary (or border dofs) up to a certain depth
      // this is needed later to make all indices public that are potentially in the overlap
      std::vector<int> dist(N_orig);
      for (size_t i=0; i<N_orig; i++)
        dist[i] = (partitiontype[i]==EPISAttribute::border) ? 0 : (overlapsize+1);
      for (int round=0; round<overlapsize; round++)
        for (size_t i=0; i<N_orig; i++)
          {
            auto cIt = A[i].begin();
            auto cEndIt = A[i].end();
            for (; cIt!=cEndIt; ++cIt)
              dist[i] = std::min(dist[i],dist[cIt.index()]+1);
          }

      // now we can set up a parallel index which is initialized as follows:
      // - it contains only interior and border dofs, ghost dofs are skipped
      // - a new local index is built up so that these are consecutive
      // - we need a map to transfer between the the new local index without ghosts and the original local index including ghosts
      //   -- new2old_localindex takes a new local index and gives an old local index
      //   -- old2new_localindex does the reverse map; this is only needed in setup phase when copying the input matrix
      // - only indices which have dist <= overlapsize are public and may be shared with indices
      pis = std::shared_ptr<ParallelIndexSet>(new ParallelIndexSet());
      pis->beginResize();
      for (size_t i=0; i<N_orig; i++)
        {
          if (partitiontype[i]==EPISAttribute::ghost)
            {
              old2new_localindex[i] = N_orig; // this index will be out of bounds to indicate that it was a ghost
              continue; // skip ghosts
            }
          EPISAttribute attr=partitiontype[i];
          bool pub = (partitiontype[i]==EPISAttribute::overlap)||(partitiontype[i]==EPISAttribute::front)||(dist[i]<=overlapsize);
          new2old_localindex.push_back(i);
          old2new_localindex[i] = new2old_localindex.size()-1;
          pis->add(globalid[i],AttributedLocalIndex(new2old_localindex.size()-1,attr,pub));
        }
      pis->endResize();
      if (verbose)
        std::cout << "rank["<<rank<<"] newsize="<<new2old_localindex.size()<<" oldsize="<<N_orig<<std::endl;

      // now let us extend the index set
      // - each round adds one layer of overlap
      for (int round=0; round<overlapsize; round++)
        {
          // build remote indices information
          RemoteIndices tempsi(*pis,*pis,comm,neighbors);
          // comm.barrier();
          // std::cout << rank << ": before building remote indices" << std::endl;
          // comm.barrier();
          tempsi.template rebuild<false>(); // find shared indices

          // build communication interface
          Dune::AllSet<EPISAttribute> allAttribute;
          Dune::Interface tempinterface;
          tempinterface.build(tempsi,allAttribute,allAttribute);

          // find new indices
          ExtenderDataHandle extdh(rank,*pis,A,new2old_localindex,old2new_localindex);
          DuneWithRank::VariableSizeCommunicator<> varcommunicator(tempinterface);
          // comm.barrier();
          // std::cout << rank << ": before finding new indices" << std::endl;
          // comm.barrier();
          varcommunicator.forward(extdh);
          // comm.barrier();
          // std::cout << rank << ": after finding new indices" << std::endl;
          // comm.barrier();
          auto newindices = extdh.getNewIndices();
          // std::cout << "rank " << rank << " has " << newindices.size() << " new indices" << std::endl;

          // and lets enlarge the index set ...
          LocalIndex localindex = pis->size(); // append new local indices at the end
          pis->beginResize();
          for (auto it=newindices.begin(); it!=newindices.end(); ++it)
            {
              //std::cout << "Proc [" << helper.rank() << "] id=" << globalidset.id(v) << " index=" << indexset.index(v) << std::endl;
              pis->add(*it,AttributedLocalIndex(localindex++,overlap,true));
            }
          pis->endResize();

          if (verbose)
            for (int i=0; i<p; i++)
              {
                comm.barrier();
                if (rank==i) {
                  std::cout <<rank<<": round " << round << " : additional indices="<<newindices.size()<<" new size now="<<pis->size()<<std::endl;
                }
              }
        }
      N_prec = pis->size(); // this now the new size after adding overlap
      //std::cout <<rank<<": input size=" << N_orig <<" new size="<<pis->size()<<std::endl;

      // now find out about shared dofs
      si = std::shared_ptr<RemoteIndices>(new RemoteIndices(*pis,*pis,comm,neighbors));
      si->template rebuild<false>(); // exchanges all sharing information
    }

    //! get number of indices in input index set taken from A.N()
    size_t originalSize () const
    {
      return N_orig;
    }

    //! get number of indices after adding overlap
    size_t extendedSize () const
    {
      return N_prec;
    }

    int overlapSize () const
    {
      return overlapsize;
    }

    //! get the ParallelIndexSet
    std::shared_ptr<ParallelIndexSet> parallelIndexSet () const
    {
      return pis;
    }

    //! get the RemoteIndices
    std::shared_ptr<RemoteIndices> remoteIndices () const
    {
      return si;
    }

    //! returns a map mapping a local index from the extended index set to the original index set
    // This only works for those degrees of freedom that are in both sets
    const std::vector<LocalIndex>& extendedToOriginalLocalIndex () const
    {
      return new2old_localindex;
    }

    //! returns a map mapping a local index from the original index set to the extended index set
    //  This map works on all indices in the original set. If the index is not present in the extended set
    //  (i.e. it corresponded to a ghost dof) then originalSize() is returned by the map.
    const std::vector<LocalIndex>& originalToExtendedLocalIndex () const
    {
      return old2new_localindex;
    }

    const CollectiveCommunication& collectiveCommunication () const
    {
      return comm;
    }
  };


  template<typename ExtendedParallelIndexSet, typename Matrix, typename GlobalId>
  std::shared_ptr<Matrix> copyAndExtendMatrix (const ExtendedParallelIndexSet& epis, const Matrix& A,
                                               int avg, const std::vector<GlobalId>& globalid) {
    return copyAndExtendMatrix<ExtendedParallelIndexSet, Matrix, GlobalId>(epis, A, A, avg, globalid);
  }


  /** Take a matrix on the original index set and return one on the extended index set
   * The input matrix is in additive form assembled on the interior elements only.
   * The returned matrix is in consistent form with Dirichlet boundary conditions.
   */
  template<typename ExtendedParallelIndexSet, typename Matrix, typename GlobalId>
  std::shared_ptr<Matrix> copyAndExtendMatrix (const ExtendedParallelIndexSet& epis, const Matrix& A_local, const Matrix& A,
                                               int avg, const std::vector<GlobalId>& globalid)
  {
    // build up communication interface using all atributes (thankfully ghosts are out of the way :-))
    Dune::AllSet<EPISAttribute> allAttribute;
    Dune::Interface allinterface;
    allinterface.build(*epis.remoteIndices(),allAttribute,allAttribute); // all to all communication
    DuneWithRank::VariableSizeCommunicator<> varcommunicator(allinterface);
    //std::cout << "interface is built" << std::endl;

    // make a copy M of the matrix A_local where rows/columns corresponding to ghosts are omitted
    //std::cout << rank << ": copy matrix" << std::endl;
    auto M = std::shared_ptr<Matrix>(new Matrix(epis.extendedSize(),epis.extendedSize(),avg,0.1,Matrix::implicit));
    auto new2old_localindex = epis.extendedToOriginalLocalIndex();
    auto old2new_localindex = epis.originalToExtendedLocalIndex();
    for (auto rIt=A_local.begin(); rIt!=A_local.end(); ++rIt) // loop over entries of A
      for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
        {
          auto i = old2new_localindex[rIt.index()];
          auto j = old2new_localindex[cIt.index()];
          if ( i<new2old_localindex.size() && j<new2old_localindex.size() )
            M->entry(i,j) = *cIt;
        }
    // do not compress yet; we need to add entries from the other processors!
    // now we need to make M a submatrix centered on the diagonal
    // add entries of M from other processors
    //std::cout << rank << ": communicate matrix" << std::endl;
    using ParallelIndexSet = typename ExtendedParallelIndexSet::ParallelIndexSet;
    OverlapTools::MatrixConstructDataHandle<Matrix,ParallelIndexSet> matconsdh(*M,A,epis.parallelIndexSet(),new2old_localindex,
                                                                               old2new_localindex,globalid);
    varcommunicator.forward(matconsdh);
    //std::cout << rank << ": matrix communication finished" << std::endl;
    M->compress();
    //std::cout << rank << ": matrix done" << std::endl;
    return M;
  }


  template<typename ExtendedParallelIndexSet, typename Matrix, typename GlobalId>
  std::pair<std::shared_ptr<Matrix>,std::shared_ptr<Matrix>> lambdaMultiCopyAndExtendMatrix (const ExtendedParallelIndexSet& epis, const Matrix& A_local, const Matrix& A2_local, std::function<std::shared_ptr<Matrix>(int)> MatrixLambda,
                                               int avg, const std::vector<GlobalId>& globalid)
  {
    // build up communication interface using all atributes (thankfully ghosts are out of the way :-))
    Dune::AllSet<EPISAttribute> allAttribute;
    Dune::Interface allinterface;
    allinterface.build(*epis.remoteIndices(),allAttribute,allAttribute); // all to all communication
    DuneWithRank::VariableSizeCommunicator<> varcommunicator(allinterface);
    //std::cout << "interface is built" << std::endl;

    // make a copy M of the matrix A_local where rows/columns corresponding to ghosts are omitted
    //std::cout << rank << ": copy matrix" << std::endl;
    std::shared_ptr<Matrix> M = std::shared_ptr<Matrix>(new Matrix(epis.extendedSize(),epis.extendedSize(),avg,0.1,Matrix::implicit));
    std::shared_ptr<Matrix> M2 = std::shared_ptr<Matrix>(new Matrix(epis.extendedSize(),epis.extendedSize(),avg,0.1,Matrix::implicit));
    auto new2old_localindex = epis.extendedToOriginalLocalIndex();
    auto old2new_localindex = epis.originalToExtendedLocalIndex();
    for (auto rIt=A_local.begin(); rIt!=A_local.end(); ++rIt) // loop over entries of A
    {
      for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
      {
        auto i = old2new_localindex[rIt.index()];
        auto j = old2new_localindex[cIt.index()];
        if ( i<new2old_localindex.size() && j<new2old_localindex.size() ) {
          M->entry(i,j) = *cIt;
        }
      }
    }
    for (auto rIt=A2_local.begin(); rIt!=A2_local.end(); ++rIt) // loop over entries of A
    {
      for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
      {
        auto i = old2new_localindex[rIt.index()];
        auto j = old2new_localindex[cIt.index()];
        if ( i<new2old_localindex.size() && j<new2old_localindex.size() ) {
          M2->entry(i,j) = *cIt;
        }
      }
    }
    // do not compress yet; we need to add entries from the other processors!
    // now we need to make M a submatrix centered on the diagonal
    // add entries of M from other processors
    //std::cout << rank << ": communicate matrix" << std::endl;
    using ParallelIndexSet = typename ExtendedParallelIndexSet::ParallelIndexSet;
    OverlapTools::LambdaMultiMatrixConstructDataHandle<Matrix,ParallelIndexSet> matconsdh(*M,*M2,MatrixLambda,epis.parallelIndexSet(),new2old_localindex,
                                                                                old2new_localindex,globalid);
    varcommunicator.forward(matconsdh);
    //std::cout << rank << ": matrix communication finished" << std::endl;
    auto stats = M->compress();
    stats = M2->compress();
    //std::cout << rank << ": matrix done" << std::endl;
    return std::make_pair(M,M2);
  }

  /** Take a matrix on the original index set and return one on the extended index set
   * The input matrix is in additive form assembled on the interior elements only.
   * The returned matrix is in consistent form with Dirichlet boundary conditions.
   */
  template<typename ExtendedParallelIndexSet, typename Matrix, typename GlobalId>
  std::pair<std::shared_ptr<Matrix>,std::shared_ptr<Matrix>> multiCopyAndExtendMatrix (const ExtendedParallelIndexSet& epis, const Matrix& A_local, const Matrix& A2_local, std::map<int,std::shared_ptr<Matrix>> A,
                                               int avg, const std::vector<GlobalId>& globalid)
  {
    // build up communication interface using all atributes (thankfully ghosts are out of the way :-))
    Dune::AllSet<EPISAttribute> allAttribute;
    Dune::Interface allinterface;
    allinterface.build(*epis.remoteIndices(),allAttribute,allAttribute); // all to all communication
    DuneWithRank::VariableSizeCommunicator<> varcommunicator(allinterface);
    //std::cout << "interface is built" << std::endl;

    // make a copy M of the matrix A_local where rows/columns corresponding to ghosts are omitted
    //std::cout << rank << ": copy matrix" << std::endl;
    std::shared_ptr<Matrix> M = std::shared_ptr<Matrix>(new Matrix(epis.extendedSize(),epis.extendedSize(),avg,0.1,Matrix::implicit));
    std::shared_ptr<Matrix> M2 = std::shared_ptr<Matrix>(new Matrix(epis.extendedSize(),epis.extendedSize(),avg,0.1,Matrix::implicit));
    auto new2old_localindex = epis.extendedToOriginalLocalIndex();
    auto old2new_localindex = epis.originalToExtendedLocalIndex();
    for (auto rIt=A_local.begin(); rIt!=A_local.end(); ++rIt) // loop over entries of A
    {
      for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
      {
        auto i = old2new_localindex[rIt.index()];
        auto j = old2new_localindex[cIt.index()];
        if ( i<new2old_localindex.size() && j<new2old_localindex.size() ) {
          M->entry(i,j) = *cIt;
        }
      }
    }
    for (auto rIt=A2_local.begin(); rIt!=A2_local.end(); ++rIt) // loop over entries of A
    {
      for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
      {
        auto i = old2new_localindex[rIt.index()];
        auto j = old2new_localindex[cIt.index()];
        if ( i<new2old_localindex.size() && j<new2old_localindex.size() ) {
          M2->entry(i,j) = *cIt;
        }
      }
    }
    // do not compress yet; we need to add entries from the other processors!
    // now we need to make M a submatrix centered on the diagonal
    // add entries of M from other processors
    //std::cout << rank << ": communicate matrix" << std::endl;
    using ParallelIndexSet = typename ExtendedParallelIndexSet::ParallelIndexSet;
    OverlapTools::MultiMatrixConstructDataHandle<Matrix,ParallelIndexSet> matconsdh(*M,*M2,A,epis.parallelIndexSet(),new2old_localindex,
                                                                                old2new_localindex,globalid);
    varcommunicator.forward(matconsdh);
    //std::cout << rank << ": matrix communication finished" << std::endl;
    auto stats = M->compress();
    stats = M2->compress();
    //std::cout << rank << ": matrix done" << std::endl;
    return std::make_pair(M,M2);
  }

  /** Take the matrix on the overlapping index set and subtract all entries from the diagonal that are not stored in this rank
   */
  template<typename ExtendedParallelIndexSet, typename Matrix, typename ScalarVector>
  void addNonlocalEntriesToDiagonal (const ExtendedParallelIndexSet& epis, Matrix& M, const ScalarVector& owner)
  {
    // build up communication interface using all atributes (thankfully ghosts are out of the way :-))
    Dune::AllSet<EPISAttribute> allAttribute;
    Dune::Interface allinterface;
    allinterface.build(*epis.remoteIndices(),allAttribute,allAttribute); // all to all communication
    DuneWithRank::VariableSizeCommunicator<> varcommunicator(allinterface);

    using ParallelIndexSet = typename ExtendedParallelIndexSet::ParallelIndexSet;
    OverlapTools::AddNonlocalEntriesToDiagonal<Matrix,ParallelIndexSet,ScalarVector> adddh(M,*epis.parallelIndexSet(),owner);
    varcommunicator.forward(adddh);
  }


  template<typename ExtendedParallelIndexSet, typename Matrix, typename Vector>
  std::shared_ptr<Vector>
  makePartitionOfUnity (const ExtendedParallelIndexSet& epis, const Matrix& M, int overlapsize)
  {
    using FieldType = typename Matrix::field_type;
    using ScalarVector = Dune::BlockVector<Dune::FieldVector<FieldType,1>>;
    using ParallelIndexSet = typename ExtendedParallelIndexSet::ParallelIndexSet;

    // communicator
    Dune::AllSet<EPISAttribute> allAttribute;
    Dune::Interface allinterface;
    allinterface.build(*epis.remoteIndices(),allAttribute,allAttribute); // all to all communication
    DuneWithRank::VariableSizeCommunicator<> varcommunicator(allinterface);
    Dune::BufferedCommunicator communicator;
    communicator.build<ScalarVector>(allinterface);

    // now lets determine an improved partition of unity ....
    // first find entries which have nonlocal edges; these are on the boundary
    ScalarVector distance(M.N());
    OverlapTools::CountNonlocalEdgesDataHandle<ParallelIndexSet,Matrix,ScalarVector> nonlocaledgesdh(M,*epis.parallelIndexSet(),distance);
    varcommunicator.forward(nonlocaledgesdh);

    // now compute distance of each vertex to the boundary up to distance 2*overlapsize
    for (typename ScalarVector::size_type i=0; i<distance.N(); i++)
      if (distance[i]!=0.0)
        distance[i] = 0.0; // this is a boundary dof
      else
        distance[i] = 2*overlapsize+1.0;

    for (int round=0; round<2*overlapsize; round++)
      for (typename ScalarVector::size_type i=0; i<distance.N(); i++)
        {
          auto cIt = M[i].begin();
          auto cEndIt = M[i].end();
          for (; cIt!=cEndIt; ++cIt)
            distance[i] = std::min(distance[i],distance[cIt.index()]+1.0);
        }

    // now we may compute the partition of unity as in the stability proof of additive Schwarz
    ScalarVector sumdistance(distance);
    for (typename ScalarVector::size_type i=0; i<sumdistance.N(); i++)
      //sumdistance[i] += 1.0; // add 1 as actually the first vertex is already inside
      sumdistance[i] += .0; // add 1 as actually the first vertex is already inside
    communicator.forward<OverlapTools::AddGatherScatter<ScalarVector>>(sumdistance,sumdistance);
    auto pu = std::shared_ptr<Vector>(new Vector(M.N()));
    for (typename Vector::size_type i=0; i<pu->N(); i++){
      (*pu)[i] = (distance[i]+.0)/sumdistance[i];
    }
    return pu;
  }

  template<typename ExtendedParallelIndexSet, typename Matrix>
  std::shared_ptr<Dune::BlockVector<Dune::FieldVector<typename Matrix::field_type,1>>>
  makeOwner (const ExtendedParallelIndexSet& epis, const Matrix& M, int overlapsize)
  {
    using FieldType = typename Matrix::field_type;
    using ScalarVector = Dune::BlockVector<Dune::FieldVector<FieldType,1>>;
    using ParallelIndexSet = typename ExtendedParallelIndexSet::ParallelIndexSet;

    // communicator
    Dune::AllSet<EPISAttribute> allAttribute;
    Dune::Interface allinterface;
    allinterface.build(*epis.remoteIndices(),allAttribute,allAttribute); // all to all communication
    DuneWithRank::VariableSizeCommunicator<> varcommunicator(allinterface);

    DuneWithRank::BufferedCommunicator communicator;
    communicator.build<ScalarVector>(allinterface);

    // now lets determine an improved partition of unity ....
    // first find entries which have nonlocal edges; these are on the boundary
    ScalarVector distance(M.N());
    OverlapTools::CountNonlocalEdgesDataHandle<ParallelIndexSet,Matrix,ScalarVector> nonlocaledgesdh(M,*epis.parallelIndexSet(),distance);
    varcommunicator.forward(nonlocaledgesdh);
    // now compute distance of each vertex to the boundary up to distance 2*overlapsize
    for (typename ScalarVector::size_type i=0; i<distance.N(); i++)
      if (distance[i]!=0.0)
        distance[i] = 0.0; // this is a boundary dof
      else
        distance[i] = 2*overlapsize+1.0;
    for (int round=0; round<2*overlapsize; round++)
      for (typename ScalarVector::size_type i=0; i<distance.N(); i++)
        {
          auto cIt = M[i].begin();
          auto cEndIt = M[i].end();
          for (; cIt!=cEndIt; ++cIt)
            distance[i] = std::min(distance[i],distance[cIt.index()]+1.0);
        }

    // build up an owner mask. Owners must be at least distance 1 from the boundary
    auto owner = std::shared_ptr<ScalarVector>(new ScalarVector(M.N()));
    OverlapTools::MakeOwnerDataHandle ownerdh(distance,epis.collectiveCommunication().rank(),epis.collectiveCommunication().size(),*owner);
    varcommunicator.forward(ownerdh);

    return owner;
  }
} // namespace Dune

// The following depends at least on Dune Grid
namespace Dune {
  namespace PDELab {

    /** find all ranks with whom we possible share neighbors in the overlap given
     */
    template<typename GV>
    std::vector<int> findNeighboringRanks (const GV& gv, int overlap)
    {
      auto& indexset = gv.indexSet();
      const int dimension = GV::dimension;

      // find out true set of neighbors in distance overlapsize
      std::vector<std::set<int>> origins1(indexset.size(dimension)), origins2(indexset.size(dimension));
      for (int i=0; i<origins1.size(); ++i) origins1[i].insert(gv.comm().rank()); // I know my indices initially
      OverlapTools::UpdateNeighborsDataHandle<GV> unbdhx(gv,origins1);
      if (gv.comm().size()>1)
        gv.communicate(unbdhx,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication); // do one update with neighbors
      std::vector<std::set<int>>* current = &origins1; // pointers for double buffering
      std::vector<std::set<int>>* update = &origins2;
      for (int round=0; round<2*overlap; ++round)
        {
          // propagate ranks by one layer of elements
          *update = *current;
          for (const auto& e : elements(gv,Dune::Partitions::interior))
            for (int I=0; I<e.geometry().corners(); I++)
              {
                auto i = indexset.subIndex(e,I,dimension);
                for (int J=0; J<e.geometry().corners(); J++)
                  if (J!=I)
                    {
                      auto j = indexset.subIndex(e,J,dimension);
                      (*update)[i].insert((*current)[j].begin(),(*current)[j].end()); // union of two sets
                    }
              }

          // update border vertices
          OverlapTools::UpdateNeighborsDataHandle<GV> unbdh(gv,*update);
          if (gv.comm().size()>1)
            gv.communicate(unbdh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication); // do one update with neighbors

          // swap pointers
          auto temp = update; update = current; current = temp;
        }

      // now unite all sets
      std::set<int> allmyneighborsset;
      for (int i=0; i<(*current).size(); ++i)
        allmyneighborsset.insert((*current)[i].begin(),(*current)[i].end());
      // ... and convert to vector
      std::vector<int> allmyneighborsvec;
      for (auto it=allmyneighborsset.begin(); it!=allmyneighborsset.end(); ++it)
        if (*it!=gv.comm().rank())
          allmyneighborsvec.push_back(*it);

      return allmyneighborsvec;
    }

  } // namespace PDELab
} // namespace Dune



namespace Dune {
  template<typename GridView, typename Vector, typename Matrix>
  class NonoverlappingOverlapAdapter {

    using GlobalId = typename GridView::Grid::GlobalIdSet::IdType;
    using CollectiveCommunication = typename GridView::CollectiveCommunication;
    using EPIS = Dune::ExtendedParallelIndexSet<CollectiveCommunication,GlobalId,Matrix>;
    using LocalIndex = typename Vector::size_type;
    using FieldType = typename Vector::field_type;
    using ScalarVector = Dune::BlockVector<Dune::FieldVector<FieldType,1>>;

    using Attribute = EPISAttribute;
    using AttributedLocalIndex = Dune::ParallelLocalIndex<Attribute>;
    using ParallelIndexSet = Dune::ParallelIndexSet<GlobalId,AttributedLocalIndex,256>;
    using RemoteIndices = Dune::RemoteIndices<ParallelIndexSet>;

  public:
    NonoverlappingOverlapAdapter(const GridView& gv,       // ranks with whom we share degrees of freedom
                                 const Matrix& A,          // Matrix used to determine matrix graph
                                 int avg,                                // average number of matrix entries
                                 int overlap)
    : gv_(gv),
    //A_(A),
    avg_(avg),
    globalid_(buildGlobalIdVector(gv)), // TODO: Avoid copy?
    overlap_(overlap),
    epis_(gv.comm(),Dune::PDELab::findNeighboringRanks(gv,overlap_),A,buildPartitionTypeVector(gv),globalid_,overlap_,false)
    {
      new2old_localindex_ = epis_.extendedToOriginalLocalIndex(); // TODO: avoid copy?
    }

    size_t getExtendedSize() const {
      return epis_.parallelIndexSet()->size();
    }
    std::shared_ptr<RemoteIndices> getRemoteIndices() const {
      return epis_.remoteIndices();
    }
    const EPIS& getEpis() const {
      return epis_;
    }
    std::vector<int> findNeighboringRanks() {
      return Dune::PDELab::findNeighboringRanks(gv_,overlap_);
    }

    std::shared_ptr<Matrix> extendMatrix(const Matrix& A) {
      return Dune::copyAndExtendMatrix(epis_,A,avg_,globalid_);
    }

    std::shared_ptr<Matrix> extendMatrix(const Matrix& A_local, const Matrix& A) {
      return Dune::copyAndExtendMatrix(epis_,A_local,A,avg_,globalid_);
    }

    std::pair<std::shared_ptr<Matrix>,std::shared_ptr<Matrix>> multiExtendMatrix(const Matrix& A_local, const Matrix& A2_local, std::map<int,std::shared_ptr<Matrix>> A) {
      return Dune::multiCopyAndExtendMatrix(epis_,A_local,A2_local,A,avg_,globalid_);
    }

    std::pair<std::shared_ptr<Matrix>,std::shared_ptr<Matrix>> lambdaMultiExtendMatrix(const Matrix& A_local, const Matrix& A2_local, std::function<std::shared_ptr<Matrix>(int)> MatrixLambda) {
      return Dune::lambdaMultiCopyAndExtendMatrix(epis_,A_local,A2_local,MatrixLambda,avg_,globalid_);
    }

    const GridView& gridView() const {
      return gv_;
    }

    void extendVector(const Vector& restricted, Vector& extended) const {
      extended = 0.0;
      for (typename Vector::size_type i=0; i<new2old_localindex_.size(); i++)
        extended[i] = restricted[new2old_localindex_[i]];
    }
    void restrictVector(const Vector& extended, Vector& restricted) const {
      for (typename ScalarVector::size_type i=0; i<new2old_localindex_.size(); i++)
        restricted[new2old_localindex_[i]] = extended[i];
    }
    std::vector<LocalIndex> get_old2new_localindex() {
      return epis_.originalToExtendedLocalIndex();
    }

    std::vector<LocalIndex> get_new2old_localindex() {
      return new2old_localindex_;
    }

    std::vector<typename GridView::Grid::GlobalIdSet::IdType> get_globalid() {
      auto& indexset = gv_.indexSet();
      auto& globalidset = gv_.grid().globalIdSet();
      using GlobalID = typename GridView::Grid::GlobalIdSet::IdType;
      std::vector<GlobalID> globalid(indexset.size(GridView::Grid::dimension));
      for (const auto& v : vertices(gv_,Dune::Partitions::all))
        globalid[indexset.index(v)] = globalidset.id(v);
      return globalid;
    }

  private:

    std::vector<EPISAttribute> buildPartitionTypeVector(const GridView& gv) {
      auto& indexset = gv.indexSet();
      std::vector<Dune::EPISAttribute> partitiontype(indexset.size(GridView::Grid::dimension));

      for (const auto& v : vertices(gv,Dune::Partitions::all))
      {
        if (v.partitionType()==Dune::InteriorEntity) partitiontype[indexset.index(v)] = Dune::EPISAttribute::interior;
        if (v.partitionType()==Dune::BorderEntity) partitiontype[indexset.index(v)] = Dune::EPISAttribute::border;
        if (v.partitionType()==Dune::OverlapEntity) partitiontype[indexset.index(v)] = Dune::EPISAttribute::overlap;
        if (v.partitionType()==Dune::GhostEntity) partitiontype[indexset.index(v)] = Dune::EPISAttribute::ghost;
      }
      return partitiontype;
    }

    //template<typename GridView>
    std::vector<typename GridView::Grid::GlobalIdSet::IdType> buildGlobalIdVector(const GridView& gv) {
      auto& indexset = gv.indexSet();
      auto& globalidset = gv.grid().globalIdSet();
      using GlobalID = typename GridView::Grid::GlobalIdSet::IdType;
      std::vector<GlobalID> globalid(indexset.size(GridView::Grid::dimension));
      for (const auto& v : vertices(gv,Dune::Partitions::all))
        globalid[indexset.index(v)] = globalidset.id(v);
      return globalid;
    }

    const GridView& gv_;
    const int avg_;
    const std::vector<GlobalId> globalid_;
    const int overlap_;
    EPIS epis_;
  public:
    std::vector<LocalIndex> new2old_localindex_;
  };

  template<typename GridView, typename Matrix, typename Vector>
  std::shared_ptr<Vector> makePartitionOfUnity(NonoverlappingOverlapAdapter<GridView, Vector, Matrix>& adapter, const Matrix& A) {
    using GlobalId = typename GridView::Grid::GlobalIdSet::IdType;
    using CollectiveCommunication = typename GridView::CollectiveCommunication;
    using EPIS = Dune::ExtendedParallelIndexSet<CollectiveCommunication,GlobalId,Matrix>;
    return Dune::makePartitionOfUnity<EPIS, Matrix, Vector>(adapter.getEpis(), A, adapter.getEpis().overlapSize());
  }
}


#endif // DUNE_PDELAB_EXTEND_OVERLAP_TOOLS_HH
