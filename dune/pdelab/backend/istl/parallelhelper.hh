// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_PARALLELHELPER_HH
#define DUNE_PDELAB_BACKEND_ISTL_PARALLELHELPER_HH

// this is here for backwards compatibility and deprecation warnings, remove after 2.5.0
#include "ensureistlinclude.hh"

#include <limits>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/typetraits.hh>

#if HAVE_UG && PDELAB_SEQUENTIAL_UG
// We need the UGGrid declaration for the assertion
#include <dune/grid/uggrid.hh>
#endif

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

#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/istl/vector.hh>
#include <dune/pdelab/backend/istl/utility.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>

namespace Dune {
  namespace PDELab {
    namespace ISTL {

      //! \addtogroup Backend
      //! \ingroup PDELab
      //! \{


      //========================================================
      // A parallel helper class providing a nonoverlapping
      // decomposition of all degrees of freedom
      //========================================================

      template<typename GFS>
      class ParallelHelper
      {

        //! Type for storing rank values.
        typedef int RankIndex;

        //! Type used to store owner rank values of all DOFs.
        using RankVector = Dune::PDELab::Backend::Vector<GFS,RankIndex>;
        //! Type used to store ghost flag of all DOFs.
        using GhostVector = Dune::PDELab::Backend::Vector<GFS,bool>;

        //! ContainerIndex of the underlying GridFunctionSpace.
        typedef typename GFS::Ordering::Traits::ContainerIndex ContainerIndex;

      public:

        ParallelHelper (const GFS& gfs, int verbose = 1)
          : _gfs(gfs)
          , _rank(gfs.gridView().comm().rank())
          , _ranks(gfs,_rank)
          , _ghosts(gfs,false)
          , _verbose(verbose)
        {

          // Let's try to be clever and reduce the communication overhead by picking the smallest
          // possible communication interface depending on the overlap structure of the GFS.
          // FIXME: Switch to simple comparison as soon as dune-grid:1b3e83ec0 is reliably available!
          if (gfs.entitySet().partitions().value == Partitions::interiorBorder.value)
            {
              // The GFS only spans the interior and border partitions, so we can skip sending or
              // receiving anything else.
              _interiorBorder_all_interface = InteriorBorder_InteriorBorder_Interface;
              _all_all_interface = InteriorBorder_InteriorBorder_Interface;
            }
          else
            {
              // In general, we have to transmit more.
              _interiorBorder_all_interface = InteriorBorder_All_Interface;
              _all_all_interface = All_All_Interface;
            }

          if (_gfs.gridView().comm().size()>1)
            {

              // find out about ghosts
              //GFSDataHandle<GFS,GhostVector,GhostGatherScatter>
              GhostDataHandle<GFS,GhostVector>
                gdh(_gfs,_ghosts,false);
              _gfs.gridView().communicate(gdh,_interiorBorder_all_interface,Dune::ForwardCommunication);

              // create disjoint DOF partitioning
              //            GFSDataHandle<GFS,RankVector,DisjointPartitioningGatherScatter<RankIndex> >
              //  ibdh(_gfs,_ranks,DisjointPartitioningGatherScatter<RankIndex>(_rank));
              DisjointPartitioningDataHandle<GFS,RankVector> pdh(_gfs,_ranks);
              _gfs.gridView().communicate(pdh,_interiorBorder_all_interface,Dune::ForwardCommunication);

            }

        }

        //! Mask out all DOFs not owned by the current process with 0.
        template<typename X>
        void maskForeignDOFs(X& x) const
        {
          using Backend::native;
          // dispatch to implementation.
          maskForeignDOFs(ISTL::container_tag(native(x)),native(x),native(_ranks));
        }

      private:

        // Implementation for block vector; recursively masks blocks.
        template<typename X, typename Mask>
        void maskForeignDOFs(ISTL::tags::block_vector, X& x, const Mask& mask) const
        {
          typename Mask::const_iterator mask_it = mask.begin();
          for (typename X::iterator it = x.begin(),
                 end_it = x.end();
               it != end_it;
               ++it, ++mask_it)
            maskForeignDOFs(ISTL::container_tag(*it),*it,*mask_it);
        }

        // Implementation for field vector, iterates over entries and masks them individually.
        template<typename X, typename Mask>
        void maskForeignDOFs(ISTL::tags::field_vector, X& x, const Mask& mask) const
        {
          typename Mask::const_iterator mask_it = mask.begin();
          for (typename X::iterator it = x.begin(),
                 end_it = x.end();
               it != end_it;
               ++it, ++mask_it)
            *it = (*mask_it == _rank ? *it : typename X::field_type(0));
        }

      public:

        //! Tests whether the given index is owned by this process.
        bool owned(const ContainerIndex& i) const
        {
          return _ranks[i] == _rank;
        }

        //! Tests whether the given index belongs to a ghost DOF.
        bool isGhost(const ContainerIndex& i) const
        {
          return _ghosts[i];
        }

        //! Calculates the (rank-local) dot product of x and y on the disjoint partition defined by the helper.
        template<typename X, typename Y>
        typename PromotionTraits<
          typename X::field_type,
          typename Y::field_type
          >::PromotedType
        disjointDot(const X& x, const Y& y) const
        {
          using Backend::native;
          return disjointDot(ISTL::container_tag(native(x)),
                             native(x),
                             native(y),
                             native(_ranks)
                             );
        }

      private:

        // Implementation for BlockVector, collects the result of recursively
        // invoking the algorithm on the vector blocks.
        template<typename X, typename Y, typename Mask>
        typename PromotionTraits<
        typename X::field_type,
        typename Y::field_type
        >::PromotedType
        disjointDot(ISTL::tags::block_vector, const X& x, const Y& y, const Mask& mask) const
        {
          typedef typename PromotionTraits<
            typename X::field_type,
            typename Y::field_type
            >::PromotedType result_type;

          result_type r(0);

          typename Y::const_iterator y_it = y.begin();
          typename Mask::const_iterator mask_it = mask.begin();
          for (typename X::const_iterator x_it = x.begin(),
                 end_it = x.end();
               x_it != end_it;
               ++x_it, ++y_it,  ++mask_it)
            r += disjointDot(ISTL::container_tag(*x_it),*x_it,*y_it,*mask_it);

          return r;
        }

        // Implementation for FieldVector, iterates over the entries and calls Dune::dot() for DOFs
        // associated with the current rank.
        template<typename X, typename Y, typename Mask>
        typename PromotionTraits<
          typename X::field_type,
          typename Y::field_type
          >::PromotedType
        disjointDot(ISTL::tags::field_vector, const X& x, const Y& y, const Mask& mask) const
        {
          typedef typename PromotionTraits<
            typename X::field_type,
            typename Y::field_type
            >::PromotedType result_type;

          result_type r(0);

          typename Y::const_iterator y_it = y.begin();
          typename Mask::const_iterator mask_it = mask.begin();
          for (typename X::const_iterator x_it = x.begin(),
                 end_it = x.end();
               x_it != end_it;
               ++x_it, ++y_it, ++mask_it)
            r += (*mask_it == _rank ? Dune::dot(*x_it,*y_it) : result_type(0));

          return r;
        }

      public:

        //! Returns the MPI rank of this process.
        RankIndex rank() const
        {
          return _rank;
        }

#if HAVE_MPI

        //! Makes the matrix consistent and creates the parallel information for AMG.
        /**
         * This function accomplishes two things:
         *
         * 1. Makes the matrix consistent w.r.t. to the disjoint partitioning of the DOF space,
         *    i.e. aggregates matrix entries for border entries from neighboring ranks.
         *
         * 2. Sets up the parallel communication information for AMG.
         *
         * \warning  This function silenty assumes that the matrix only has a single level
         *           of blocking and will not work correctly otherwise. Also note that AMG
         *           will only work correctly for P1 discretisations.
         *
         * \param m  The PDELab matrix container.
         * \param c  The parallel information object providing index set, interfaces and
         *           communicators.
         */
        template<typename MatrixType, typename Comm>
        void createIndexSetAndProjectForAMG(MatrixType& m, Comm& c);

      private:

        // Checks whether a matrix block is owned by the current process. Used for the AMG
        // construction and thus assumes a single level of blocking and blocks with ownership
        // restricted to a single DOF.
        bool owned_for_amg(std::size_t i) const
        {
          return Backend::native(_ranks)[i][0] == _rank;
        }

#endif // HAVE_MPI

      private:

        const GFS& _gfs;
        const RankIndex _rank;
        RankVector _ranks; // vector to identify unique decomposition
        GhostVector _ghosts; //vector to identify ghost dofs
        int _verbose; //verbosity

        //! The actual communication interface used when algorithm requires InteriorBorder_All_Interface.
        InterfaceType _interiorBorder_all_interface;

        //! The actual communication interface used when algorithm requires All_All_Interface.
        InterfaceType _all_all_interface;
      };

#if HAVE_MPI

      template<typename GFS>
      template<typename M, typename C>
      void ParallelHelper<GFS>::createIndexSetAndProjectForAMG(M& m, C& c)
      {

        using Backend::native;

        const bool is_bcrs_matrix =
          std::is_same<
            typename ISTL::tags::container<
              Backend::Native<M>
              >::type::base_tag,
          ISTL::tags::bcrs_matrix
          >::value;

        const bool block_type_is_field_matrix =
          std::is_same<
            typename ISTL::tags::container<
              typename Backend::Native<M>::block_type
              >::type::base_tag,
          ISTL::tags::field_matrix
          >::value;

        // We assume M to be a BCRSMatrix in the following, so better check for that
        static_assert(is_bcrs_matrix && block_type_is_field_matrix, "matrix structure not compatible with AMG");

        // ********************************************************************************
        // In the following, the code will always assume that all DOFs stored in a single
        // block of the BCRSMatrix are attached to the same entity and can be handled
        // identically. For that reason, the code often restricts itself to inspecting the
        // first entry of the blocks in the diverse BlockVectors.
        // ********************************************************************************

        typedef typename GFS::Traits::GridViewType GV;
        typedef typename RankVector::size_type size_type;
        const GV& gv = _gfs.gridView();

        // Do we need to communicate at all?
        const bool need_communication = _gfs.gridView().comm().size() > 1;

        // First find out which dofs we share with other processors
        using BoolVector = Backend::Vector<GFS,bool>;
        BoolVector sharedDOF(_gfs, false);

        if (need_communication)
          {
            SharedDOFDataHandle<GFS,BoolVector> data_handle(_gfs,sharedDOF,false);
            _gfs.gridView().communicate(data_handle,_all_all_interface,Dune::ForwardCommunication);
          }

        // Count shared dofs that we own
        typedef typename C::ParallelIndexSet::GlobalIndex GlobalIndex;
        GlobalIndex count = 0;

        for (size_type i = 0; i < sharedDOF.N(); ++i)
          if (owned_for_amg(i) && native(sharedDOF)[i][0])
            ++count;

        dverb << gv.comm().rank() << ": shared block count is " << count.touint() << std::endl;

        // Communicate per-rank count of owned and shared DOFs to all processes.
        std::vector<GlobalIndex> counts(_gfs.gridView().comm().size());
        _gfs.gridView().comm().allgather(&count, 1, &(counts[0]));

        // Compute start index start_p = \sum_{i=0}^{i<p} counts_i
        GlobalIndex start = std::accumulate(counts.begin(),counts.begin() + _rank,GlobalIndex(0));

        using GIVector = Dune::PDELab::Backend::Vector<GFS,GlobalIndex>;
        GIVector scalarIndices(_gfs, std::numeric_limits<GlobalIndex>::max());

        for (size_type i = 0; i < sharedDOF.N(); ++i)
          if (owned_for_amg(i) && native(sharedDOF)[i][0])
            {
              native(scalarIndices)[i][0] = start;
              ++start;
            }

        // Publish global indices for the shared DOFS to other processors.
        if (need_communication)
          {
            MinDataHandle<GFS,GIVector> data_handle(_gfs,scalarIndices);
            _gfs.gridView().communicate(data_handle,_interiorBorder_all_interface,Dune::ForwardCommunication);
          }

        // Setup the index set
        c.indexSet().beginResize();
        for (size_type i=0; i<scalarIndices.N(); ++i)
          {
            Dune::OwnerOverlapCopyAttributeSet::AttributeSet attr;
            if(native(scalarIndices)[i][0] != std::numeric_limits<GlobalIndex>::max())
              {
                // global index exist in index set
                if (owned_for_amg(i))
                  {
                    // This dof is managed by us.
                    attr = Dune::OwnerOverlapCopyAttributeSet::owner;
                  }
                else
                  {
                    attr = Dune::OwnerOverlapCopyAttributeSet::copy;
                  }
                c.indexSet().add(native(scalarIndices)[i][0], typename C::ParallelIndexSet::LocalIndex(i,attr));
              }
          }
        c.indexSet().endResize();

        // Compute neighbors using communication
        std::set<int> neighbors;

        if (need_communication)
          {
            GFSNeighborDataHandle<GFS,int> data_handle(_gfs,_rank,neighbors);
            _gfs.gridView().communicate(data_handle,_all_all_interface,Dune::ForwardCommunication);
          }

        c.remoteIndices().setNeighbours(neighbors);
        c.remoteIndices().template rebuild<false>();
      }

#endif // HAVE_MPI

      template<int s, bool isFakeMPIHelper>
      struct CommSelector
      {
        typedef Dune::Amg::SequentialInformation type;
      };


#if HAVE_MPI

      // Need MPI for OwnerOverlapCopyCommunication
      template<int s>
      struct CommSelector<s,false>
      {
        typedef OwnerOverlapCopyCommunication<bigunsignedint<s>,int> type;
      };

#endif // HAVE_MPI

      template<typename T>
      void assertParallelUG(T comm)
      {}

#if HAVE_UG && PDELAB_SEQUENTIAL_UG
      template<int dim>
      void assertParallelUG(Dune::CollectiveCommunication<Dune::UGGrid<dim> > comm)
      {
        static_assert(Dune::AlwaysFalse<Dune::UGGrid<dim> >::value, "Using sequential UG in parallel environment");
      };
#endif
      //! \} group Backend

    } // namespace ISTL
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_PARALLELHELPER_HH
