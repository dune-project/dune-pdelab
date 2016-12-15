// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_SUBSPACELOCALFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_SUBSPACELOCALFUNCTIONSPACE_HH

/** \file
 *  \brief Support infrastructure to make LocalFunctionSpaces of GridFunctionSubSpace work.
 */

// nothing in here is of interest to our users
#ifndef DOXYGEN

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

namespace Dune {
  namespace PDELab {

    namespace gfs {

      // ********************************************************************************
      // utility functions for copying a tree path to an iterator range in reverse order
      // ********************************************************************************

      template<typename It>
      void extract_tree_path_elements(TypeTree::TreePath<>, It it)
      {
        // end of recursion
      }

      template<typename TP, typename It>
      void extract_tree_path_elements(TP, It it)
      {
        *it = TypeTree::TreePathBack<TP>::value;
        extract_tree_path_elements(typename TypeTree::TreePathPopBack<TP>::type(),++it);
      }

      // ********************************************************************************
      // intermediate base class for LocalFunctionSpaces of GridFunctionSubSpace
      // ********************************************************************************

      // This class works for subspaces of any kind of underlying space (leaf, Power, Composite)
      // thanks to the magic of perfect forwarding
      template<typename GFS, typename LFS>
      class SubSpaceLocalFunctionSpaceNode
        : public LFS
      {

      public:

        typedef typename LFS::Traits Traits;

        template<typename... T>
        SubSpaceLocalFunctionSpaceNode(T&&... t)
          : LFS(std::forward<T>(t)...)
        {
          extract_tree_path_elements(typename GFS::SubSpacePath(),_tree_path.begin());
        }

        // modify bind to fill up the DOFIndices with our subspace path
        template<typename E>
        void bind(const E& e)
        {
          LFS::bind(e);
          for (auto di= this->_dof_index_storage.begin(), end=this->_dof_index_storage.end();
               di!=end; ++di)
            complete_dof_index(*di);
        }

        std::size_t subSpaceDepth() const
        {
          return this->gridFunctionSpace().subSpaceDepth();
        }

      private:

        //! Extends DOFIndex up to the root of the original DOFIndex tree.
        void complete_dof_index(typename Traits::DOFIndex& di) const
        {
          std::copy(_tree_path.begin(),_tree_path.end(),std::back_inserter(di.treeIndex()));
        }

        std::array<std::size_t,TypeTree::TreePathSize<typename GFS::SubSpacePath>::value> _tree_path;

      };

      // forward declaration for use in LocalFunctionSpace specialization.
      template<typename GFS, typename TreePath>
      class GridFunctionSubSpace;

    } // namespace gfs

    // specialized version of LocalFunctionSpace interface class
    // This class injects SubSpaceLocalFunctionSpaceNode as an intermediate base class
    // to fix the DOFIndex handling in bind().
    template <typename BaseGFS, typename SubSpaceTreePath>
    class LocalFunctionSpace<gfs::GridFunctionSubSpace<BaseGFS,SubSpaceTreePath>, AnySpaceTag>
      : public gfs::SubSpaceLocalFunctionSpaceNode<gfs::GridFunctionSubSpace<BaseGFS,SubSpaceTreePath>,
                                                   typename TypeTree::TransformTree<
                                                     gfs::GridFunctionSubSpace<BaseGFS,SubSpaceTreePath>,
                                                     gfs_to_lfs<gfs::GridFunctionSubSpace<
                                                                  BaseGFS,
                                                                  SubSpaceTreePath
                                                                  >
                                                                >
                                                     >::Type
                                                   >
    {

      typedef gfs::GridFunctionSubSpace<BaseGFS,SubSpaceTreePath> GFS;

      typedef gfs::SubSpaceLocalFunctionSpaceNode<
        GFS,
        typename TypeTree::TransformTree<
          GFS,
          gfs_to_lfs<GFS>
          >::Type
        > BaseT;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ClearSizeVisitor;

      template<typename,bool>
      friend struct ComputeSizeVisitor;

      template<typename,bool>
      friend struct FillIndicesVisitor;

    public:

      LocalFunctionSpace(const GFS & gfs)
        : BaseT(TypeTree::TransformTree<GFS,gfs_to_lfs<GFS> >::transform(gfs))
      {
        this->_dof_indices = &(this->_dof_index_storage);
        this->setup(*this);
      }

      LocalFunctionSpace(std::shared_ptr<const GFS> pgfs)
        : BaseT(TypeTree::TransformTree<GFS,gfs_to_lfs<GFS> >::transform_storage(pgfs))
      {
        this->_dof_indices = &(this->_dof_index_storage);
        this->setup(*this);
      }

      LocalFunctionSpace(const LocalFunctionSpace & lfs)
        : BaseT(lfs)
      {
        // We need to reset the DOFIndex storage pointers in the new LFS tree,
        // as they are still pointing to the _dof_index_storage of the
        // old tree.
        this->_dof_indices = &(this->_dof_index_storage);
        this->setup(*this);
      }

    };

  } // namespace PDELab
} // namespace Dune

#endif // DOXYGEN

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_SUBSPACELOCALFUNCTIONSPACE_HH
