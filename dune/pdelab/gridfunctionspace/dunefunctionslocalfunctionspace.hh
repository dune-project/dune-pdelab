// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_DUNEFUNCTIONSLOCALFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_DUNEFUNCTIONSLOCALFUNCTIONSPACE_HH

#include<vector>

#include <dune/common/stdstreams.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/typetree/typetree.hh>

#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/gridfunctionspace/localvector.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    namespace Experimental {

      template<typename LFS>
      struct LeafLFSMixin
        : public TypeTree::LeafNode
      {

        const auto& finiteElement() const
        {
          return static_cast<const LFS*>(this)->tree().finiteElement();
        }

        template<typename Tree>
        struct Traits
        {
          using FiniteElement = typename Tree::FiniteElement;
          using FiniteElementType = FiniteElement;
        };
      };

      template<typename GFS, typename TreePath = TypeTree::HybridTreePath<>>
      class LocalFunctionSpace
        : public LeafLFSMixin<LocalFunctionSpace<GFS,TreePath>>
      {

      public:

        using Basis = typename GFS::Basis;
        using LocalView = typename Basis::LocalView;
        using Tree = TypeTree::ChildForTreePath<typename LocalView::Tree,TreePath>;
        using LocalIndexSet = typename Basis::LocalIndexSet;
        using DOFIndex = typename Basis::MultiIndex;

        template<typename LFS, typename C, typename Tag, bool fast>
        friend class LFSIndexCacheBase;

        struct Traits
          : public LeafLFSMixin<LocalFunctionSpace<GFS,TreePath>>::template Traits<Tree>
        {

          using GridFunctionSpace = GFS;
          using GridView          = typename GFS::Traits::GridView;
          using SizeType          = std::size_t;
          using DOFIndex          = typename Basis::MultiIndex;
          using ConstraintsType   = typename GFS::Traits::ConstraintsType;

        };

        using size_type = std::size_t;

        LocalFunctionSpace(std::shared_ptr<const GFS> gfs, TreePath tree_path = TreePath(), size_type offset = 0)
          : _gfs(gfs)
          , _local_view(gfs->basis())
          , _tree_path(tree_path)
          , _tree(TypeTree::child(_local_view.tree(),tree_path))
          , _local_index_set(gfs->basis().localIndexSet())
        {}

        size_type subSpaceDepth() const
        {
          return 0;
        }

        //! \brief get current size
        size_type size () const
        {
          return _local_view.size();
        }

        size_type maxSize () const
        {
          // _dof_indices is always as large as the max local size of the root GFS
          return _local_view.maxSize();
        }

        //! \brief map index in this local function space to root local function space
        size_type localIndex (size_type index) const
        {
          return _tree.localIndex(index);
        }

        DOFIndex dofIndex(size_type index) const
        {
          return _local_index_set.index(_tree.localIndex(index));
        }

        //! Returns the GridFunctionSpace underlying this LocalFunctionSpace.
        const GFS& gridFunctionSpace() const
        {
          return *_gfs;
        }

        void bind(const typename GFS::Traits::EntitySet::template Codim<0>::Entity& e)
        {
          _local_view.bind(e);
          _local_index_set.bind(_local_view);
        }

        const typename Traits::ConstraintsType& constraints() const
        {
          return _gfs->constraints();
        }

        const Tree& tree() const
        {
          return _tree;
        }

      private:

        std::shared_ptr<const GFS> _gfs;
        LocalView _local_view;
        TreePath _tree_path;
        const Tree& _tree;
        LocalIndexSet _local_index_set;

      };

      // forward declare GridFunctionSpace
      template<typename DFBasis, typename V, typename CE=NoConstraints>
      class GridFunctionSpace;


    } // namespace Experimental


    template<typename DFBasis, typename V, typename CE, typename TAG>
    class LocalFunctionSpace<Experimental::GridFunctionSpace<DFBasis,V,CE>,TAG>
      : public Experimental::LocalFunctionSpace<Experimental::GridFunctionSpace<DFBasis,V,CE>>
    {

      using GFS = Experimental::GridFunctionSpace<DFBasis,V,CE>;

    public:

      LocalFunctionSpace(std::shared_ptr<const GFS> gfs)
        : Experimental::LocalFunctionSpace<GFS>(gfs)
      {}

      LocalFunctionSpace(const GFS& gfs)
        : Experimental::LocalFunctionSpace<GFS>(stackobject_to_shared_ptr(gfs))
      {}

    };

    template<typename DFBasis, typename V, typename CE>
    class LocalFunctionSpace<Experimental::GridFunctionSpace<DFBasis,V,CE>,AnySpaceTag>
      : public Experimental::LocalFunctionSpace<Experimental::GridFunctionSpace<DFBasis,V,CE>>
    {

      using GFS = Experimental::GridFunctionSpace<DFBasis,V,CE>;

    public:

      LocalFunctionSpace(std::shared_ptr<const GFS> gfs)
        : Experimental::LocalFunctionSpace<GFS>(gfs)
      {}

      LocalFunctionSpace(const GFS& gfs)
        : Experimental::LocalFunctionSpace<GFS>(stackobject_to_shared_ptr(gfs))
      {}

    };

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_DUNEFUNCTIONSLOCALFUNCTIONSPACE_HH
