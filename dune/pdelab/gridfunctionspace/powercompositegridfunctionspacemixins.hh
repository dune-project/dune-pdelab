// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_POWERCOMPOSITEGRIDFUNCTIONSPACEMIXINS_HH
#define DUNE_PDELAB_POWERCOMPOSITEGRIDFUNCTIONSPACEMIXINS_HH

#include "../common/typetree.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{


    //! Mixin class providing common functionality of PowerGridFunctionSpace and CompositeGridFunctionSpace
    template<typename GridFunctionSpace, typename TreeNode, typename Mapper>
    class PowerCompositeGridFunctionSpaceBase
    {

#ifndef DOXYGEN

      //! We put the actual method in a base class because we want to use it with different tree iteration patterns
      struct SetupVisitorBase
        : public TypeTree::DefaultVisitor
      {

        static const TypeTree::TreePathType::Type treePathType = TypeTree::TreePathType::dynamic;

        template<typename CompositeGFS, typename Child, typename TreePath, typename ChildIndex>
        void afterChild(CompositeGFS& cgfs, const Child& child, TreePath treePath, ChildIndex childIndex)
        {
          cgfs.childGlobalSize[childIndex] = child.globalSize();
          cgfs.childLocalSize[childIndex] = child.maxLocalSize();
        }

      };

      //! Visitor for setting up the GFS from pre-initialized children
      struct SetupVisitor
        : public SetupVisitorBase
        , public TypeTree::VisitDirectChildren
      {};

      //! Visitor for updating the complete GFS tree
      struct UpdateVisitor
        : public SetupVisitorBase
        , public TypeTree::VisitTree
      {

        template<typename LeafNode, typename TreePath>
        void leaf(const LeafNode& node, TreePath treePath)
        {
          node.update();
        }

        template<typename LeafNode, typename TreePath>
        void post(const LeafNode& node, TreePath treePath)
        {
          node.calculateSizes();
        }

      };

      //! Functor for dataHandleFixedSize()
      struct DataHandleFixedSize
      {

        template<typename Node, typename TreePath>
        bool operator()(const Node& node, TreePath treePath) const
        {
          return node.dataHandleFixedSize(dim,codim);
        }

        DataHandleFixedSize(int dimension, int codimension)
          : dim(dimension)
          , codim(codimension)
        {}

        const int dim;
        const int codim;

      };

      //! Functor for dataHandleContains()
      struct DataHandleContains
      {

        template<typename Node, typename TreePath>
        bool operator()(const Node& node, TreePath treePath) const
        {
          return node.dataHandleContains(dim,codim);
        }

        DataHandleContains(int dimension, int codimension)
          : dim(dimension)
          , codim(codimension)
        {}

        const int dim;
        const int codim;

      };

      //! Functor for dataHandleSize()
      template<typename Entity>
      struct DataHandleSize
      {

        template<typename Node, typename TreePath>
        bool operator()(const Node& node, TreePath treePath) const
        {
          return node.dataHandleSize(e);
        }

        DataHandleSize(const Entity& entity)
          : e(entity)
        {}

        const Entity& e;

      };

      //! Visitor for retrieving the global DOF indices of a given entity.
      template<typename Entity, typename Container>
      struct DataHandleGlobalIndicesVisitor
        : public TypeTree::TreeVisitor
      {

        static const TypeTree::TreePathType::Type treePathType = TypeTree::TreePathType::dynamic;

        template<typename Node, typename TreePath>
        void leaf(const Node& node, TreePath treePath)
        {
          pos += node.dataHandleGlobalIndices(e,g,pos,false);
        }

        template<typename CompositeGFS, typename Child, typename TreePath, typename ChildIndex>
        void beforeChild(CompositeGFS& cgfs, const Child& child, TreePath treePath, ChildIndex childIndex)
        {
          offsets.push_back(pos);
        }

        template<typename CompositeGFS, typename Child, typename TreePath, typename ChildIndex>
        void afterChild(CompositeGFS& cgfs, const Child& child, TreePath treePath, ChildIndex childIndex)
        {
          std::size_t offset = offsets.back();
          offsets.pop_back();
          for (std::size_t i = offset; i < pos; ++i)
            g[i] = cgfs.subMap(childIndex,g[i]);
        }

        DataHandleGlobalIndicesVisitor(const Entity& entity, Container& global)
          : e(entity)
          , g(global)
          , pos(0)
        {
          // a reasonable upper bound for the tree depth - this way, we avoid reallocations
          offsets.reserve(16);
        }

        const Entity& e;
        Container& g;
        std::size_t pos;
        std::vector<std::size_t> offsets;

      };

      const GridFunctionSpace& gfs() const
      {
        return static_cast<const GridFunctionSpace&>(*this);
      }

      GridFunctionSpace& gfs()
      {
        return static_cast<GridFunctionSpace&>(*this);
      }

#endif // DOXYGEN

    public:

      //! export traits class
      typedef PowerCompositeGridFunctionSpaceTraits<typename TreeNode::template Child<0>::Type::Traits::GridViewType,
                                                    typename TreeNode::template Child<0>::Type::Traits::BackendType,
                                                    Mapper,
                                                    TreeNode::CHILDREN> Traits;

      //! extract type of container storing Es
      template<typename E>
      struct VectorContainer
      {
        //! \brief define Type as the Type of a container of E's
        typedef typename Traits::BackendType::template VectorContainer<GridFunctionSpace,E> Type;
      private:
        VectorContainer ();
      };

      //! extract type for storing constraints
      template<typename E>
      struct ConstraintsContainer
      {
        //! \brief define Type as the Type of a container of E's
        typedef ConstraintsTransformation<typename Traits::SizeType,E> Type;
      private:
        ConstraintsContainer ();
      };


      //! recalculate sizes
      void update ()
      {
        TypeTree::applyToTree(gfs(),UpdateVisitor());
      }

      //! get dimension of root finite element space
      typename Traits::SizeType globalSize () const
      {
        return offset[TreeNode::CHILDREN];
      }

      //! get dimension of this finite element space
      typename Traits::SizeType size () const
      {
        return offset[TreeNode::CHILDREN];
      }

      //! get max dimension of shape function space
      typename Traits::SizeType maxLocalSize () const
      {
        // this is bullshit !
        return maxlocalsize;
      }

      //! get grid view
      const typename Traits::GridViewType& gridview () const
      {
        return gfs().template child<0>().gridview();
      }

      //------------------------------
      // generic data handle interface
      //------------------------------

      //! returns true if data for this codim should be communicated
      bool dataHandleContains (int dim, int codim) const
      {
        return TypeTree::reduceOverLeafs(gfs(),DataHandleContains(dim,codim),std::logical_or<bool>(),false);
      }

      //! returns true if size per entity of given dim and codim is a constant
      bool dataHandleFixedSize (int dim, int codim) const
      {
        return TypeTree::reduceOverLeafs(gfs(),DataHandleFixedSize(dim,codim),std::logical_and<bool>(),true);
      }

      /*! how many objects of type DataType have to be sent for a given entity

        Note: Only the sender side needs to know this size.
      */
      template<typename EntityType>
      size_t dataHandleSize (const EntityType& e) const
      {
        return TypeTree::reduceOverLeafs(gfs(),DataHandleSize<EntityType>(e),std::plus<size_t>(),size_t(0));
      }

      //! return vector of global indices associated with the given entity
      template<typename EntityType, typename SizeType>
      void dataHandleGlobalIndices (const EntityType& e,
                                    std::vector<SizeType>& global) const
      {
        global.resize(dataHandleSize(e));
        DataHandleGlobalIndicesVisitor<EntityType,std::vector<SizeType> > visitor(e,global);
        TypeTree::applyToTree(gfs(),visitor);
      }

    protected:

      typename Traits::SizeType childGlobalSize[TreeNode::CHILDREN];
      typename Traits::SizeType childLocalSize[TreeNode::CHILDREN];
      typename Traits::SizeType offset[TreeNode::CHILDREN+1];
      typename Traits::SizeType maxlocalsize;

       void setup ()
      {
        TypeTree::applyToTree(gfs(),SetupVisitor());
        gfs().calculateSizes();
      }

    };

  }

}
//! \}
#endif // DUNE_PDELAB_POWERCOMPOSITEGRIDFUNCTIONSPACEMIXINS_HH
