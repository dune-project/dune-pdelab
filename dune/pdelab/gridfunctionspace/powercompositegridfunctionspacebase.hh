// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_POWERCOMPOSITEGRIDFUNCTIONSPACEBASE_HH
#define DUNE_PDELAB_POWERCOMPOSITEGRIDFUNCTIONSPACEBASE_HH

#include <cstddef>
#include <functional>
#include <vector>

#include <dune/pdelab/common/typetree/fixedcapacitystack.hh>
#include <dune/pdelab/common/typetree/traversal.hh>
#include <dune/pdelab/common/typetree/traversalutilities.hh>
#include <dune/pdelab/common/typetree/utility.hh>
#include <dune/pdelab/common/typetree/visitor.hh>
#include <dune/pdelab/constraints/constraintstransformation.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

#ifndef DOXYGEN // don't use an anomyous namespace - it breaks friend declarations

    //! We put the actual method in a base class because we want to use it with different tree iteration patterns
    struct SetupVisitorBase
      : public TypeTree::DefaultVisitor // do not fix the traversal depth yet
      , public TypeTree::DynamicTraversal
    {

      template<typename CompositeGFS, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(CompositeGFS& cgfs, const Child& child, TreePath treePath, ChildIndex childIndex) const
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
      void leaf(const LeafNode& node, TreePath treePath) const
      {
        node.update();
      }

      template<typename LeafNode, typename TreePath>
      void post(const LeafNode& node, TreePath treePath) const
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
      std::size_t operator()(const Node& node, TreePath treePath) const
      {
        return node.dataHandleSize(e);
      }

      DataHandleSize(const Entity& entity)
        : e(entity)
      {}

      const Entity& e;

    };

    //! Visitor for retrieving the global DOF indices of a given entity.
    template<typename Entity, typename Container, std::size_t gfs_depth>
    struct DataHandleGlobalIndicesVisitor
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

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
      }

      const Entity& e;
      Container& g;
      std::size_t pos;
      TypeTree::FixedCapacityStack<std::size_t,gfs_depth> offsets;

    };

#endif // DOXYGEN

    //! Trait class for the multi component grid function spaces
    template<typename G, typename B, typename M, std::size_t k>
    struct PowerCompositeGridFunctionSpaceTraits
    {
      enum{
        //! \brief True if this grid function space is composed of others.
        isComposite = 1,
        //! \brief number of child spaces
        noChilds = k
      };

      const static std::size_t CHILDREN = k;

      //! \brief the grid view where grid function is defined upon
      typedef G GridViewType;

      //! \brief vector backend
      typedef B BackendType;

      //! \brief mapper
      typedef M MapperType;

      //! \brief short cut for size type exported by Backend
      typedef typename B::size_type SizeType;
    };

    //! Mixin class providing common functionality of PowerGridFunctionSpace and CompositeGridFunctionSpace
    template<typename GridFunctionSpace, typename GV, typename B, typename Mapper, std::size_t k>
    class PowerCompositeGridFunctionSpaceBase
    {

#ifndef DOXYGEN

      friend struct SetupVisitorBase;

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
      typedef PowerCompositeGridFunctionSpaceTraits<GV,B,Mapper,k> Traits;

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
        return offset[Traits::CHILDREN];
      }

      //! get dimension of this finite element space
      typename Traits::SizeType size () const
      {
        return offset[Traits::CHILDREN];
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
        DataHandleGlobalIndicesVisitor<EntityType,std::vector<SizeType>,TypeTree::TreeInfo<GridFunctionSpace>::depth> visitor(e,global);
        TypeTree::applyToTree(gfs(),visitor);
      }

    protected:

      typename Traits::SizeType childGlobalSize[Traits::CHILDREN];
      typename Traits::SizeType childLocalSize[Traits::CHILDREN];
      typename Traits::SizeType offset[Traits::CHILDREN+1];
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
#endif // DUNE_PDELAB_POWERCOMPOSITEGRIDFUNCTIONSPACEBASE_HH
