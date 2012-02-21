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
#include <dune/pdelab/common/typetree/transformation.hh>
#include <dune/pdelab/common/typetree/visitor.hh>
#include <dune/pdelab/constraints/constraintstransformation.hh>

#include <dune/pdelab/gridfunctionspace/lexicographicordering.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

#ifndef DOXYGEN // don't use an anomyous namespace - it breaks friend declarations

    //! Visitor for updating the complete GFS tree
    struct UpdateVisitor
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename LeafNode, typename TreePath>
      void leaf(LeafNode& node, TreePath treePath) const
      {
        node.shallowUpdate();
      }

      template<typename Node, typename TreePath>
      void post(Node& node, TreePath treePath) const
      {
        node.shallowUpdate();
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

      typedef G GridView;

      //! \brief vector backend
      typedef B BackendType;

      typedef B Backend;

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

      typedef Mapper OrderingTag;

      //! extract type for storing constraints
      template<typename E>
      struct ConstraintsContainer
      {
        //! \brief define Type as the Type of a container of E's
        typedef ConstraintsTransformation<typename GridFunctionSpace::Ordering::Traits::DOFIndex,E> Type;
      private:
        ConstraintsContainer ();
      };

      //! assumes all children are up to date.
      void shallowUpdate() { gfs().ordering().update(); }

      //! recalculate sizes
      void update ()
      {
        TypeTree::applyToTree(gfs(),UpdateVisitor());
      }



      //! get dimension of root finite element space
      typename Traits::SizeType globalSize () const
      {
        return gfs().ordering()->size();
      }

      //! get dimension of this finite element space
      typename Traits::SizeType size () const
      {
        return gfs().ordering()->size();
      }

      //! get max dimension of shape function space
      typename Traits::SizeType maxLocalSize () const
      {
        // this is bullshit !
        return gfs().ordering()->maxLocalSize();
      }

      //! get grid view
      const typename Traits::GridViewType& gridview () const DUNE_DEPRECATED
      {
        return gfs().template child<0>().gridView();
      }

      //! get grid view
      const typename Traits::GridViewType& gridView () const
      {
        return gfs().template child<0>().gridView();
      }

      //! map index from our index set [0,size()-1] to root index set
      typename Traits::SizeType upMap (typename Traits::SizeType i) const
      { return i; }

      //! map index from child i's index set into our index set
      template<int i>
      typename Traits::SizeType subMap(typename Traits::SizeType j) const
      { return subMap(i,j); }

      typename Traits::SizeType subMap(typename Traits::SizeType i,
                                       typename Traits::SizeType j) const
      { return gfs().ordering().subMap(i, j); }

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

      B& backend()
      {
        return _backend;
      }

      const B& backend() const
      {
        return _backend;
      }

      PowerCompositeGridFunctionSpaceBase(const B& backend)
        : _backend(backend)
      {}

      const std::string& name() const
      {
        return _name;
      }

      void name(const std::string& name)
      {
        _name = name;
      }

    private:

      B _backend;
      std::string _name;

    };

  }

}
//! \}
#endif // DUNE_PDELAB_POWERCOMPOSITEGRIDFUNCTIONSPACEBASE_HH
