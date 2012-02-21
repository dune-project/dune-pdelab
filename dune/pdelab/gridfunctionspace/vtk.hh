#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_VTK_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_VTK_HH

#include <vector>
#include <sstream>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/gridfunctionspace/lfscontainerindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/pdelab/common/typetree/visitor.hh>
#include <dune/pdelab/common/typetree/traversal.hh>

namespace Dune {
  namespace PDELab {

    template<typename LFS, typename Data>
    class DGFTreeLeafFunction;

    template<typename LFS, typename Data>
    class DGFTreeVectorFunction;

    template<typename VTKWriter, typename GFS, typename X>
    struct vtk_output_collector;


    //! Helper class for common data of a DGFTree.
    template<typename GFS, typename X>
    class DGFTreeCommonData
    {

      template<typename LFS, typename Data>
      friend class DGFTreeLeafFunction;

      template<typename LFS, typename Data>
      friend class DGFTreeVectorFunction;

      template<typename, typename, typename>
      friend struct vtk_output_collector;

      typedef LocalFunctionSpace<GFS> LFS;
      typedef LFSContainerIndexCache<LFS> LFSCache;
      typedef typename X::template ConstLocalView<LFSCache> XView;
      typedef LocalVector<typename X::ElementType> XLocalVector;
      typedef typename GFS::Traits::GridView::template Codim<0>::Entity Cell;
      typedef typename GFS::Traits::SizeType size_type;
      typedef typename GFS::Traits::GridView::IndexSet IndexSet;

    public:

      DGFTreeCommonData(const GFS& gfs, const X& x)
        : _lfs(gfs)
        , _lfs_cache(_lfs)
        , _x_view(x)
        , _x_local(_lfs.maxSize())
        , _current_cell_index(std::numeric_limits<size_type>::max())
        , _index_set(gfs.gridView().indexSet())
      {}

    public:

      void bind(const Cell& cell)
      {
        const size_type cell_index = _index_set.index(cell);
        if (_current_cell_index == cell_index)
          return;

        _lfs.bind(cell);
        _lfs_cache.update();
        _x_view.bind(_lfs_cache);
        _x_view.read(_x_local);
        _x_view.unbind();
        _current_cell_index = cell_index;
      }

      LFS _lfs;
      LFSCache _lfs_cache;
      XView _x_view;
      XLocalVector _x_local;
      size_type _current_cell_index;
      const IndexSet& _index_set;

    };



    template<typename LFS, typename Data>
    class DGFTreeLeafFunction
      : public TypeTree::LeafNode
      , public GridFunctionInterface<GridFunctionTraits<
                                       typename LFS::Traits::GridView,
                                       typename BasisInterfaceSwitch<
                                         typename FiniteElementInterfaceSwitch<
                                           typename LFS::Traits::FiniteElement
                                           >::Basis
                                         >::RangeField,
                                           BasisInterfaceSwitch<
                                           typename FiniteElementInterfaceSwitch<
                                             typename LFS::Traits::FiniteElement
                                             >::Basis
                                           >::dimRange,
                                             typename BasisInterfaceSwitch<
                                             typename FiniteElementInterfaceSwitch<
                                               typename LFS::Traits::FiniteElement
                                               >::Basis
                                             >::Range
                                       >,
                                     DGFTreeLeafFunction<LFS,Data>
                                     >
    {

      typedef BasisInterfaceSwitch<
        typename FiniteElementInterfaceSwitch<
          typename LFS::Traits::FiniteElement
          >::Basis
        > BasisSwitch;

      typedef GridFunctionInterface<
        GridFunctionTraits<
          typename LFS::Traits::GridView,
          typename BasisSwitch::RangeField,
          BasisSwitch::dimRange,
          typename BasisSwitch::Range
          >,
        DGFTreeLeafFunction<LFS,Data>
        > BaseT;

    public:
      typedef typename BaseT::Traits Traits;

      DGFTreeLeafFunction (const LFS& lfs, const shared_ptr<Data>& data)
        : _lfs(lfs)
        , _data(data)
        , _basis(lfs.maxSize())
      {}

      // Evaluate
      void evaluate (const typename Traits::ElementType& e,
                     const typename Traits::DomainType& x,
                     typename Traits::RangeType& y) const
      {
        _data->bind(e);

        typedef FiniteElementInterfaceSwitch<
          typename LFS::Traits::FiniteElement
          > FESwitch;

        y = 0;

        FESwitch::basis(_lfs.finiteElement()).evaluateFunction(x,_basis);
        for (std::size_t i = 0; i < _lfs.size(); ++i)
          y.axpy(_data->_x_local(_lfs,i),_basis[i]);
      }

      //! get a reference to the GridView
      const typename Traits::GridViewType& gridView() const
      {
        return _lfs.gridFunctionSpace().gridView();
      }

      const LFS& localFunctionSpace() const
      {
        return _lfs;
      }

    private:

      const LFS& _lfs;
      const shared_ptr<Data> _data;
      mutable std::vector<typename Traits::RangeType> _basis;

    };



    template<typename LFS, typename Data>
    class DGFTreeVectorFunction
      : public TypeTree::LeafNode
      , public GridFunctionInterface<GridFunctionTraits<
                                       typename LFS::Traits::GridView,
                                       typename BasisInterfaceSwitch<
                                         typename FiniteElementInterfaceSwitch<
                                           typename LFS::ChildType::Traits::FiniteElement
                                           >::Basis
                                         >::RangeField,
                                       LFS::CHILDREN,
                                       Dune::FieldVector<
                                         typename BasisInterfaceSwitch<
                                           typename FiniteElementInterfaceSwitch<
                                             typename LFS::ChildType::Traits::FiniteElement
                                             >::Basis
                                           >::RangeField,
                                         LFS::CHILDREN
                                         >
                                       >,
                                     DGFTreeVectorFunction<LFS,Data>
                                     >
    {

      typedef BasisInterfaceSwitch<
        typename FiniteElementInterfaceSwitch<
          typename LFS::ChildType::Traits::FiniteElement
          >::Basis
        > BasisSwitch;

      dune_static_assert(BasisSwitch::dimRange == 1,
                         "Automatic conversion to vector-valued function only supported for scalar components");

      typedef GridFunctionInterface<
        GridFunctionTraits<
          typename LFS::Traits::GridView,
          typename BasisSwitch::RangeField,
          LFS::CHILDREN,
          Dune::FieldVector<
            typename BasisSwitch::RangeField,
            LFS::CHILDREN
            >
          >,
        DGFTreeVectorFunction<LFS,Data>
        > BaseT;

    public:

      typedef typename BaseT::Traits Traits;
      typedef typename LFS::ChildType ChildLFS;
      typedef typename ChildLFS::Traits::FiniteElement::Traits::LocalBasisType::Traits::RangeFieldType RF;
      typedef typename ChildLFS::Traits::FiniteElement::Traits::LocalBasisType::Traits::RangeType RT;

      DGFTreeVectorFunction (const LFS& lfs, const shared_ptr<Data>& data)
        : _lfs(lfs)
        , _data(data)
        , _basis(lfs.maxSize())
      {}

      void evaluate (const typename Traits::ElementType& e,
                     const typename Traits::DomainType& x,
                     typename Traits::RangeType& y) const
      {
        _data->bind(e);

        typedef FiniteElementInterfaceSwitch<
          typename ChildLFS::Traits::FiniteElement
          > FESwitch;

        y = 0;

        for (std::size_t k = 0; k < LFS::CHILDREN; ++k)
          {
            const ChildLFS& child_lfs = _lfs.child(k);
            FESwitch::basis(child_lfs.finiteElement()).evaluateFunction(x,_basis);

            for (std::size_t i = 0; i < child_lfs.size(); ++i)
              y[k] += _data->_x_local(child_lfs,i) * _basis[i];
          }
      }

      //! get a reference to the GridView
      const typename Traits::GridViewType& gridView() const
      {
        return _lfs.gridFunctionSpace().gridView();
      }

      const LFS& localFunctionSpace() const
      {
        return _lfs;
      }

    private:

      const LFS& _lfs;
      const shared_ptr<Data> _data;
      mutable std::vector<typename BasisSwitch::Range> _basis;

    };




    template<typename VTKWriter, typename Data>
    struct add_solution_to_vtk_writer_visitor
      : public TypeTree::DefaultVisitor
      , public TypeTree::DynamicTraversal
    {


      template<typename LFS, typename Child, typename TreePath>
      struct VisitChild
      {

        //! Do not descend into children of VectorGridFunctionSpace
        static const bool value = !is_same<
          typename LFS::Traits::GridFunctionSpace::ImplementationTag,
          VectorGridFunctionSpaceTag
          >::value;

      };

      //! Helper function for extracting (or building) the component name and adding
      //! the component to the VTKWriter.
      template<typename DGF, typename TreePath>
      void add_to_vtk_writer(const shared_ptr<DGF>& dgf, TreePath tp)
      {
        std::string name(dgf->localFunctionSpace().gridFunctionSpace().name());
        if (name.empty())
          {
            std::stringstream name_stream;

            // Build a simple name based on the component's TreePath (e.g. 0_2_3)
            for (std::size_t i = 0; i < tp.size(); ++i)
              name_stream << (i > 0 ? "_" : "") << tp.element(i);
            name = name_stream.str();
          }
        vtk_writer.addCellData(new VTKGridFunctionAdapter<DGF>(dgf,name.c_str()));
      }

      //! Tag dispatch-based switch that creates a vector-valued function for a VectorGridFunctionSpace.
      /**
       * This version handles an actual VectorGridFunctionSpace.
       */
      template<typename LFS, typename TreePath>
      void add_vector_solution(const LFS& lfs, TreePath tp, VectorGridFunctionSpaceTag tag)
      {
        add_to_vtk_writer(make_shared<DGFTreeVectorFunction<LFS,Data> >(lfs,data),tp);
      }

      //! Tag dispatch-based switch that creates a vector-valued function for a VectorGridFunctionSpace.
      /**
       * This is the default version for different types of spaces that does nothing.
       */
      template<typename LFS, typename TreePath, typename Tag>
      void add_vector_solution(const LFS& lfs, TreePath tp, Tag tag)
      {
        // do nothing here - not a vector space
      }

      //! Handle VectorGridFunctionSpace components in here.
      template<typename LFS, typename TreePath>
      void post(const LFS& lfs, TreePath tp)
      {
        add_vector_solution(lfs,tp,typename LFS::Traits::GridFunctionSpace::ImplementationTag());
      }

      //! Create a standard leaf function for leaf GridFunctionSpaces.
      template<typename LFS, typename TreePath>
      void leaf(const LFS& lfs, TreePath tp)
      {
        add_to_vtk_writer(make_shared<DGFTreeLeafFunction<LFS,Data> >(lfs,data),tp);
      }


      add_solution_to_vtk_writer_visitor(VTKWriter& vtk_writer_, shared_ptr<Data> data_)
        : vtk_writer(vtk_writer_)
        , data(data_)
      {}

      VTKWriter& vtk_writer;
      shared_ptr<Data> data;

    };

    template<typename VTKWriter, typename GFS, typename X>
    struct vtk_output_collector
    {

      //! Common data container (hierarchic LFS, global solution data etc.)
      typedef DGFTreeCommonData<GFS,X> Data;

      vtk_output_collector& add_solution()
      {
        add_solution_to_vtk_writer_visitor<VTKWriter,Data> visitor(_vtk_writer,_data);
        TypeTree::applyToTree(_data->_lfs,visitor);
        return *this;
      }

      template<typename Factory, typename TreePath>
      vtk_output_collector& add_cell_function(Factory factory, TreePath tp, std::string name)
      {
        typedef typename TypeTree::extract_child_type<typename Data::LFS,TreePath>::type LFS;
        typedef typename Factory::template create_type<LFS,Data>::type DGF;
        _vtk_writer.addCellData(new VTKGridFunctionAdapter<DGF>(factory.create(TypeTree::extract_child(_data->_lfs,tp),_data),name));
        return *this;
      }

      vtk_output_collector(VTKWriter& vtk_writer, const GFS& gfs, const X& x)
        : _vtk_writer(vtk_writer)
        , _data(make_shared<Data>(gfs,x))
      {}

      VTKWriter& _vtk_writer;
      shared_ptr<Data> _data;

    };


    template<typename VTKWriter, typename GFS, typename X>
    vtk_output_collector<VTKWriter,GFS,X> add_solution_to_vtk_writer(VTKWriter& vtk_writer, const GFS& gfs, const X& x)
    {
      vtk_output_collector<VTKWriter,GFS,X> collector(vtk_writer,gfs,x);
      collector.add_solution();
      return std::move(collector);
    }


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_VTK_HH
