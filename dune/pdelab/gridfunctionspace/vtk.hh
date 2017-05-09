#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_VTK_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_VTK_HH

#include <vector>
#include <sstream>

#include <dune/common/exceptions.hh>

#include <dune/geometry/typeindex.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/typetree/visitor.hh>
#include <dune/typetree/traversal.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localvector.hh>

namespace Dune {

  template<typename GV>
  class VTKWriter;

  template<typename GV>
  class SubsamplingVTKWriter;

  template<typename GV>
  class VTKSequenceWriter;

  namespace PDELab {

    namespace vtk {

      namespace {

        template<typename VTKWriter>
        struct vtk_writer_traits;

        template<typename GV>
        struct vtk_writer_traits<Dune::VTKWriter<GV> >
        {
          typedef GV GridView;
        };

        template<typename GV>
        struct vtk_writer_traits<Dune::SubsamplingVTKWriter<GV> >
        {
          typedef GV GridView;
        };

        template<typename GV>
        struct vtk_writer_traits<Dune::VTKSequenceWriter<GV> >
        {
          typedef GV GridView;
        };

      }

      template<typename LFS, typename Data>
      class DGFTreeLeafFunction;

      template<typename LFS, typename Data>
      class DGFTreeVectorFunction;

      template<typename VTKWriter, typename Data>
      struct OutputCollector;


      //! Helper class for common data of a DGFTree.
      template<typename GFS, typename X, typename Pred>
      class DGFTreeCommonData
      {

        template<typename LFS, typename Data>
        friend class DGFTreeLeafFunction;

        template<typename LFS, typename Data>
        friend class DGFTreeVectorFunction;

        template<typename, typename>
        friend struct OutputCollector;

        typedef LocalFunctionSpace<GFS> LFS;
        typedef LFSIndexCache<LFS> LFSCache;
        typedef typename X::template ConstLocalView<LFSCache> XView;
        typedef LocalVector<typename X::ElementType> XLocalVector;
        using EntitySet = typename GFS::Traits::EntitySet;
        using Cell = typename EntitySet::Traits::Element;
        using IndexSet = typename EntitySet::Traits::IndexSet;
        typedef typename IndexSet::IndexType size_type;

        static const auto dim = EntitySet::dimension;

      public:

        typedef GFS GridFunctionSpace;
        typedef X Vector;
        typedef Pred Predicate;

        DGFTreeCommonData(const GFS& gfs, const X& x)
          : _lfs(gfs)
          , _lfs_cache(_lfs)
          , _x_view(x)
          , _x_local(_lfs.maxSize())
          , _index_set(gfs.entitySet().indexSet())
          , _current_cell_index(std::numeric_limits<size_type>::max())
        {}

      public:

        void bind(const Cell& cell)
        {
          auto cell_index = _index_set.uniqueIndex(cell);
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
        const IndexSet& _index_set;
        size_type _current_cell_index;

      };



      template<typename LFS, typename Data>
      class DGFTreeLeafFunction
        : public GridFunctionBase<GridFunctionTraits<
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

        typedef GridFunctionBase<
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
          : BaseT(lfs.gridFunctionSpace().dataSetType())
          , _lfs(lfs)
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
        : public GridFunctionBase<GridFunctionTraits<
                                         typename LFS::Traits::GridView,
                                         typename BasisInterfaceSwitch<
                                           typename FiniteElementInterfaceSwitch<
                                             typename LFS::ChildType::Traits::FiniteElement
                                             >::Basis
                                           >::RangeField,
                                         TypeTree::StaticDegree<LFS>::value,
                                         Dune::FieldVector<
                                           typename BasisInterfaceSwitch<
                                             typename FiniteElementInterfaceSwitch<
                                               typename LFS::ChildType::Traits::FiniteElement
                                               >::Basis
                                             >::RangeField,
                                           TypeTree::StaticDegree<LFS>::value
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

        static_assert(BasisSwitch::dimRange == 1,
                      "Automatic conversion to vector-valued function only supported for scalar components");

        typedef GridFunctionBase<
          GridFunctionTraits<
            typename LFS::Traits::GridView,
            typename BasisSwitch::RangeField,
            TypeTree::StaticDegree<LFS>::value,
            Dune::FieldVector<
              typename BasisSwitch::RangeField,
              TypeTree::StaticDegree<LFS>::value
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
          : BaseT(lfs.gridFunctionSpace().dataSetType())
          , _lfs(lfs)
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

          for (std::size_t k = 0; k < TypeTree::degree(_lfs); ++k)
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


      class DefaultFunctionNameGenerator
      {

      public:

        template<typename TreePath>
        std::string operator()(std::string component_name, TreePath tp) const
        {
          if (component_name.empty())
            {

              if (_prefix.empty() && _suffix.empty())
                {
                  DUNE_THROW(IOError,
                             "You need to either name all GridFunctionSpaces "
                             "written to the VTK file or provide a prefix / suffix.");
                }

              std::stringstream name_stream;

              if (!_prefix.empty())
                name_stream << _prefix << _separator;

              // Build a simple name based on the component's TreePath (e.g. 0_2_3)
              for (std::size_t i = 0; i < tp.size(); ++i)
                name_stream << (i > 0 ? _separator : "") << tp.element(i);

              if (!_suffix.empty())
                name_stream << _separator << _suffix;
              return name_stream.str();
            }
          else
            {
              // construct name from prefix, component name and suffix
              return _prefix + component_name + _suffix;
            }
        }

        DefaultFunctionNameGenerator& prefix(std::string prefix)
        {
          _prefix = prefix;
          return *this;
        }

        DefaultFunctionNameGenerator& suffix(std::string suffix)
        {
          _suffix = suffix;
          return *this;
        }

        DefaultFunctionNameGenerator& separator(std::string separator)
        {
          _separator = separator;
          return *this;
        }

        DefaultFunctionNameGenerator(std::string prefix = "",
                                        std::string suffix = "",
                                        std::string separator = "_")
          : _prefix(prefix)
          , _suffix(suffix)
          , _separator(separator)
        {}

      private:

        std::string _prefix;
        std::string _suffix;
        std::string _separator;

      };

      inline DefaultFunctionNameGenerator defaultNameScheme()
      {
        return DefaultFunctionNameGenerator();
      }


      template<typename VTKWriter, typename Data, typename NameGenerator>
      struct add_solution_to_vtk_writer_visitor
        : public TypeTree::DefaultVisitor
        , public TypeTree::DynamicTraversal
      {


        template<typename LFS, typename Child, typename TreePath>
        struct VisitChild
        {

          static const bool value =
            // Do not descend into children of VectorGridFunctionSpace
            !std::is_convertible<
              TypeTree::ImplementationTag<typename LFS::Traits::GridFunctionSpace>,
              VectorGridFunctionSpaceTag
            >::value;

        };

        //! Helper function for extracting (or building) the component name and adding
        //! the component to the VTKWriter.
        template<typename DGF, typename TreePath>
        void add_to_vtk_writer(const shared_ptr<DGF>& dgf, TreePath tp)
        {
          std::string name = name_generator(dgf->localFunctionSpace().gridFunctionSpace().name(),tp);
          switch (dgf->dataSetType())
            {
            case DGF::Output::vertexData:
              vtk_writer.addVertexData(std::make_shared<VTKGridFunctionAdapter<DGF> >(dgf,name.c_str()));
              break;
            case DGF::Output::cellData:
              vtk_writer.addCellData(std::make_shared<VTKGridFunctionAdapter<DGF> >(dgf,name.c_str()));
              break;
            default:
              DUNE_THROW(NotImplemented,"Unsupported data set type");
            }
        }

        //! Tag dispatch-based switch that creates a vector-valued function for a VectorGridFunctionSpace.
        /**
         * This version handles an actual VectorGridFunctionSpace.
         */
        template<typename LFS, typename TreePath>
        void add_vector_solution(const LFS& lfs, TreePath tp, VectorGridFunctionSpaceTag tag)
        {
          add_to_vtk_writer(std::make_shared<DGFTreeVectorFunction<LFS,Data> >(lfs,data),tp);
        }

        //! Tag dispatch-based switch that creates a vector-valued function for a VectorGridFunctionSpace.
        /**
         * This is the default version for different types of spaces that does nothing.
         */
        template<typename LFS, typename TreePath>
        void add_vector_solution(const LFS& lfs, TreePath tp, GridFunctionSpaceTag tag)
        {
          // do nothing here - not a vector space
        }

        // **********************************************************************
        // Visitor functions for adding DiscreteGridFunctions to VTKWriter
        //
        // The visitor functions contain a switch that will make them ignore
        // function spaces with a different underlying GridView type than
        // the VTKWriter.
        // This cannot happen in vanilla PDELab, but is required for MultiDomain
        // support
        // **********************************************************************

        // don't do anything if GridView types differ
        template<typename LFS, typename TreePath>
        typename std::enable_if<
          !std::is_same<
            typename LFS::Traits::GridFunctionSpace::Traits::GridView,
            typename vtk_writer_traits<VTKWriter>::GridView
            >::value
          >::type
        post(const LFS& lfs, TreePath tp)
        {
        }

        // don't do anything if GridView types differ
        template<typename LFS, typename TreePath>
        typename std::enable_if<
          !std::is_same<
            typename LFS::Traits::GridFunctionSpace::Traits::GridView,
            typename vtk_writer_traits<VTKWriter>::GridView
            >::value
          >::type
        leaf(const LFS& lfs, TreePath tp)
        {
        }

        //! Handle VectorGridFunctionSpace components in here.
        template<typename LFS, typename TreePath>
        typename std::enable_if<
          std::is_same<
            typename LFS::Traits::GridFunctionSpace::Traits::GridView,
            typename vtk_writer_traits<VTKWriter>::GridView
            >::value
          >::type
        post(const LFS& lfs, TreePath tp)
        {
          if (predicate(lfs, tp))
            add_vector_solution(lfs,tp,TypeTree::ImplementationTag<typename LFS::Traits::GridFunctionSpace>());
        }

        //! Create a standard leaf function for leaf GridFunctionSpaces.
        template<typename LFS, typename TreePath>
        typename std::enable_if<
          std::is_same<
            typename LFS::Traits::GridFunctionSpace::Traits::GridView,
            typename vtk_writer_traits<VTKWriter>::GridView
            >::value
          >::type
        leaf(const LFS& lfs, TreePath tp)
        {
          if (predicate(lfs, tp))
            add_to_vtk_writer(std::make_shared<DGFTreeLeafFunction<LFS,Data> >(lfs,data),tp);
        }


        add_solution_to_vtk_writer_visitor(VTKWriter& vtk_writer_, shared_ptr<Data> data_, const NameGenerator& name_generator_, const typename Data::Predicate& predicate_)
          : vtk_writer(vtk_writer_)
          , data(data_)
          , name_generator(name_generator_)
          , predicate(predicate_)
        {}

        VTKWriter& vtk_writer;
        shared_ptr<Data> data;
        const NameGenerator& name_generator;
        typename Data::Predicate predicate;

      };

      struct DefaultPredicate
      {
        template<typename LFS, typename TP>
        bool operator()(const LFS& lfs, TP tp) const
        {
          return true;
        }
      };

      template<typename VTKWriter, typename Data_>
      struct OutputCollector
      {

        //! Common data container (hierarchic LFS, global solution data etc.)
        typedef Data_ Data;

        typedef typename Data::GridFunctionSpace GFS;
        typedef typename Data::Vector Vector;
        typedef typename Data::Predicate Predicate;

        template<typename NameGenerator>
        OutputCollector& addSolution(const NameGenerator& name_generator)
        {

          add_solution_to_vtk_writer_visitor<VTKWriter,Data,NameGenerator> visitor(_vtk_writer,_data,name_generator,_predicate);
          TypeTree::applyToTree(_data->_lfs,visitor);
          return *this;
        }

        template<typename Factory, typename TreePath>
        OutputCollector& addCellFunction(Factory factory, TreePath tp, std::string name)
        {
          typedef typename std::remove_reference<decltype(*factory.create(_data->_lfs.child(tp),_data))>::type DGF;
          _vtk_writer.addCellData(std::make_shared<VTKGridFunctionAdapter<DGF> >(factory.create(_data->_lfs.child(tp),_data),name));
          return *this;
        }

        template<template<typename...> class Function, typename TreePath, typename... Params>
        OutputCollector& addCellFunction(TreePath tp, std::string name, Params&&... params)
        {
          using LFS = TypeTree::ChildForTreePath<typename Data::LFS,TreePath>;
          typedef Function<LFS,Data,Params...> DGF;
          _vtk_writer.addCellData(
            std::make_shared<VTKGridFunctionAdapter<DGF> >(
              std::make_shared<DGF>(
                TypeTree::child(_data->_lfs,tp)
                ),
                _data,
                std::forward<Params>(params)...
              ),
              name
            );
          return *this;
        }

        template<typename Factory, typename TreePath>
        OutputCollector& addVertexFunction(Factory factory, TreePath tp, std::string name)
        {
          typedef typename std::remove_reference<decltype(*factory.create(_data->_lfs.child(tp),_data))>::type DGF;
          _vtk_writer.addVertexData(std::make_shared<VTKGridFunctionAdapter<DGF> >(factory.create(_data->_lfs.child(tp),_data),name));
          return *this;
        }

        template<template<typename...> class Function, typename TreePath, typename... Params>
        OutputCollector& addVertexFunction(TreePath tp, std::string name, Params&&... params)
        {
          using LFS = TypeTree::ChildForTreePath<typename Data::LFS,TreePath>;
          typedef Function<LFS,Data,Params...> DGF;
          _vtk_writer.addVertexData(
            std::make_shared<VTKGridFunctionAdapter<DGF> >(
              std::make_shared<DGF>(
                TypeTree::child(_data->_lfs,tp)
                ),
                _data,
                std::forward<Params>(params)...
              ),
              name
            );
          return *this;
        }

        OutputCollector(VTKWriter& vtk_writer, const shared_ptr<Data>& data, const Predicate& predicate = Predicate())
          : _vtk_writer(vtk_writer)
          , _data(data)
          , _predicate(predicate)
        {}

        VTKWriter& _vtk_writer;
        shared_ptr<Data> _data;
        Predicate _predicate;

      };

    } // namespace vtk

    template<typename VTKWriter,
             typename GFS,
             typename X,
             typename NameGenerator = vtk::DefaultFunctionNameGenerator,
             typename Predicate = vtk::DefaultPredicate>
    vtk::OutputCollector<
      VTKWriter,
      vtk::DGFTreeCommonData<GFS,X,Predicate>
      >
    addSolutionToVTKWriter(VTKWriter& vtk_writer,
                           const GFS& gfs,
                           const X& x,
                           const NameGenerator& name_generator = vtk::defaultNameScheme(),
                           const Predicate& predicate = Predicate())
    {
      typedef vtk::DGFTreeCommonData<GFS,X,Predicate> Data;
      vtk::OutputCollector<VTKWriter,Data> collector(vtk_writer,std::make_shared<Data>(gfs,x),predicate);
      collector.addSolution(name_generator);
      return collector;
    }


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_VTK_HH
