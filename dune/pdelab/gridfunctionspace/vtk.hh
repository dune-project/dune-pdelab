#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_VTK_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_VTK_HH

#include <vector>
#include <sstream>
#include <functional>

#include <dune/common/exceptions.hh>
#include <dune/common/std/functional.hh>

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

      template<typename LFS, typename Data, typename GV>
      class DGFTreeLeafFunction;

      template<typename LFS, typename Data, typename GV>
      class DGFTreeVectorFunction;

      template<typename VTKWriter, typename Data>
      struct OutputCollector;

      /**
       * @brief Helper class for common data of a DGFTree.
       *
       * @tparam GFS      GridFunctionSpace type
       * @tparam X        Vector type
       * @tparam Pred     Predicate for deciding which nodes will be written
       * @tparam GV       GridView type used for entity binding
       * @tparam ET       Entity transformation
       *
       * @note If GV is not the same as the GFS grid view, the entity transformation
       * must be able to take an entity from the GFS grid view and return an entity
       * from the VTK grid view. For instance, for a multidomain grid, the following
       * etity transformation would make possible to write data from GFS with a host
       * domain grid view into a one of the subdomains of the multidomain grid:
       * @code{.cpp}
       *   auto etity_transformation = [&](auto e){return grid->multiDomainEntity(e);};
       * @endcode
       */
      template<typename GFS, typename X, typename Pred, typename GV = typename GFS::Traits::GridView, typename ET = Std::identity>
      struct DGFTreeCommonData
      {

        template<typename LFS, typename Data,typename _GV>
        friend class DGFTreeLeafFunction;

        template<typename LFS, typename Data, typename _GV>
        friend class DGFTreeVectorFunction;

        template<typename, typename>
        friend struct OutputCollector;

        using LFS = LocalFunctionSpace<GFS>;
        using LFSCache = LFSIndexCache<LFS>;
        using XView = typename X::template ConstLocalView<LFSCache>;
        using XLocalVector = LocalVector<typename X::ElementType>;
        using Cell = typename GV::template Codim<0>::Entity;
        using GFSCell = decltype(std::declval<ET>()(std::declval<const Cell&>()));
        using EntitySet = typename GFS::Traits::EntitySet;
        using IndexSet = typename EntitySet::Traits::IndexSet;
        using size_type = typename IndexSet::IndexType;

        //! cache entities only when they are a temporary object
        static constexpr bool cache_entity = not std::is_reference_v<GFSCell>;
        using CellCache = std::decay_t<GFSCell>;

        static const auto dim = GV::dimension;

      public:

        typedef GFS GridFunctionSpace;
        typedef X Vector;
        typedef Pred Predicate;
        typedef GV GridView;
        typedef ET EntityTransformation;

        /**
         * @brief Construct a new DGFTreeCommonData object
         *
         * @param gfs             grid function space
         * @param x               coefficient vector associated with the gfs
         * @param gv              grid view used for entity binding
         * @param entity_transf   entity transformation to the grid view used in the gfs
         */
        DGFTreeCommonData(const GFS& gfs, const X& x, const GV& gv, const ET& entity_transf)
          : _entity_transf(entity_transf)
          , _gv(gv)
          , _lfs(gfs)
          , _lfs_cache(_lfs)
          , _x_view(x)
          , _x_local(_lfs.maxSize())
          , _index_set(gfs.entitySet().indexSet())
          , _current_cell_index(std::numeric_limits<size_type>::max())
        {}

        /**
         * @brief Construct a new DGFTreeCommonData object
         *
         * @param gfs             grid function space
         * @param x               coefficient vector associated with the gfs
         * @param gv              grid view used for entity binding
         *
         * @note Only available if entity transformation is default constructible
         */
        template<class T = int, class = std::enable_if_t<std::is_default_constructible_v<ET>,T>>
        DGFTreeCommonData(const GFS& gfs, const X& x, const GV& gv)
          : DGFTreeCommonData(gfs,x,gv,ET{})
        {}

        //! Returns the grid view used for data binding
        const GridView& gridView() const {return _gv;}

      private:

        //! cache entity when necessary
        template<class E>
        const CellCache& cacheCell(const E& cell)
        {
          if (cache_entity)
          {
            _cached_cell = cell;
            return _cached_cell;
          }
          else
            return cell;
        }

      public:

        void bind(const Cell& cell)
        {
          // transform cell to the gfs cell and cache it if necessary
          const CellCache& gfs_cell = cacheCell(_entity_transf(cell));

          // avoid duplicated binds
          auto cell_index = _index_set.uniqueIndex(gfs_cell);
          if (_current_cell_index == cell_index)
            return;

          // actual databinding
          _lfs.bind(gfs_cell);
          _lfs_cache.update();
          _x_view.bind(_lfs_cache);
          _x_view.read(_x_local);
          _x_view.unbind();
          _current_cell_index = cell_index;
        }

        ET _entity_transf;
        GV _gv;
        CellCache _cached_cell;
        LFS _lfs;
        LFSCache _lfs_cache;
        XView _x_view;
        XLocalVector _x_local;
        const IndexSet& _index_set;
        size_type _current_cell_index;

      };



      template<typename LFS, typename Data, typename GV>
      class DGFTreeLeafFunction
        : public GridFunctionBase<GridFunctionTraits<
                                         GV,
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
                                       DGFTreeLeafFunction<LFS,Data,GV>
                                       >
      {

        typedef BasisInterfaceSwitch<
          typename FiniteElementInterfaceSwitch<
            typename LFS::Traits::FiniteElement
            >::Basis
          > BasisSwitch;

        typedef GridFunctionBase<
          GridFunctionTraits<
            GV,
            typename BasisSwitch::RangeField,
            BasisSwitch::dimRange,
            typename BasisSwitch::Range
            >,
          DGFTreeLeafFunction<LFS,Data,GV>
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
        const typename Traits::GridViewType& getGridView() const
        {
          return _data->gridView();
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



      template<typename LFS, typename Data, typename GV>
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
                                       DGFTreeVectorFunction<LFS,Data,GV>
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
          DGFTreeVectorFunction<LFS,Data,GV>
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
          return _data->gridView();
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

        using GV = typename Data::GridView;

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
          add_to_vtk_writer(std::make_shared<DGFTreeVectorFunction<LFS,Data,GV> >(lfs,data),tp);
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

        //! Handle VectorGridFunctionSpace components in here.
        template<typename LFS, typename TreePath>
        void post(const LFS& lfs, TreePath tp)
        {
          if (predicate(lfs, tp))
            add_vector_solution(lfs,tp,TypeTree::ImplementationTag<typename LFS::Traits::GridFunctionSpace>());
        }

        //! Create a standard leaf function for leaf GridFunctionSpaces.
        template<typename LFS, typename TreePath>
        void leaf(const LFS& lfs, TreePath tp)
        {
          if (predicate(lfs, tp))
            add_to_vtk_writer(std::make_shared<DGFTreeLeafFunction<LFS,Data,GV> >(lfs,data),tp);
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
        typedef typename Data::GridView GridView;

        template<typename LFS, typename NameGenerator>
        OutputCollector& addSolution(const LFS& lfs, const NameGenerator& name_generator)
        {

          add_solution_to_vtk_writer_visitor<VTKWriter,Data,NameGenerator> visitor(_vtk_writer,_data,name_generator,_predicate);
          TypeTree::applyToTree(lfs,visitor);
          return *this;
        }

        template<typename NameGenerator>
        OutputCollector& addSolution(const NameGenerator& name_generator)
        {
          return addSolution(_data->_lfs,name_generator);
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
    auto
    addSolutionToVTKWriter(VTKWriter& vtk_writer,
                           const GFS& gfs,
                           const X& x,
                           const NameGenerator& name_generator = vtk::defaultNameScheme(),
                           const Predicate& predicate = Predicate())
    {
      typedef vtk::DGFTreeCommonData<GFS,X,Predicate> Data;
      vtk::OutputCollector<VTKWriter,Data> collector(vtk_writer,std::make_shared<Data>(gfs,x,gfs.gridView()),predicate);
      collector.addSolution(name_generator);
      return collector;
    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_VTK_HH
