// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_ADAPTIVITY_HH
#define DUNE_PDELAB_ADAPTIVITY_HH

#include<dune/common/exceptions.hh>

#include<limits>
#include<vector>
#include<map>
#include<unordered_map>
#include<dune/common/dynmatrix.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/localfunctionspace.hh>

#include<dune/pdelab/common/function.hh>
// for InterpolateBackendStandard
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
// for intersectionoperator
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/flags.hh>

#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

namespace Dune {
  namespace PDELab {


    template<typename GFS>
    struct LeafOffsetCache
    {

      typedef typename GFS::Traits::GridView::template Codim<0>::Entity Cell;
      typedef LocalFunctionSpace<GFS> LFS;

        // we need an additional entry because we store offsets and we also want the
        // offset after the last leaf for size calculations
      typedef array<std::size_t,TypeTree::TreeInfo<GFS>::leafCount + 1> LeafOffsets;

      const LeafOffsets& operator[](GeometryType gt) const
      {
        const LeafOffsets& leaf_offsets = _leaf_offset_cache[GlobalGeometryTypeIndex::index(gt)];
        // make sure we have data for this geometry type
        assert(leaf_offsets.back() > 0);
        return leaf_offsets;
      }

      void update(const Cell& e)
      {
        LeafOffsets& leaf_offsets = _leaf_offset_cache[GlobalGeometryTypeIndex::index(e.type())];
        if (leaf_offsets.back() == 0)
          {
            _lfs.bind(e);
            extract_lfs_leaf_sizes(_lfs,leaf_offsets.begin()+1);
            // convert to offsets
            std::partial_sum(leaf_offsets.begin(),leaf_offsets.end(),leaf_offsets.begin());
            // sanity check
            assert(leaf_offsets.back() == _lfs.size());
          }
      }

      explicit LeafOffsetCache(const GFS& gfs)
        : _lfs(gfs)
        , _leaf_offset_cache(GlobalGeometryTypeIndex::size(Cell::dimension))
      {}

      LFS _lfs;
      std::vector<LeafOffsets> _leaf_offset_cache;

    };


    namespace {

      template<typename MassMatrices,typename Cell>
      struct inverse_mass_matrix_calculator
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        static const int dim  = Cell::Geometry::dimension;
        typedef std::size_t size_type;
        typedef typename MassMatrices::value_type MassMatrix;
        typedef typename MassMatrix::field_type DF;
        typedef typename Dune::QuadratureRule<DF,dim>::const_iterator QRIterator;

        template<typename GFS, typename TreePath>
        void leaf(const GFS& gfs, TreePath treePath)
        {
          typedef typename GFS::Traits::FiniteElementMap FEM;
          const FEM& fem = gfs.finiteElementMap();
          typedef typename FEM::Traits::FiniteElement FiniteElement;
          const FiniteElement& fe = fem.find(_element);
          size_type local_size = fe.localBasis().size();

          MassMatrix& mass_matrix = _mass_matrices[_leaf_index];
          mass_matrix.resize(local_size,local_size);


          std::vector<typename FiniteElement::Traits::LocalBasisType::Traits::RangeType> phi;
          phi.resize(std::max(phi.size(),local_size));

          for (QRIterator it = _quadrature_rule.begin(); it != _quadrature_rule.end(); ++it)
            {
              //std::fill(yVector.begin(),yVector.end(),typename Traits::RangeType(0));
              fe.localBasis().evaluateFunction(it->position(),phi);
              const DF factor = it->weight();

              for (int i = 0; i < local_size; ++i)
                for (int j = 0; j < local_size; ++j)
                  mass_matrix[i][j] += phi[i] * phi[j] * factor;
            }

          mass_matrix.invert();
          ++_leaf_index;

        }

        inverse_mass_matrix_calculator(MassMatrices& mass_matrices, const Cell& element, size_type intorder)
          : _element(element)
          , _mass_matrices(mass_matrices)
          , _quadrature_rule(QuadratureRules<DF,dim>::rule(element.type(),intorder))
          , _leaf_index(0)
        {}

        const Cell& _element;
        MassMatrices& _mass_matrices;
        const QuadratureRule<DF,dim>& _quadrature_rule;
        size_type _leaf_index;

      };

    } // anonymous namespace


    /*! @class L2Projection
     *
     * @brief @todo
     *
     * @tparam GFS Type of ansatz space
     * @tparam U   Container class for the solution
     */
    template<class GFS, class U>
    class L2Projection
    {
      typedef typename GFS::Traits::GridViewType::Grid Grid;
      typedef typename Grid::template Codim<0>::Entity Element;
      typedef LocalFunctionSpace<GFS> LFS;
      typedef typename U::ElementType DF;

    public:

      typedef DynamicMatrix<typename U::ElementType> MassMatrix;
      typedef array<MassMatrix,TypeTree::TreeInfo<GFS>::leafCount> MassMatrices;

      /*! @brief The constructor.
       *
       * @todo Doc params!
       */
      explicit L2Projection(const GFS& gfs, int intorder = 2)
        : _gfs(gfs)
        , _intorder(intorder)
        , _inverse_mass_matrices(GlobalGeometryTypeIndex::size(Element::dimension))
      {}

      /*! @brief Calculate the inverse local mass matrix, used in the local L2 projection
       *
       * @todo Doc template params
       * @todo Doc params
       */
      const MassMatrices& inverseMassMatrices(const Element& e)
      {
        const GeometryType gt = e.geometry().type();
        MassMatrices& inverse_mass_matrices = _inverse_mass_matrices[GlobalGeometryTypeIndex::index(gt)];
        // if the matrix isn't empty, it has already been cached
        if (inverse_mass_matrices[0].N() > 0)
          return inverse_mass_matrices;

        inverse_mass_matrix_calculator<MassMatrices,Element> calculate_mass_matrices(
          inverse_mass_matrices,
          e,
          _intorder
          );

        TypeTree::applyToTree(_gfs,calculate_mass_matrices);

        return inverse_mass_matrices;
      }

    private:

      const GFS& _gfs;
      int _intorder;
      std::vector<MassMatrices> _inverse_mass_matrices;
    };


    template<typename GFS, typename DOFVector, typename TransferMap>
    struct backup_visitor
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      typedef LocalFunctionSpace<GFS> LFS;
      typedef LFSIndexCache<LFS> LFSCache;
      typedef Dune::PDELab::LeafOffsetCache<GFS> LeafOffsetCache;

      typedef typename GFS::Traits::GridView::Grid::LocalIdSet IDSet;
      typedef typename GFS::Traits::GridView::template Codim<0>::Entity Cell;
      typedef typename GFS::Traits::GridView::template Codim<0>::EntityPointer CellPointer;
      typedef typename Cell::Geometry Geometry;
      static const int dim = Geometry::dimension;
      typedef typename Cell::HierarchicIterator HierarchicIterator;
      typedef typename DOFVector::ElementType RF;
      typedef typename TransferMap::mapped_type LocalDOFVector;


      typedef L2Projection<typename LFS::Traits::GridFunctionSpace,DOFVector> Projection;
      typedef typename Projection::MassMatrices MassMatrices;
      typedef typename Projection::MassMatrix MassMatrix;

      typedef std::size_t size_type;
      typedef typename GFS::Traits::GridView::ctype DF;

      template<typename LFSLeaf, typename TreePath>
      void leaf(const LFSLeaf& leaf_lfs, TreePath treePath)
      {

        typedef typename LFSLeaf::Traits::GridFunctionSpace::Traits::FiniteElementMap FEM;
        typedef typename FEM::Traits::FiniteElement FE;
        const FEM& fem = leaf_lfs.gridFunctionSpace().finiteElementMap();
        size_type fine_offset = _leaf_offset_cache[_current->type()][_leaf_index];
        size_type coarse_offset = _leaf_offset_cache[_ancestor->type()][_leaf_index];

        typedef typename FE::Traits::LocalBasisType::Traits::RangeType Range;

        const MassMatrix& inverse_mass_matrix = _projection.inverseMassMatrices(*_element)[_leaf_index];

        std::vector<Range> coarse_phi;
        std::vector<Range> fine_phi;

        Geometry fine_geometry = _current->geometry();
        Geometry coarse_geometry = _ancestor->geometry();

        const QuadratureRule<DF,dim>& rule = QuadratureRules<DF,dim>::rule(_current->type(),_int_order);
        // iterate over quadrature points
        for (typename QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            typename Geometry::LocalCoordinate coarse_local = coarse_geometry.local(fine_geometry.global(it->position()));
            const FE* fe = &fem.find(*_current);
            fe->localBasis().evaluateFunction(it->position(),fine_phi);
            fe = &fem.find(*_ancestor);
            fe->localBasis().evaluateFunction(coarse_local,coarse_phi);
            const DF factor = it->weight()
              * fine_geometry.integrationElement(it->position())
              / coarse_geometry.integrationElement(coarse_local);

            Range val(0.0);
            for (size_type i = 0; i < fine_phi.size(); ++i)
              {
                val.axpy(_u_fine[fine_offset + i],fine_phi[i]);
              }

            for (size_type i = 0; i < coarse_phi.size(); ++i)
              {
                Range x(0.0);
                for (size_type j = 0; j < inverse_mass_matrix.M(); ++j)
                  x.axpy(inverse_mass_matrix[i][j],coarse_phi[j]);
                (*_u_coarse)[coarse_offset + i] += factor * (x * val);
              }
          }

        ++_leaf_index;
      }

      void operator()(const Cell& element)
      {
        _element = &element;

        _lfs.bind(element);
        _lfs_cache.update();
        _u_view.bind(_lfs_cache);
        _u_coarse = &_transfer_map[_id_set.id(element)];
        _u_coarse->resize(_lfs.size());
        _u_view.read(*_u_coarse);
        _u_view.unbind();

        _leaf_offset_cache.update(element);

        size_type max_level = _lfs.gridFunctionSpace().gridView().grid().maxLevel();

        CellPointer ancestor(element);
        while (ancestor->mightVanish())
          {
            // work around UG bug!
            if (!ancestor->hasFather())
              break;

            ancestor = ancestor->father();
            _ancestor = &(*ancestor);

            _u_coarse = &_transfer_map[_id_set.id(*_ancestor)];
            // don't project more than once
            if (_u_coarse->size() > 0)
              continue;
            _u_coarse->resize(_leaf_offset_cache[_ancestor->type()].back());
            std::fill(_u_coarse->begin(),_u_coarse->end(),RF(0));

            for (HierarchicIterator hit = _ancestor->hbegin(max_level),
                   hend = _ancestor->hend(max_level);
                 hit != hend;
                 ++hit)
              {
                // only evaluate on entities with data
                if (hit->isLeaf())
                  {
                    _current = &(*hit);
                    // reset leaf_index for next run over tree
                    _leaf_index = 0;
                    // load data
                    _lfs.bind(*hit);
                    _leaf_offset_cache.update(*hit);
                    _lfs_cache.update();
                    _u_view.bind(_lfs_cache);
                    _u_fine.resize(_lfs_cache.size());
                    _u_view.read(_u_fine);
                    _u_view.unbind();
                    // do projection on all leafs
                    TypeTree::applyToTree(_lfs,*this);
                  }
              }
          }
      }

      backup_visitor(const GFS& gfs,
                     Projection& projection,
                     const DOFVector& u,
                     LeafOffsetCache& leaf_offset_cache,
                     TransferMap& transfer_map,
                     std::size_t int_order = 2)
        : _lfs(gfs)
        , _lfs_cache(_lfs)
        , _id_set(gfs.gridView().grid().localIdSet())
        , _element(nullptr)
        , _ancestor(nullptr)
        , _current(nullptr)
        , _projection(projection)
        , _u_view(u)
        , _transfer_map(transfer_map)
        , _u_coarse(nullptr)
        , _leaf_offset_cache(leaf_offset_cache)
        , _int_order(int_order)
        , _leaf_index(0)
      {}

      LFS _lfs;
      LFSCache _lfs_cache;
      const IDSet& _id_set;
      const Cell* _element;
      const Cell* _ancestor;
      const Cell* _current;
      Projection& _projection;
      typename DOFVector::template ConstLocalView<LFSCache> _u_view;
      TransferMap& _transfer_map;
      LocalDOFVector* _u_coarse;
      LeafOffsetCache& _leaf_offset_cache;
      size_type _int_order;
      size_type _leaf_index;
      LocalDOFVector _u_fine;

    };



    template<typename GFS, typename DOFVector, typename CountVector>
    struct replay_visitor
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      typedef LocalFunctionSpace<GFS> LFS;
      typedef LFSIndexCache<LFS> LFSCache;
      typedef Dune::PDELab::LeafOffsetCache<GFS> LeafOffsetCache;

      typedef typename LFS::Traits::GridFunctionSpace::Traits::GridView::template Codim<0>::Entity Cell;
      typedef typename Cell::Geometry Geometry;
      typedef typename DOFVector::ElementType RF;
      typedef std::vector<RF> LocalDOFVector;
      typedef std::vector<typename CountVector::ElementType> LocalCountVector;

      typedef std::size_t size_type;

      template<typename FiniteElement>
      struct coarse_function
      {

        template<typename X, typename Y>
        void evaluate(const X& x, Y& y) const
        {
          _phi.resize(_finite_element.localBasis().size());
          _finite_element.localBasis().evaluateFunction(_coarse_geometry.local(_fine_geometry.global(x)),_phi);
          y = 0;
          for (size_type i = 0; i < _phi.size(); ++i)
            y.axpy(_dofs[_offset + i],_phi[i]);
        }

        coarse_function(const FiniteElement& finite_element, Geometry coarse_geometry, Geometry fine_geometry, const LocalDOFVector& dofs, size_type offset)
          : _finite_element(finite_element)
          , _coarse_geometry(coarse_geometry)
          , _fine_geometry(fine_geometry)
          , _dofs(dofs)
          , _offset(offset)
        {}

        const FiniteElement& _finite_element;
        Geometry _coarse_geometry;
        Geometry _fine_geometry;
        const LocalDOFVector& _dofs;
        mutable std::vector<typename FiniteElement::Traits::LocalBasisType::Traits::RangeType> _phi;
        size_type _offset;

      };


      template<typename LeafLFS, typename TreePath>
      void leaf(const LeafLFS& leaf_lfs, TreePath treePath)
      {

        typedef typename LeafLFS::Traits::GridFunctionSpace::Traits::FiniteElementMap FEM;
        const FEM& fem = leaf_lfs.gridFunctionSpace().finiteElementMap();
        size_type element_offset = _leaf_offset_cache[_element->type()][_leaf_index];
        size_type ancestor_offset = _leaf_offset_cache[_ancestor->type()][_leaf_index];

        coarse_function<typename FEM::Traits::FiniteElement> f(fem.find(*_ancestor),_ancestor->geometry(),_element->geometry(),*_u_coarse,ancestor_offset);
        const typename FEM::Traits::FiniteElement& fe = fem.find(*_element);

        _u_tmp.resize(fe.localBasis().size());
        std::fill(_u_tmp.begin(),_u_tmp.end(),RF(0.0));
        fe.localInterpolation().interpolate(f,_u_tmp);
        std::copy(_u_tmp.begin(),_u_tmp.end(),_u_fine.begin() + element_offset);

        ++_leaf_index;
      }

      void operator()(const Cell& element, const Cell& ancestor, const LocalDOFVector& u_coarse)
      {
        _element = &element;
        _ancestor = &ancestor;
        _u_coarse = &u_coarse;
        _lfs.bind(*_element);
        _leaf_offset_cache.update(*_element);
        _lfs_cache.update();
        _u_view.bind(_lfs_cache);

        // test identity using ids
        if (_lfs.gridFunctionSpace().gridView().grid().localIdSet().id(element) ==
            _lfs.gridFunctionSpace().gridView().grid().localIdSet().id(ancestor))
          {
            // no interpolation necessary, just copy the saved data
            _u_view.add(*_u_coarse);
          }
        else
          {
            _u_fine.resize(_lfs_cache.size());
            std::fill(_u_fine.begin(),_u_fine.end(),RF(0));
            _leaf_index = 0;
            TypeTree::applyToTree(_lfs,*this);
            _u_view.add(_u_fine);
          }
        _u_view.commit();

        _uc_view.bind(_lfs_cache);
        _counts.resize(_lfs_cache.size(),1);
        _uc_view.add(_counts);
        _uc_view.commit();
      }

      replay_visitor(const GFS& gfs, DOFVector& u, CountVector& uc, LeafOffsetCache& leaf_offset_cache)
        : _lfs(gfs)
        , _lfs_cache(_lfs)
        , _element(nullptr)
        , _ancestor(nullptr)
        , _u_view(u)
        , _uc_view(uc)
        , _leaf_offset_cache(leaf_offset_cache)
        , _leaf_index(0)
      {}

      LFS _lfs;
      LFSCache _lfs_cache;
      const Cell* _element;
      const Cell* _ancestor;
      typename DOFVector::template LocalView<LFSCache> _u_view;
      typename CountVector::template LocalView<LFSCache> _uc_view;
      const LocalDOFVector* _u_coarse;
      LeafOffsetCache& _leaf_offset_cache;
      size_type _leaf_index;
      LocalDOFVector _u_fine;
      LocalDOFVector _u_tmp;
      LocalCountVector _counts;

    };


    /*! @class GridAdaptor
     *
     * @brief Class for automatic adaptation of the grid.
     *
     *        The GridAdaptor capsules the act of deciding which Elems to refine and coarsen,
     *        adapting the grid, and transfering the solution from the old grid to the new one.
     *        Currrently this only works for scalar solutions.
     *
     * @tparam Grid       Type of the grid we want to adapt
     * @tparam GFSU       Type of ansatz space, we need to update it after adaptation
     * @tparam U          Container class of the solution
     * @tparam Projection Projection used when Elems vanish
     */
    template<class Grid, class GFSU, class U, class Projection>
    class GridAdaptor
    {
      typedef typename Grid::LeafGridView LeafGridView;
      typedef typename LeafGridView::template Codim<0>
      ::template Partition<Dune::Interior_Partition>::Iterator LeafIterator;
      typedef typename Grid::template Codim<0>::Entity Element;
      typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
      typedef typename Grid::LocalIdSet IDSet;
      typedef typename IDSet::IdType ID;

    public:
      typedef std::unordered_map<ID,std::vector<typename U::ElementType> > MapType;


      /*! @brief The constructor.
       *
       * @param grid_       The grid we want to adapt
       * @param gfsu_       The ansatz space, we need to update it
       * @param projection_ The Projection used when Elems vanish
       */
      explicit GridAdaptor(const GFSU& gfs)
        : _leaf_offset_cache(gfs)
      {}

      /* @brief @todo
       *
       * @param[in]  u           The solution that will be saved
       * @param[out] transferMap The map containing the solution during adaptation
       */
      void backupData(Grid& grid, GFSU& gfsu, Projection& projection, U& u, MapType& transfer_map)
      {
        typedef backup_visitor<GFSU,U,MapType> Visitor;

        Visitor visitor(gfsu,projection,u,_leaf_offset_cache,transfer_map);

        // iterate over all elems
        LeafGridView leafView = grid.leafGridView();
        for (LeafIterator it = leafView.template begin<0,Dune::Interior_Partition>();
             it!=leafView.template end<0,Dune::Interior_Partition>(); ++it)
          {
            visitor(*it);
          }
      }

      /* @brief @todo
       *
       * @param[out] u           The solution after adaptation
       * @param[in]  transferMap The map that contains the information for the rebuild of u
       */
      void replayData(Grid& grid, GFSU& gfsu, Projection& projection, U& u, const MapType& transfer_map)
      {
        const IDSet& id_set = grid.globalIdSet();

        typedef typename BackendVectorSelector<GFSU,int>::Type CountVector;
        CountVector uc(gfsu,0);

        typedef replay_visitor<GFSU,U,CountVector> Visitor;
        Visitor visitor(gfsu,u,uc,_leaf_offset_cache);

        // iterate over all elems
        LeafGridView leafView = grid.leafGridView();
        for (LeafIterator it = leafView.template begin<0,Dune::Interior_Partition>();
             it!=leafView.template end<0,Dune::Interior_Partition>(); ++it)
          {
            const Element& e = *it;

            ElementPointer ancestor(e);

            typename MapType::const_iterator map_it;
            while ((map_it = transfer_map.find(id_set.id(*ancestor))) == transfer_map.end())
              {
                if (!ancestor->hasFather())
                  DUNE_THROW(Exception,
                             "transferMap of GridAdaptor didn't contain ancestor of element with id " << id_set.id(*ancestor));
                ancestor = ancestor->father();
              }

            visitor(e,*ancestor,map_it->second);
          }

        typedef Dune::PDELab::AddDataHandle<GFSU,U> DOFHandle;
        DOFHandle addHandle1(gfsu,u);
        leafView.communicate (addHandle1,
                              Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
        typedef Dune::PDELab::AddDataHandle<GFSU,CountVector> CountHandle;
        CountHandle addHandle2(gfsu,uc);
        leafView.communicate (addHandle2,
                              Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);

        // normalize multiple-interpolated DOFs by taking the arithmetic average
        typename CountVector::iterator ucit = uc.begin();
        for (typename U::iterator uit = u.begin(), uend = u.end(); uit != uend; ++uit, ++ucit)
          (*uit) /= ((*ucit) > 0 ? (*ucit) : 1.0);
      }

    private:

      LeafOffsetCache<GFSU> _leaf_offset_cache;

    };

    /*! grid adaptation as a function
     *
     * @brief adapt a grid, corresponding function space and solution vectors
     *
     * Assumes that the grid's elements have been marked for refinement and coarsening appropriately before
     *
     * @tparam Grid       Type of the grid we want to adapt
     * @tparam GFS        Type of ansatz space, we need to update it after adaptation
     * @tparam X          Container class for DOF vectors
     */
    template<class Grid, class GFS, class X>
    void adapt_grid (Grid& grid, GFS& gfs, X& x1, int int_order)
    {
      typedef L2Projection<GFS,X> Projection;
      Projection projection(gfs,int_order);

      GridAdaptor<Grid,GFS,X,Projection> grid_adaptor(gfs);

      // prepare the grid for refinement
      grid.preAdapt();

      // save u
      typename GridAdaptor<Grid,GFS,X,Projection>::MapType transferMap1;
      grid_adaptor.backupData(grid,gfs,projection,x1,transferMap1);

      // adapt the grid
      grid.adapt();

      // update the function spaces
      gfs.update();

      // reset u
      x1 = X(gfs,0.0);
      grid_adaptor.replayData(grid,gfs,projection,x1,transferMap1);

      // clean up
      grid.postAdapt();
    }

    /*! grid adaptation as a function
     *
     * @brief adapt a grid, corresponding function space and solution vectors
     *
     * Assumes that the grid's elements have been marked for refinement and coarsening appropriately before
     *
     * @tparam Grid       Type of the grid we want to adapt
     * @tparam GFS        Type of ansatz space, we need to update it after adaptation
     * @tparam X          Container class for DOF vectors
     * @tparam Projection Projection used when Elems vanish
     */
    template<class Grid, class GFS, class X>
    void adapt_grid (Grid& grid, GFS& gfs, X& x1, X& x2, int int_order)
    {
      typedef L2Projection<GFS,X> Projection;
      Projection projection(gfs,int_order);

      GridAdaptor<Grid,GFS,X,Projection> grid_adaptor(gfs);

      // prepare the grid for refinement
      grid.preAdapt();

      // save solution
      typename GridAdaptor<Grid,GFS,X,Projection>::MapType transferMap1;
      grid_adaptor.backupData(grid,gfs,projection,x1,transferMap1);
      typename GridAdaptor<Grid,GFS,X,Projection>::MapType transferMap2;
      grid_adaptor.backupData(grid,gfs,projection,x2,transferMap2);

      // adapt the grid
      grid.adapt();

      // update the function spaces
      gfs.update();

      // interpolate solution
      x1 = X(gfs,0.0);
      grid_adaptor.replayData(grid,gfs,projection,x1,transferMap1);
      x2 = X(gfs,0.0);
      grid_adaptor.replayData(grid,gfs,projection,x2,transferMap2);

      // clean up
      grid.postAdapt();
    }

    // deprecated versions which always force the mass matrix integration order to 2
    // function attributes are only allowed on function declarations, not defitions, so we have to do the double
    // dance of first declaring and then immediately defining those functions...
    template<class Grid, class GFS, class X>
    void adapt_grid (Grid& grid, GFS& gfs, X& x1) DUNE_DEPRECATED_MSG("Please use the version of adapt_grid() that explicity specifies the integration order instead");

    template<class Grid, class GFS, class X>
    void adapt_grid (Grid& grid, GFS& gfs, X& x1, X& x2) DUNE_DEPRECATED_MSG("Please use the version of adapt_grid() that explicity specifies the integration order instead");

    template<class Grid, class GFS, class X>
    void adapt_grid (Grid& grid, GFS& gfs, X& x1)
    {
      adapt_grid(grid,gfs,x1,2);
    }

    template<class Grid, class GFS, class X>
    void adapt_grid (Grid& grid, GFS& gfs, X& x1, X& x2)
    {
      adapt_grid(grid,gfs,x1,x2,2);
    }



    template<typename T>
    void error_fraction(const T& x, typename T::ElementType alpha, typename T::ElementType beta,
                        typename T::ElementType& eta_alpha, typename T::ElementType& eta_beta, int verbose=0)
    {
      if (verbose>0)
        std::cout << "+++ error fraction: alpha=" << alpha << " beta=" << beta << std::endl;
      const int steps=20; // max number of bisection steps
      typedef typename T::ElementType NumberType;
      NumberType total_error = x.one_norm();
      NumberType max_error = x.infinity_norm();
      NumberType eta_alpha_left = 0.0;
      NumberType eta_alpha_right = max_error;
      NumberType eta_beta_left = 0.0;
      NumberType eta_beta_right = max_error;
      for (int j=1; j<=steps; j++)
        {
          eta_alpha = 0.5*(eta_alpha_left+eta_alpha_right);
          eta_beta = 0.5*(eta_beta_left+eta_beta_right);
          NumberType sum_alpha=0.0;
          NumberType sum_beta=0.0;
          unsigned int alpha_count = 0;
          unsigned int beta_count = 0;
          for (typename T::const_iterator it = x.begin(),
                 end = x.end();
               it != end;
               ++it)
            {
              if (*it >=eta_alpha) { sum_alpha += *it; alpha_count++;}
              if (*it < eta_beta) { sum_beta += *it; beta_count++;}
            }
          if (verbose>1)
            {
              std::cout << "+++ " << j << " eta_alpha=" << eta_alpha << " alpha_fraction=" << sum_alpha/total_error
                        << " elements: " << alpha_count << " of " << x.N() << std::endl;
              std::cout << "+++ " << j << " eta_beta=" << eta_beta << " beta_fraction=" << sum_beta/total_error
                        << " elements: " << beta_count << " of " << x.N() << std::endl;
            }
          if (std::abs(alpha-sum_alpha/total_error) <= 0.01 && std::abs(beta-sum_beta/total_error) <= 0.01) break;
          if (sum_alpha>alpha*total_error)
            eta_alpha_left = eta_alpha;
          else
            eta_alpha_right = eta_alpha;
          if (sum_beta>beta*total_error)
            eta_beta_right = eta_beta;
          else
            eta_beta_left = eta_beta;
        }
      if (verbose>0)
        {
          std::cout << "+++ refine_threshold=" << eta_alpha
                    << " coarsen_threshold=" << eta_beta << std::endl;
        }
    }


    template<typename T>
    void element_fraction(const T& x, typename T::ElementType alpha, typename T::ElementType beta,
                          typename T::ElementType& eta_alpha, typename T::ElementType& eta_beta, int verbose=0)
    {
      const int steps=20; // max number of bisection steps
      typedef typename T::ElementType NumberType;
      NumberType total_error =x.N();
      NumberType max_error = x.infinity_norm();
      NumberType eta_alpha_left = 0.0;
      NumberType eta_alpha_right = max_error;
      NumberType eta_beta_left = 0.0;
      NumberType eta_beta_right = max_error;
      for (int j=1; j<=steps; j++)
        {
          eta_alpha = 0.5*(eta_alpha_left+eta_alpha_right);
          eta_beta = 0.5*(eta_beta_left+eta_beta_right);
          NumberType sum_alpha=0.0;
          NumberType sum_beta=0.0;
          unsigned int alpha_count = 0;
          unsigned int beta_count = 0;

          for (typename T::const_iterator it = x.begin(),
                 end = x.end();
               it != end;
               ++it)
            {
              if (*it>=eta_alpha) { sum_alpha += 1.0; alpha_count++;}
              if (*it< eta_beta) { sum_beta +=1.0; beta_count++;}
            }
          if (verbose>1)
            {
              std::cout << j << " eta_alpha=" << eta_alpha << " alpha_fraction=" << sum_alpha/total_error
                        << " elements: " << alpha_count << " of " << x.N() << std::endl;
              std::cout << j << " eta_beta=" << eta_beta << " beta_fraction=" << sum_beta/total_error
                        << " elements: " << beta_count << " of " << x.N() << std::endl;
            }
          if (std::abs(alpha-sum_alpha/total_error) <= 0.01 && std::abs(beta-sum_beta/total_error) <= 0.01) break;
          if (sum_alpha>alpha*total_error)
            eta_alpha_left = eta_alpha;
          else
            eta_alpha_right = eta_alpha;
          if (sum_beta>beta*total_error)
            eta_beta_right = eta_beta;
          else
            eta_beta_left = eta_beta;
        }
      if (verbose>0)
        {
          std::cout << "+++ refine_threshold=" << eta_alpha
                    << " coarsen_threshold=" << eta_beta << std::endl;
        }
    }

    /** Compute error distribution
     */
    template<typename T>
    void error_distribution(const T& x, int bins)
    {
      const int steps=30; // max number of bisection steps
      typedef typename T::ElementType NumberType;
      NumberType total_error = x.one_norm();
      NumberType total_elements = x.N();
      NumberType max_error = x.infinity_norm();
      std::vector<NumberType> left(bins,0.0);
      std::vector<NumberType> right(bins,max_error*(1.0+1e-8));
      std::vector<NumberType> eta(bins);
      std::vector<NumberType> target(bins);
      for (unsigned int k=0; k<bins; k++)
        target[k]= (k+1)/((NumberType)bins);
      for (int j=1; j<=steps; j++)
        {
          for (unsigned int k=0; k<bins; k++)
            eta[k]= 0.5*(left[k]+right[k]);
          std::vector<NumberType> sum(bins,0.0);
          std::vector<int> count(bins,0);

          for (typename T::const_iterator it = x.begin(),
                 end = x.end();
               it != end;
               ++it)
            {
              for (unsigned int k=0; k<bins; k++)
                if (*it<=eta[k])
                  {
                    sum[k] += *it;
                    count[k] += 1;
                  }
            }
          // std::cout << std::endl;
          // std::cout << "// step " << j << std::endl;
          // for (unsigned int k=0; k<bins; k++)
          //    std::cout << k+1 << " " << count[k] << " " << eta[k] << " " << right[k]-left[k]
          //          << " " << sum[k]/total_error << " " << target[k] << std::endl;
          for (unsigned int k=0; k<bins; k++)
            if (sum[k]<=target[k]*total_error)
              left[k] = eta[k];
            else
              right[k] = eta[k];
        }
      std::vector<NumberType> sum(bins,0.0);
      std::vector<int> count(bins,0);
      for (unsigned int i=0; i<x.N(); i++)
        for (unsigned int k=0; k<bins; k++)
          if (x[i]<=eta[k])
            {
              sum[k] += x[i];
              count[k] += 1;
            }
      std::cout << "+++ error distribution" << std::endl;
      std::cout << "+++ number of elements: " << x.N() << std::endl;
      std::cout << "+++ max element error:  " << max_error << std::endl;
      std::cout << "+++ total error:        " << total_error << std::endl;
      std::cout << "+++ bin #elements eta sum/total " << std::endl;
      for (unsigned int k=0; k<bins; k++)
        std::cout << "+++ " << k+1 << " " << count[k] << " " << eta[k] << " " << sum[k]/total_error << std::endl;
    }

    template<typename Grid, typename X>
    void mark_grid (Grid &grid, const X& x, typename X::ElementType refine_threshold,
                    typename X::ElementType coarsen_threshold, int min_level = 0, int max_level = std::numeric_limits<int>::max(), int verbose=0)
    {
      //typedef typename Grid::template Codim<0>::template Partition<Dune::All_Partition>::LeafIterator
      //  Iterator;
      typedef typename Grid::template Partition<Dune::All_Partition>::LeafGridView GV;
      typedef typename GV::template Codim<0>::Iterator Iterator;

      const GV& gv=grid.template leafGridView<Dune::All_Partition>();
      //Iterator it = grid.template leafbegin<0,Dune::All_Partition>();
      //Iterator eit = grid.template leafend<0,Dune::All_Partition>();
      Iterator it = gv.template begin<0>();
      Iterator eit = gv.template end<0>();

      unsigned int refine_cnt=0;
      unsigned int coarsen_cnt=0;

      typedef typename X::GridFunctionSpace GFS;
      typedef LocalFunctionSpace<GFS> LFS;
      typedef LFSIndexCache<LFS> LFSCache;
      typedef typename X::template ConstLocalView<LFSCache> XView;

      LFS lfs(x.gridFunctionSpace());
      LFSCache lfs_cache(lfs);
      XView x_view(x);

      for(;it!=eit;++it)
        {
          lfs.bind(*it);
          lfs_cache.update();
          x_view.bind(lfs_cache);

          if (x_view[0]>=refine_threshold && it->level() < max_level)
            {
              grid.mark(1,*(it));
              refine_cnt++;
            }
          if (x_view[0]<=coarsen_threshold && it->level() > min_level)
            {
              grid.mark(-1,*(it));
              coarsen_cnt++;
            }
          x_view.unbind();
        }
      if (verbose>0)
        std::cout << "+++ mark_grid: " << refine_cnt << " marked for refinement, "
                  << coarsen_cnt << " marked for coarsening" << std::endl;
    }


    template<typename Grid, typename X>
    void mark_grid_for_coarsening (Grid &grid, const X& x, typename X::ElementType refine_threshold,
                                   typename X::ElementType coarsen_threshold, int verbose=0)
    {
      typedef typename Grid::template Codim<0>::template Partition<Dune::All_Partition>::LeafIterator
        Iterator;
      typedef typename Grid::LeafGridView GV;
      typedef typename GV::IndexSet IndexSet;

      const GV& gv=grid.leafGridView();
      const IndexSet& is(gv.indexSet());
      Iterator it = grid.template leafbegin<0,Dune::All_Partition>();
      Iterator eit = grid.template leafend<0,Dune::All_Partition>();

      unsigned int coarsen_cnt=0;

      typedef typename X::GridFunctionSpace GFS;
      typedef LocalFunctionSpace<GFS> LFS;
      typedef LFSIndexCache<LFS> LFSCache;
      typedef typename X::template ConstLocalView<LFSCache> XView;

      LFS lfs(x.gridFunctionSpace());
      LFSCache lfs_cache(lfs);
      XView x_view(x);

      for(;it!=eit;++it)
        {
          lfs.bind(*it);
          lfs_cache.update();
          x_view.bind(lfs_cache);

          if (x_view[0]>=refine_threshold)
            {
              grid.mark(-1,*(it));
              coarsen_cnt++;
            }
          if (x_view[0]<=coarsen_threshold)
            {
              grid.mark(-1,*(it));
              coarsen_cnt++;
            }
        }
      if (verbose>0)
        std::cout << "+++ mark_grid_for_coarsening: "
                  << coarsen_cnt << " marked for coarsening" << std::endl;
    }


    class TimeAdaptationStrategy
    {
      // strategy parameters
      double scaling;
      double optimistic_factor;
      double coarsen_limit;
      double balance_limit;
      double tol;
      double T;
      int verbose;
      bool no_adapt;
      double refine_fraction_while_refinement;
      double coarsen_fraction_while_refinement;
      double coarsen_fraction_while_coarsening;
      double timestep_decrease_factor;
      double timestep_increase_factor;

      // results to be reported to the user after evaluating the error
      bool accept;
      bool adapt_dt;
      bool adapt_grid;
      double newdt;
      double q_s, q_t;

      // state variables
      bool have_decreased_time_step;
      bool have_refined_grid;

      // the only state variable: accumulated error
      double accumulated_estimated_error_squared;
      double minenergy_rate;

    public:
      TimeAdaptationStrategy (double tol_, double T_, int verbose_=0)
        : scaling(16.0), optimistic_factor(1.0), coarsen_limit(0.5), balance_limit(0.33333),
          tol(tol_), T(T_), verbose(verbose_), no_adapt(false),
          refine_fraction_while_refinement(0.7),
          coarsen_fraction_while_refinement(0.2),
          coarsen_fraction_while_coarsening(0.2),
          timestep_decrease_factor(0.5), timestep_increase_factor(1.5),
          accept(false), adapt_dt(false), adapt_grid(false), newdt(1.0),
          have_decreased_time_step(false), have_refined_grid(false),
          accumulated_estimated_error_squared(0.0),
          minenergy_rate(0.0)
      {
      }

      void setTimeStepDecreaseFactor (double s)
      {
        timestep_decrease_factor=s;
      }

      void setTimeStepIncreaseFactor (double s)
      {
        timestep_increase_factor=s;
      }

      void setRefineFractionWhileRefinement (double s)
      {
        refine_fraction_while_refinement=s;
      }

      void setCoarsenFractionWhileRefinement (double s)
      {
        coarsen_fraction_while_refinement=s;
      }

      void setCoarsenFractionWhileCoarsening (double s)
      {
        coarsen_fraction_while_coarsening=s;
      }

      void setMinEnergyRate (double s)
      {
        minenergy_rate=s;
      }

      void setCoarsenLimit (double s)
      {
        coarsen_limit=s;
      }

      void setBalanceLimit (double s)
      {
        balance_limit=s;
      }

      void setTemporalScaling (double s)
      {
        scaling=s;
      }

      void setOptimisticFactor (double s)
      {
        optimistic_factor=s;
      }

      void setAdaptationOn ()
      {
        no_adapt = false;
      }

      void setAdaptationOff ()
      {
        no_adapt = true;
      }

      bool acceptTimeStep () const
      {
        return accept;
      }

      bool adaptDT () const
      {
        return adapt_dt;
      }

      bool adaptGrid () const
      {
        return adapt_grid;
      }

      double newDT () const
      {
        return newdt;
      }

      double qs () const
      {
        return q_s;
      }

      double qt () const
      {
        return q_t;
      }

      double endT() const
      {
        return T;
      }

      double accumulatedErrorSquared () const
      {
        return accumulated_estimated_error_squared;
      }

      // to be called when new time step is done
      void startTimeStep ()
      {
        have_decreased_time_step = false;
        have_refined_grid = false;
      }

      template<typename GM, typename X>
      void evaluate_estimators (GM& grid, double time, double dt, const X& eta_space,  const X& eta_time, double energy_timeslab)
      {
        accept=false;
        adapt_dt=false;
        adapt_grid=false;
        newdt=dt;

        double spatial_error = eta_space.one_norm();
        double temporal_error = scaling*eta_time.one_norm();
        double sum = spatial_error + temporal_error;
        //double allowed = optimistic_factor*(tol*tol-accumulated_estimated_error_squared)*dt/(T-time);
        double allowed = tol*tol*(energy_timeslab+minenergy_rate*dt);
        q_s = spatial_error/sum;
        q_t = temporal_error/sum;

        // print some statistics
        if (verbose>0)
          std::cout << "+++"
                    << " q_s=" << q_s
                    << " q_t=" << q_t
                    << " sum/allowed=" << sum/allowed
                    // << " allowed=" << allowed
                    << " estimated error=" << sqrt(accumulated_estimated_error_squared+sum)
                    << " energy_rate=" << energy_timeslab/dt
                    << std::endl;

        // for simplicity: a mode that does no adaptation at all
        if (no_adapt)
          {
            accept = true;
            accumulated_estimated_error_squared += sum;
            if (verbose>1) std::cout << "+++ no adapt mode" << std::endl;
            return;
          }

        // the adaptation strategy
        if (sum<=allowed)
          {
            // we will accept this time step
            accept = true;
            if (verbose>1) std::cout << "+++ accepting time step" << std::endl;
            accumulated_estimated_error_squared += sum;

            // check if grid size or time step needs to be adapted
            if (sum<coarsen_limit*allowed)
              {
                // the error is too small, i.e. the computation is inefficient
                if (q_t<balance_limit)
                  {
                    // spatial error is dominating => increase time step
                    newdt = timestep_increase_factor*dt;
                    adapt_dt = true;
                    if (verbose>1) std::cout << "+++ spatial error dominates: increase time step" << std::endl;
                  }
                else
                  {
                    if (q_s>balance_limit)
                      {
                        // step sizes balanced: coarsen in time
                        newdt = timestep_increase_factor*dt;
                        adapt_dt = true;
                        if (verbose>1) std::cout << "+++ increasing time step" << std::endl;
                      }
                    // coarsen grid in space
                    double eta_refine, eta_coarsen;
                    if (verbose>1) std::cout << "+++ mark grid for coarsening" << std::endl;
                    //error_distribution(eta_space,20);
                    Dune::PDELab::error_fraction(eta_space,coarsen_fraction_while_coarsening,
                                                 coarsen_fraction_while_coarsening,eta_refine,eta_coarsen);
                    Dune::PDELab::mark_grid_for_coarsening(grid,eta_space,eta_refine,eta_coarsen,verbose);
                    adapt_grid = true;
                  }
              }
            else
              {
                // modify at least the time step
                if (q_t<balance_limit)
                  {
                    // spatial error is dominating => increase time step
                    newdt = timestep_increase_factor*dt;
                    adapt_dt = true;
                    if (verbose>1) std::cout << "+++ spatial error dominates: increase time step" << std::endl;
                  }
              }
          }
        else
          {
            // error is too large, we need to do something
            if (verbose>1) std::cout << "+++ will redo time step" << std::endl;
            if (q_t>1-balance_limit)
              {
                // temporal error is dominating => decrease time step only
                newdt = timestep_decrease_factor*dt;
                adapt_dt = true;
                have_decreased_time_step = true;
                if (verbose>1) std::cout << "+++ decreasing time step only" << std::endl;
              }
            else
              {
                if (q_t<balance_limit)
                  {
                    if (!have_decreased_time_step)
                      {
                        // time steps size not balanced (too small)
                        newdt = timestep_increase_factor*dt;
                        adapt_dt = true;
                        if (verbose>1) std::cout << "+++ increasing time step" << std::endl;
                      }
                  }
                else
                  {
                    // step sizes balanced: refine in time as well
                    newdt = timestep_decrease_factor*dt;
                    adapt_dt = true;
                    have_decreased_time_step = true;
                    if (verbose>1) std::cout << "+++ decreasing time step" << std::endl;
                  }
                // refine grid in space
                double eta_refine, eta_coarsen;
                if (verbose>1) std::cout << "+++ BINGO mark grid for refinement and coarsening" << std::endl;
                //error_distribution(eta_space,20);
                Dune::PDELab::error_fraction(eta_space,refine_fraction_while_refinement,
                                             coarsen_fraction_while_refinement,eta_refine,eta_coarsen,0);
                Dune::PDELab::mark_grid(grid,eta_space,eta_refine,eta_coarsen,verbose);
                adapt_grid = true;
              }
          }
      }
    };



  } // namespace PDELab
} // namespace Dune

#endif
