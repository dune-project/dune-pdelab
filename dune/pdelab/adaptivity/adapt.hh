// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_ADAPT_HH
#define DUNE_PDELAB_ADAPT_HH

#include<dune/common/exceptions.hh>

#include<limits>
#include<vector>
#include<map>
#include<dune/geometry/quadraturerules.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/localfunctionspace.hh>

// for CTLFA
#include<dune/pdelab/common/function.hh>
// for InterpolateBackendStandard
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
// for ErrorEstimator
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
// for intersectionoperator
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/flags.hh>

#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

namespace Dune {
  namespace PDELab {

    /*! @class CoeffsToLocalFunctionAdapter
     *
     *  @brief A wrapper class interpreting coefficients as a LocalFunction
     *
     *         A wrapper class interpreting a vector of coefficients (and a local finite element)
     *         as a LocalFunction living on an Elem. The function may be evaluated on another Elem,
     *         e.g. a child or grandchild. Useful if the GridFunctionSpace the coefficients originate from is
     *         no longer available (think adaptation).
     *
     *  @tparam CoeffType Type of the interpreted coefficients
     *  @tparam DGF       The DiscreteGridFunction this should mimic
     *  @tparam FEM       The FiniteElementMap used
     *  @tparam E         Type for Elems
     */
    template<class CoeffType, class DGF, class FEM, class E>
      class CoeffsToLocalFunctionAdapter : public FunctionInterface<typename DGF::Traits,
      CoeffsToLocalFunctionAdapter<CoeffType,DGF,FEM,E> >
    {

      public:

        typedef typename DGF::Traits Traits;
        typedef typename FEM::Traits::FiniteElementType FiniteElement;

        /** @brief The constructor.
         *
         *  @param[in] coeff_ A vector of coefficients
         *  @param[in] fem_ A FiniteElementMap used for interpretation of the coefficients
         *  @param[in] from_ The Elem on which the LocalFunction lives
         *  @param[in] to_ The Elem on which the LocalFunction is evaluated (should be e.g. a child
         *             of from_ to make the values of the shape functions meaningful)
         */
        CoeffsToLocalFunctionAdapter (const std::vector<CoeffType>& coeff_, const FEM& fem_, const E& from_, const E& to_)
          : coeff(coeff_), fem(fem_), from(from_), to(to_), fe(fem.find(from)) {}

        /** @brief Special version of the constructor with to = from.
         *
         *  @param[in] coeff_ A vector of coefficients
         *  @param[in] fem_ A FiniteElementMap used for interpretation of the coefficients
         *  @param[in] elem_ The Elem on which the LocalFunction lives
         */
        CoeffsToLocalFunctionAdapter (const std::vector<CoeffType>& coeff_, const FEM& fem_, const E& elem_)
          : coeff(coeff_), fem(fem_), from(elem_), to(elem_), fe(fem.find(from)) {}

        /** @brief Evaluate the LocalFunction at the given position
         *
         *  @param[in]  x The position in local coordinates
         *  @param[out] y The result of the evaluation
         */
        inline void evaluate (const typename Traits::DomainType& x, typename Traits::RangeType& y) const
        {
          std::vector<typename Traits::RangeType> yVector;
          fe.localBasis().evaluateFunction(from.geometry().local(to.geometry().global(x)),yVector);
          if (yVector.size() != coeff.size()) DUNE_THROW(Dune::Exception,"Coefficient vector has wrong length in CoeffsToLocalFunctionAdapter");
          typename Traits::RangeType sum=0.0;
          for (unsigned int i = 0; i < yVector.size(); ++i)
            sum += yVector[i]*coeff[i];
          y = sum;
        }

      private:
        const std::vector<CoeffType>& coeff;
        const FEM& fem;
        const E& from;
        const E& to;
        const FiniteElement& fe;
    };

    /*! @class GradientSmoothnessOperator
     *
     *  @brief An Operator evaluating the jumps of the gradient over the
     *         intersections between Elems.
     */
    class GradientSmoothnessOperator :
      public Dune::PDELab::NumericalJacobianSkeleton<GradientSmoothnessOperator>,
      public Dune::PDELab::LocalOperatorDefaultFlags
    {
      public:

        enum {doAlphaSkeleton = true};

        GradientSmoothnessOperator() {}

        /*! @brief The skeleton term of the operator.
         *
         * @param[in]  ig     The Geometry of the Intersection
         * @param[in]  lfsu_s The ansatz space on the "inside" Elem
         * @param[in]  x_s    The coefficients on the "inside" Elem
         * @param[in]  lfsv_s The trial space  on the "inside" Elem
         * @param[in]  lfsu_n The ansatz space on the "outside" Elem
         * @param[in]  x_n    The coefficients on the "outside" Elem
         * @param[in]  lfsv_n The trial space  on the "outside" Elem
         * @param[out] r_s    The residual on the "inside" Elem
         * @param[out] r_n    The residual on the "outside" Elem
         */
        template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
          void alpha_skeleton (const IG& ig,
              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
              R& r_s, R& r_n) const
          {
            typedef typename LFSU::Traits::FiniteElementType::
              Traits::LocalBasisType::Traits::DomainFieldType DF;
            typedef typename LFSU::Traits::FiniteElementType::
              Traits::LocalBasisType::Traits::RangeFieldType RF;
            typedef typename LFSU::Traits::FiniteElementType::
              Traits::LocalBasisType::Traits::JacobianType JacobianType;
            typedef typename LFSU::Traits::FiniteElementType::
              Traits::LocalBasisType::Traits::RangeType RangeType;
            typedef typename LFSU::Traits::SizeType size_type;

            const int intorder = 2;
            const int dim = IG::dimension;
            const int dimw = IG::dimensionworld;

            const Dune::QuadratureRule<DF,dim-1>&
              rule = Dune::QuadratureRules<DF,dim-1>::rule(ig.geometry().type(),intorder);

            RF integrand(0.);

            for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator
                it=rule.begin(); it!=rule.end(); ++it)
            {
              // position in inside elem
              Dune::FieldVector<DF,dim> inside_pos = ig.geometryInInside().global(it->position());

              // gradients of basis on reference elem
              std::vector<JacobianType> js_s(lfsu_s.size());
              lfsu_s.finiteElement().localBasis().evaluateJacobian(inside_pos,js_s);

              // gradients of basis on inside elem
              const Dune::FieldMatrix<DF,dim,dim>
                jac_s = ig.inside()->geometry().jacobianInverseTransposed(inside_pos);
              std::vector<Dune::FieldVector<RF,dim> > gradphi_s(lfsu_s.size());
              for (size_type i=0; i<lfsu_s.size(); i++)
                jac_s.mv(js_s[i][0],gradphi_s[i]);

              // compute gradient of u on inside elem
              Dune::FieldVector<RF,dim> gradu_s(0.0);
              for (size_type i=0; i<lfsu_s.size(); i++)
                gradu_s.axpy(x_s(lfsu_s,i),gradphi_s[i]);

              // position in outside elem
              Dune::FieldVector<DF,dim> outside_pos = ig.geometryInOutside().global(it->position());

              // gradient on reference elem
              std::vector<JacobianType> js_n(lfsu_n.size());
              lfsu_n.finiteElement().localBasis().evaluateJacobian(outside_pos,js_n);

              // gradient on outside elem
              const Dune::FieldMatrix<DF,dim,dim>
                jac_n = ig.outside()->geometry().jacobianInverseTransposed(outside_pos);
              std::vector<Dune::FieldVector<RF,dimw> > gradphi_n(lfsu_n.size());
              for (size_type i=0; i<lfsu_n.size(); i++)
                jac_n.mv(js_n[i][0],gradphi_n[i]);

              // compute gradient of u on outside elem
              Dune::FieldVector<RF,dim> gradu_n(0.0);
              for (size_type i=0; i<lfsu_n.size(); i++)
                gradu_n.axpy(x_n(lfsu_n,i),gradphi_n[i]);

              // jump of gradient
              const Dune::FieldVector<DF,dim> outer_normal = ig.unitOuterNormal(it->position());
              RF grad_normal(0.0);
              for (int i=0; i<dim; i++) grad_normal += (gradu_s[i]-gradu_n[i])*outer_normal[i];

              // integrate
              RF factor = it->weight()*ig.inside()->geometry().integrationElement(inside_pos);
              integrand += grad_normal * grad_normal * factor;
            }

            // compute estimate for diameter of intersection
            DF hmax = -1.0E00;
            if (dim==1)
              {
                hmax = 1.0;
              }
            else
              {
                for (int i=1; i<ig.geometry().corners(); i++)
                  {
                    Dune::FieldVector<DF,dim> x = ig.geometry().corner(0);
                    Dune::FieldVector<DF,dim> y = ig.geometry().corner(i);
                    x -= y;
                    hmax = std::max(hmax,x.two_norm());
                  }
              }
            r_s.accumulate(lfsu_s,0,sqrt(hmax*integrand));
            r_n.accumulate(lfsu_n,0,sqrt(hmax*integrand));
          }
    };

    /*! @class ResidualErrorEstimation
     *
     * @brief Estimation of the error
     *
     *        Estimation of the error by evaluation of the alpha_volume member of the LocalOperator
     *        (residual on the Elems) and e.g. the jumps of the gradient (residual on the Intersections)
     *
     * @tparam GFSU Type of GridFunctionSpace used
     * @tparam U    Container class for the solution
     * @tparam LOP  Type of the LocalOperator calculating the error estimate
     */
    template<class GFSU, class U, class LOP = Dune::PDELab::GradientSmoothnessOperator, bool nonoverlapping_mode = false>
      class ResidualErrorEstimation
      {
        typedef typename GFSU::Traits::GridViewType GV;
        typedef typename GV::Grid::template Codim<0>::Entity Element;
        typedef typename GV::template Codim<0>::template Partition<Dune::Interior_Partition>::Iterator LeafIterator;
        typedef typename GV::IntersectionIterator IntersectionIterator;
        typedef typename IntersectionIterator::Intersection Intersection;
        typedef typename GV::ctype Coord;
        typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,double,GV::dimension> P0FEM;
        typedef GridFunctionSpace<GV,P0FEM> GFSV;
        typedef typename Dune::PDELab::BackendVectorSelector<GFSV,double>::Type V;

        public:

        /*! @brief The constructor.
         *
         * @param[in] gfsu_ The GridFunctionSpace of the problem
         * @param[in] lop_  The LocalOperator of the problem
         * @param[in] gsop_ A LocalOperator to evaluate on the Intersections, defaults to GradientSmoothnessOperator
         */
        ResidualErrorEstimation(const GFSU& gfsu_, const LOP& lop_ = GradientSmoothnessOperator())
          : gfsu(gfsu_), lop(lop_) {}

        /*! @brief Calculate an estimate of the error.
         *
         * @param[in]  u        The current solution, tested against the operators
         * @param[out] estimate The error estimate, one value per Elem (i.e. P0)
         */
        void apply(const U& u, V& estimate)
        {
          //! @todo allgemein
          static Dune::GeometryType simplex(Dune::GeometryType::simplex,GV::dimension);
          P0FEM p0fem(simplex);
          GFSV gfsv(gfsu.gridView(),p0fem);

          // make local function spaces
          typedef LocalFunctionSpace<GFSU> LFSU;
          LFSU lfsu(gfsu);
          LFSU lfsu_n(gfsu);
          typedef LocalFunctionSpace<GFSV> LFSV;
          LFSV lfsv(gfsv);
          LFSV lfsv_n(gfsv);

          LocalVector<typename U::ElementType, TrialSpaceTag> ul;
          LocalVector<typename U::ElementType, TrialSpaceTag> ul_n;

          typedef LocalVector<typename V::ElementType, TestSpaceTag> LV;
          typedef typename LV::WeightedAccumulationView LVView;
          LV vl;
          LV vl_n;

          // traverse grid view
          for (LeafIterator it = gfsu.gridView().template begin<0,Dune::Interior_Partition>();
              it!=gfsu.gridView().template end<0,Dune::Interior_Partition>(); ++it)
          {
            const Element& e = *it;
            // bind local function spaces to element
            lfsu.bind(e);
            lfsv.bind(e);

            ul.resize(lfsu.size());
            vl.assign(lfsv.size(),0.0);

            // read coefficents
            lfsu.vread(u,ul);

            LVView vlview = vl.weightedAccumulationView(1.0);
            // volume evaluation
            LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolume>::
              alpha_volume(lop,ElementGeometry<Element>(e),lfsu,ul,lfsv,vlview);
            LocalAssemblerCallSwitch<LOP,LOP::doLambdaVolume>::
              lambda_volume(lop,ElementGeometry<Element>(e),lfsv,vlview);

            lfsv.vadd(vl,estimate);

            // skeleton term
            IntersectionIterator endit = gfsu.gridView().iend(*it);
            for (IntersectionIterator iit = gfsu.gridView().ibegin(*it); iit!=endit; ++iit)
            {
              if (iit->neighbor())
              {
                // bind local function spaces to neighbor
                lfsu_n.bind(*(iit->outside()));
                lfsv_n.bind(*(iit->outside()));

                ul_n.resize(lfsu_n.size());
                vl.assign(lfsv.size(),0.0);
                vl_n.assign(lfsv_n.size(),0.0);

                // only assemble  where we have a solution
                if (!nonoverlapping_mode || iit->outside()->partitionType()==Dune::InteriorEntity)
                {
                  // read coefficients
                  lfsu_n.vread(u,ul_n);

                  LVView vlview_n = vl_n.weightedAccumulationView(1.0);
                  LocalAssemblerCallSwitch<LOP,LOP::doAlphaSkeleton>::
                    alpha_skeleton(lop,IntersectionGeometry<Intersection>(*iit,0),
                                   lfsu,ul,lfsv,lfsu_n,ul_n,lfsv_n,vlview,vlview_n);
                }

                lfsv.vadd(vl,estimate);
                lfsv_n.vadd(vl_n,estimate);
              }
            } // end of intersection
          } // end of element
        }

        private:

        const GFSU& gfsu;
        const LOP& lop;
      };

    /*! @class AdaptationInterface
     *
     * @brief Interface base class for grid adaptation schemes
     *
     * @tparam Grid Type of the grid we want to adapt
     * @tparam U    Container class for the solution
     */
    template<class Grid, class U>
      class AdaptationInterface
      {
        typedef typename Grid::template Codim<0>::Entity Element;

        public:

        /*! @brief Prepare information before marking any of the Elements
         *
         * @param[in] u The current solution
         */
        inline virtual void prepare (const U& u) {}; // only if needed

        /*! @brief Estimate the error and mark Elems for refinement
         *
         * @param[in] e The Elem which is marked (or not)
         * @param[in] u The current solution
         */
        virtual void mark (const Element& e, const U& u) = 0;

        //! every abstract base class has a virtual destructor
        virtual ~AdaptationInterface () {}

      };

    /*! @class EstimationAdaptation
     *
     * @brief An adaptation scheme based on the estimation of the error.
     *
     * @tparam Grid       Type of the Grid we want to adapt
     * @tparam GFSU       Type of GridFunctionSpace used
     * @tparam U          Container class for the solution
     * @tparam Estimation Class that does the actual estimation of the error
     */
    template<class Grid, class GFSU, class U, class Estimation>
      class EstimationAdaptation : public AdaptationInterface<Grid,U>
    {
      typedef typename Grid::template Codim<0>::Entity Element;
      typedef typename Grid::LeafGridView LeafGridView;
      typedef typename LeafGridView::template Codim<0>::template Partition<Dune::Interior_Partition>::Iterator LeafIterator;
      typedef typename Grid::LeafIndexSet IndexSet;
      typedef typename IndexSet::IndexType IndexType;
      typedef typename GFSU::Traits::GridViewType GV;
      typedef typename GV::ctype Coord;
      typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,double,GV::dimension> P0FEM;
      typedef GridFunctionSpace<GV,P0FEM> GFSV;
      typedef LocalFunctionSpace<GFSV> LFSV;
      typedef typename Dune::PDELab::BackendVectorSelector<GFSV,double>::Type V;
      typedef std::multimap<typename V::ElementType, const IndexType> MapType;

      public:

      EstimationAdaptation(Grid& grid_, const GFSU& gfsu_, Estimation& estimation_,
          double refine_, double coarsen_ = 0., int min_ = 0, int max_ = std::numeric_limits<int>::max(), bool doCommunicate_ = true)
        : grid(grid_), gfsu(gfsu_), estimation(estimation_),refine(refine_), coarsen(coarsen_),
        min(min_), max(max_), doCommunicate(doCommunicate_), refinementMap(), localEstimate(0.) {}

      /*! @brief Prepare information before marking any of the Elements
       *
       * @param[in] u The current solution
       */
      void prepare (const U& u)
      {
        //! @todo allgemein
        static Dune::GeometryType simplex(Dune::GeometryType::simplex,GV::dimension);
        P0FEM p0fem(simplex);
        GFSV gfsv(gfsu.gridView(),p0fem);
        V estimate(gfsv,0.);
        estimation.apply(u,estimate);
        LFSV lfsv(gfsv);
        const LeafGridView leafView = grid.leafView();
        std::vector<typename V::ElementType> vl;
        const IndexSet& indexset = grid.leafIndexSet();
        refinementMap.clear();
        localEstimate = 0.;

        std::multimap<typename V::ElementType, const IndexType> tempMultiMap;
        for (LeafIterator it = leafView.template begin<0,Dune::Interior_Partition>();
            it!=leafView.template end<0,Dune::Interior_Partition>(); ++it)
        {
          const Element& e = *it;
          const IndexType& index = indexset.index(e);
          lfsv.bind(e);
          lfsv.vread(estimate,vl);
          tempMultiMap.insert(std::pair<typename V::ElementType,const IndexType>(vl[0],index));
          localEstimate += vl[0] * e.geometry().volume();
        }

        coarsenNumber = coarsen * leafView.size(0);
        refineNumber = (1. - refine) * leafView.size(0);

        unsigned int count = 0;
        for (typename std::multimap<typename V::ElementType, const IndexType>::const_iterator it = tempMultiMap.begin();
            it!=tempMultiMap.end(); ++it)
        {
          refinementMap.insert(std::pair<const IndexType, unsigned int>((*it).second,count++));
        }

        if (doCommunicate && grid.comm().size() > 1) communicate();
      }

      /*! @brief Estimate the error and mark Elems for refinement
       *
       * @param[in] e The Elem which is marked (or not)
       * @param[in] u The current solution
       */
      void mark (const Element& e, const U& u)
      {
        const int level = e.level();
        const IndexSet& indexset = grid.leafIndexSet();
        const IndexType index = indexset.index(e);
        const int i = (refinementMap.find(index))->second;

        if      (i > refineNumber && level < max)
        {
          grid.mark( 1, e);
        }
        else if (i < coarsenNumber && level > min)
        {
          grid.mark(-1, e);
        }
      }

      private:

      /*! @brief Communicate the number of Elems each processor wants to mark,
       * so they can shift the quota of refined Elems to where they
       * are the most useful.
       */
      void communicate ()
      {
        double globalEstimate = grid.comm().sum(localEstimate);
        //! @todo get rid of ghosts and overlap
        int    localNumElem = grid.leafView().size(0);
        int    globalNumElem  = grid.comm().sum(localNumElem);

        const double estimateQuotient = localEstimate / globalEstimate;
        refineNumber = localNumElem - estimateQuotient * refine * globalNumElem;
      }

      Grid& grid;
      const GFSU& gfsu;
      Estimation& estimation;
      const double refine, coarsen;
      int refineNumber, coarsenNumber;
      const int min, max;
      const bool doCommunicate;
      MapType refinementMap;
      double localEstimate;
    };

    /*! @class TestingAdaptation
     *
     * @brief Adaptation scheme @todo
     *
     * @tparam Grid       Type of the grid we want to adapt
     * @tparam GFSU       Type of ansatz space
     * @tparam U          Container class for the solution
     * @tparam Projection Projection used for testing
     */
    template<class Grid, class GFSU, class U, class Projection>
      class TestingAdaptation : public AdaptationInterface<Grid,U>
    {
      typedef typename Grid::template Codim<0>::Entity Element;
      typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
      typedef typename Grid::LeafGridView LeafGridView;
      typedef typename LeafGridView::template Codim<0>::template Partition<Dune::Interior_Partition>::Iterator LeafIterator;
      typedef typename Grid::LeafIndexSet IndexSet;
      typedef typename IndexSet::IndexType IndexType;
      typedef std::multimap<typename U::ElementType, const IndexType> MapType;

      typedef DiscreteGridFunction<GFSU, U> DGF;
      typedef typename GFSU::Traits::FiniteElementMapType FEM;
      typedef GridFunctionToLocalFunctionAdapter<DGF> GFTLFA;
      typedef CoeffsToLocalFunctionAdapter<typename U::ElementType,DGF,FEM,Element> CTLFA;

      public:

      /*! @brief The constructor.
       *
       * @todo Doc params
       */
      TestingAdaptation (Grid& grid_, const GFSU& gfsu_, Projection& projection_,
          double refine_, double coarsen_ = 0., int min_ = 0, int max_ = std::numeric_limits<int>::max())
        : grid(grid_), gfsu(gfsu_), projection(projection_), refine(refine_), coarsen(coarsen_),
        min(min_), max(max_), refinementMap() {}

      /*! @brief Prepare information before marking any of the Elements
       *
       * @param[in] u The current solution
       */
      void prepare (const U& u)
      {
        const LeafGridView leafView = grid.leafView();
        const IndexSet& indexset = grid.leafIndexSet();
        refinementMap.clear();

        std::multimap<typename U::ElementType, const IndexType> tempMultiMap;
        for (LeafIterator it = leafView.template begin<0,Dune::Interior_Partition>();
            it!=leafView.template end<0,Dune::Interior_Partition>(); ++it)
        {
          const Element& e = *it;
          const IndexType& index = indexset.index(e);
          tempMultiMap.insert(std::pair<typename U::ElementType,const IndexType>(test(e,u),index));
        }

        unsigned int count = 0;
        for (typename std::multimap<typename U::ElementType, const IndexType>::const_iterator it = tempMultiMap.begin(); it!=tempMultiMap.end(); ++it)
        {
          refinementMap.insert(std::pair<const IndexType, unsigned int>((*it).second,count++));
        }

        coarsenNumber = coarsen * leafView.size(0);
        refineNumber = (1. - refine) * leafView.size(0);
      }

      /*! @brief Estimate the error and mark Elems for refinement
       *
       * @param[in] e The Elem which is marked (or not)
       * @param[in] u The current solution
       */
      void mark (const Element& e, const U& u)
      {
        const int level = e.level();
        const IndexSet& indexset = grid.leafIndexSet();
        const IndexType index = indexset.index(e);
        const int i = (refinementMap.find(index))->second;

        if      (i > refineNumber && level < max)
        {
          grid.mark( 1, e);
        }
        else if (i < coarsenNumber && level > min)
        {
          grid.mark(-1, e);
        }
      }

      private:

      typename U::ElementType test (const Element& e, const U& u)
      {
        const FEM& fem = gfsu.finiteElementMap();
        DGF dgf(gfsu,u);

        ElementPointer pAncestor = e.father();
        while (e.geometry().volume() == pAncestor->geometry().volume())
        {
          if (pAncestor->level() == 0) DUNE_THROW(Dune::Exception,
              "search for ancestor of element with different volume (i.e. not a direct copy) walked over coarsest level.");
          pAncestor = pAncestor->father();
        }

        const Element& father = *pAncestor;
        const int localSize = (fem.find(father)).localBasis().size();
        std::vector<typename U::ElementType> uCoarse(localSize,0.);
        std::vector<typename U::ElementType> coarseBasis(localSize,0.);
        const typename Element::HierarchicIterator& hbegin = father.hbegin(grid.maxLevel());
        const typename Element::HierarchicIterator& hend   = father.hend(grid.maxLevel());

        for (typename Element::HierarchicIterator hit = hbegin; hit != hend; ++hit)
        {
          // only evaluate on entities with data
          if ((*hit).isLeaf())
          {
            typedef GridFunctionToLocalFunctionAdapter<DGF> GFTLFA;
            GFTLFA gftlfa(dgf,*hit);
            CTLFA  ctlfa(coarseBasis,fem,father,*hit);

            // iterate over canonical basis vectors
            for (unsigned int i = 0; i < coarseBasis.size(); ++i)
            {
              coarseBasis = projection.template inverseMassMatrix<CTLFA,FEM>(e,fem,i);
              projection.template apply<GFTLFA,CTLFA>(father,*hit,gftlfa,ctlfa,uCoarse[i]);
            }
          }
        }

        GFTLFA gftlfa(dgf,e);
        CTLFA  ctlfa(uCoarse,fem,father,e);

        typename U::ElementType norm = 0.;
        projection.template error<GFTLFA,CTLFA>(father,e,gftlfa,ctlfa,norm);

        return norm;
      }

      Grid& grid;
      const GFSU& gfsu;
      Projection& projection;
      const double refine, coarsen;
      int refineNumber, coarsenNumber;
      const int min, max;
      MapType refinementMap;
    };

    /*! @class CoarsenIfPossible
     *
     * @brief @todo
     *
     * @tparam Grid       Type of the grid we want to adapt
     * @tparam GFSU       Type of ansatz space
     * @tparam U          Container class for the solution
     * @tparam Projection Projection used for testing
     */
    template<class Grid, class GFSU, class U, class Projection>
      class CoarsenIfPossible : public AdaptationInterface<Grid,U>
    {
      typedef typename Grid::template Codim<0>::Entity Element;

      typedef DiscreteGridFunction<GFSU, U> DGF;
      typedef typename GFSU::Traits::FiniteElementMapType FEM;
      typedef GridFunctionToLocalFunctionAdapter<DGF> GFTLFA;
      typedef CoeffsToLocalFunctionAdapter<typename U::ElementType,DGF,FEM,Element> CTLFA;

      public:

      /*! @brief The constructor.
       *
       * @todo Doc params!
       */
      CoarsenIfPossible (Grid& grid_, const GFSU& gfsu_, Projection& projection_,
          double threshold_ = 0.1, int min_ = 0)
        : grid(grid_), gfsu(gfsu_), projection(projection_),
        threshold(threshold_), min(min_) {}

      /*! @brief Estimate the error and mark elems for refinement
       *
       * @todo Doc params!
       */
      void mark (const Element& e, const U& u)
      {
        const int level = e.level();

        if (level > min)
        {
          // based on residual this would be coarsened
          grid.mark(-1, e);

          //check if coarsening retains information
          const bool keep = test(e,u);
          if (keep) grid.mark(0, e);
        }
      }

      private:

      bool test (const Element& e, const U& u)
      {
        const FEM& fem = gfsu.finiteElementMap();
        DGF dgf(gfsu,u);

        const Element& father = *(e.father());
        const int localSize = (fem.find(father)).localBasis().size();
        std::vector<typename U::ElementType> uCoarse(localSize,0.);
        std::vector<typename U::ElementType> coarseBasis(localSize,0.);
        const typename Element::HierarchicIterator& hbegin = father.hbegin(grid.maxLevel());
        const typename Element::HierarchicIterator& hend   = father.hend(grid.maxLevel());

        for (typename Element::HierarchicIterator hit = hbegin; hit != hend; ++hit)
        {
          // only evaluate on entities with data
          if ((*hit).isLeaf())
          {
            typedef GridFunctionToLocalFunctionAdapter<DGF> GFTLFA;
            GFTLFA gftlfa(dgf,*hit);
            CTLFA  ctlfa(coarseBasis,fem,father,*hit);

            // iterate over canonical basis vectors
            for (int i = 0; i < coarseBasis.size(); ++i)
            {
              coarseBasis = projection.template inverseMassMatrix<CTLFA,FEM>(e,fem,i);
              projection.template apply<GFTLFA,CTLFA>(father,*hit,gftlfa,ctlfa,uCoarse[i]);
            }
          }
        }

        GFTLFA gftlfa(dgf,e);
        CTLFA  ctlfa(uCoarse,fem,father,e);

        typename U::ElementType norm = 0.;
        projection.template error<GFTLFA,CTLFA>(father,e,gftlfa,ctlfa,norm);

        return (norm > threshold);
      }

      Grid& grid;
      const GFSU& gfsu;
      Projection& projection;
      const double threshold;
      const int min;
    };

    /*! @class FunctionAdaptation
     *
     * @brief An adaptation scheme based on the value of a given function.
     *
     *        An adaptation scheme based on the value of a given function.
     *        Elems are marked in an effort to minimize the difference between
     *        their level and the average of the function on them.
     *
     * @tparam Grid Type of the grid we want to adapt
     * @tparam U    Container class for the solution
     * @tparam F    Type of Function used
     */
    template<class Grid, class U, class F>
      class FunctionAdaptation : public AdaptationInterface<Grid,U>
    {
      typedef typename Grid::template Codim<0>::Entity Element;

      public:

      /*! @brief The constructor.
       *
       * @param[in] grid_ The grid we want to adapt
       * @param[in] f_    The function used as guidance
       */
      FunctionAdaptation (Grid& grid_, F& f_) : grid(grid_), f(f_) {}

      /*! \brief Estimate the error and mark Elems for refinement
       *
       * @param[in] e The Elem which is marked (or not)
       * @param[in] u The current solution
       */
      void mark (const Element& e, const U& u)
      {
        const int level = e.level();

        //! @todo integrate over the Elem
        double diff = std::numeric_limits<double>::min();
        for (int i = 0; i < e.geometry().corners(); ++i)
        {
          //const typename F::Traits::RangeType& value = f.evaluate(e.geometry().center(),time);
          const double value = f.evaluate(e.geometry().corner(i));

          diff = std::max(diff,round(value) - level);
        }

        // only coarsen to level 0
        if (level > 0 && diff < 0) grid.mark(-1, e);
        if (diff > 0) grid.mark( 1, e);
      }

      private:

      Grid& grid;
      F& f;
    };

    /*! @class GlobalRefine
     *
     * @brief A pseudo adaptation scheme refining every element
     *
     * @tparam Grid Type of the grid we want to adapt
     * @tparam U    Container class for the solution
     */
    template<class Grid, class U>
      class GlobalRefine : public AdaptationInterface<Grid,U>
    {
      typedef typename Grid::template Codim<0>::Entity Element;

      public:

      /*! @brief The constructor.
       *
       * @todo Doc params
       */
      GlobalRefine (Grid& grid_, int max_ = std::numeric_limits<int>::max()) : grid(grid_), max(max_) {}

      /*! @brief Mark every Elem for refinement
       *
       * @param[in] e The Elem which is marked
       * @param[in] u The current solution
       */
      void mark (const Element& e, const U& u)
      {
        const int level = e.level();
        if (level < max) grid.mark(1, e);
      }

      private:

      Grid& grid;
      const int max;
    };

    /*! @class GlobalCoarsen
     *
     * @brief A pseudo adaptation scheme coarsening every element.
     *
     * @tparam Grid Type of the grid we want to adapt
     * @tparam U    Container class for the solution
     */
    template<class Grid, class U>
      class GlobalCoarsen : public AdaptationInterface<Grid,U>
    {
      typedef typename Grid::template Codim<0>::Entity Element;

      public:

      /*! @brief The constructor.
       *
       * @todo Doc params!
       */
      GlobalCoarsen (Grid& grid_, int min_ = 0) : grid(grid_), min(min_) {}

      /*! @brief Mark every Elem for coarsening
       *
       * @param[in] e The Elem which is marked
       * @param[in] u The current solution
       */
      void mark (const Element& e, const U& u)
      {
        const int level = e.level();
        if (level > min) grid.mark(-1, e);
      }

      private:

      Grid& grid;
      const int min;
    };

    /*! @class RandomAdaptation
     *
     * @brief A random grid adaptation scheme; @todo
     *
     * @tparam Grid Type of the grid we want to adapt
     * @tparam U    Container class for the solution
     */
    template<class Grid, class U>
      class RandomAdaptation : public AdaptationInterface<Grid,U>
    {
      typedef typename Grid::template Codim<0>::Entity Element;

      public:

      /*! @brief The constructor.
       *
       * @todo Doc params
       */
      RandomAdaptation (Grid& grid_, double refine_, double coarsen_ = 1., int min_ = 0, int max_ = std::numeric_limits<int>::max())
        : grid(grid_), refine(refine_), coarsen(coarsen_), min(min_), max(max_)
      {
        // initialize random numbers
        srand((unsigned)time(0));
      }

      /*! @brief Randomly mark Elems for refinement
       *
       * @param[in] e The Elem which is marked (or not)
       * @param[in] u The current solution
       */
      void mark (const Element& e, const U& u)
      {
        const double random = (double)rand()/RAND_MAX;
        const int level = e.level();

        if (random < refine && level < max)
        {
          grid.mark(1, e);
        }
        else if (random > 1 - coarsen && level > min)
        {
          if (e.level() > 0) grid.mark(-1, e);
        }
      }

      private:

      Grid& grid;
      const double refine, coarsen;
      const int min, max;
    };

    /*! @class L2Projection
     *
     * @brief @todo
     *
     * @tparam GFSU Type of ansatz space
     * @tparam U    Container class for the solution
     */
    template<class GFSU, class U>
      class L2Projection
      {
        typedef typename GFSU::Traits::GridViewType::Grid Grid;
        typedef typename Grid::template Codim<0>::Entity Element;
        typedef LocalFunctionSpace<GFSU> LFSU;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;

        public:

        /*! @brief The constructor.
         *
         * @todo Doc params!
         */
        explicit L2Projection(int intorder_=2) : intorder(intorder_), haveMatrix(), matrix() {}

        /*! @brief Calculate the L2 scalar product of functions X and Y on e, but in the geometry of its ancestor
         *
         * @todo Doc template params
         * @todo params
         */
        template<typename F1, typename F2>
          void apply (const Element& father, const Element& e, F1& X, F2& Y, typename U::ElementType& u) const
          {
            Dune::GeometryType gt = e.geometry().type();

            const int dim  = Element::Geometry::dimension;
            const QuadratureRule<DF,dim>& rule = QuadratureRules<DF,dim>::rule(gt,intorder);
            typename F1::Traits::RangeType x;
            typename F2::Traits::RangeType y;

            // iterate over quadrature points
            for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
            {
              X.evaluate(it->position(),x);
              Y.evaluate(it->position(),y);
              const typename U::ElementType factor = it->weight()
                * e.geometry().integrationElement(it->position())
                / father.geometry().integrationElement(father.geometry().local(e.geometry().global(it->position())));

              u += x * y * factor;
            }
          }

        /*! @brief Calculate the L2 norm of X - Y on e, but in the geometry of its ancestor
         *
         * @todo Doc template params
         * @todo Doc params
         */
        template<typename F1, typename F2>
          void error (const Element& father, const Element& e, F1& X, F2& Y, typename U::ElementType& u) const
          {
            Dune::GeometryType gt = e.geometry().type();

            const int dim  = Element::Geometry::dimension;
            const QuadratureRule<DF,dim>& rule = QuadratureRules<DF,dim>::rule(gt,intorder);
            typename F1::Traits::RangeType x;
            typename F2::Traits::RangeType y;

            typename U::ElementType temp = 0.;

            // iterate over quadrature points
            for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
            {
              X.evaluate(it->position(),x);
              Y.evaluate(it->position(),y);
              const typename U::ElementType factor = it->weight();
              //  * e.geometry().integrationElement(it->position());
              //  / father.geometry().integrationElement(father.geometry().local(e.geometry().global(it->position())));

              temp += (x-y) * (x-y) * factor;
            }

            u += sqrt(temp);
          }

        /*! @brief Calculate the inverse local mass matrix, used in the local L2 projection
         *
         * @todo Doc template params
         * @todo Doc params
         */
        template <class CTLFA, class FEM>
          const std::vector<double>& inverseMassMatrix(const Element& e, const FEM& fem, int k)
          {
            const Dune::GeometryType gt = e.geometry().type();
            if (haveMatrix.find(gt) != haveMatrix.end())
            {
              return matrix[gt][k];
            }

            const int localSize = (fem.find(e)).localBasis().size();
            typename CTLFA::Traits::RangeType x;
            typename CTLFA::Traits::RangeType y;

            std::vector<std::vector<typename U::ElementType> > massMatrix(localSize,std::vector<typename U::ElementType>(localSize,0.));
            std::vector<std::vector<typename U::ElementType> > inverseMatrix(localSize,std::vector<typename U::ElementType>(localSize,0.));

            std::vector<double> phi_i(localSize,0.), phi_j(localSize,0.);
            CTLFA ctlfa_i(phi_i,fem,e,e);
            CTLFA ctlfa_j(phi_j,fem,e,e);

            for (int i = 0; i < localSize; ++i)
            {
              ++phi_i[i];
              for (int j = 0; j < localSize; ++j)
              {
                ++phi_j[j];
                (*this).template apply<CTLFA,CTLFA>(e,e,ctlfa_i,ctlfa_j,massMatrix[i][j]);
                --phi_j[j];
              }
              --phi_i[i];
            }

            for (int i = 0; i < localSize; ++i)
            {
              inverseMatrix[i][i] = 1.;
            }

            for (int i = 0; i < localSize; ++i)
            {
              for (int j = 0; j < localSize; ++j)
              {
                if (i != j)
                {
                  const typename U::ElementType factor = massMatrix[j][i]/massMatrix[i][i];
                  for (int l = 0; l < localSize; ++l)
                  {
                    massMatrix[j][l] -= factor * massMatrix[i][l];
                    inverseMatrix[j][l] -= factor * inverseMatrix[i][l];
                  }
                }
              }
            }
            for (int i = localSize - 1; i >= 0; --i)
            {
              for (int j = localSize - 1; j >= 0; --j)
              {
                if (i != j)
                {
                  const typename U::ElementType factor = massMatrix[j][i]/massMatrix[i][i];
                  for (int l = localSize - 1; l >= 0; --l)
                  {
                    massMatrix[j][l] -= factor * massMatrix[i][l];
                    inverseMatrix[j][l] -= factor * inverseMatrix[i][l];
                  }
                }
              }
            }
            for (int i = 0; i < localSize; ++i)
            {
              const typename U::ElementType factor = massMatrix[i][i];
              for (int j = 0; j < localSize; ++j)
              {
                massMatrix[i][j] /= factor;
                inverseMatrix[i][j] /= factor;
              }
            }

            matrix.insert(std::pair<Dune::GeometryType,std::vector<std::vector<typename U::ElementType> > >(gt,inverseMatrix));
            haveMatrix.insert(gt);
            return matrix[gt][k];
          }

        private:

        const int intorder;

        // only for inverseMassMatrix
        std::set<Dune::GeometryType> haveMatrix;
        std::map<Dune::GeometryType, std::vector<std::vector<typename U::ElementType> > > matrix;
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
     * @tparam Adaptation Class that decides what to mark
     * @tparam Projection Projection used when Elems vanish
     */
    template<class Grid, class GFSU, class U,
      class Adaptation, class Projection>
        class GridAdaptor
        {
          typedef typename Grid::LeafGridView LeafGridView;
          typedef typename LeafGridView::template Codim<0>
            ::template Partition<Dune::Interior_Partition>::Iterator LeafIterator;
          typedef typename Grid::template Codim<0>::Entity Element;
          typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
          typedef typename Grid::GlobalIdSet IdSet;
          typedef typename IdSet::IdType IdType;
          typedef typename Grid::LeafIndexSet IndexSet;
          typedef typename IndexSet::IndexType IndexType;
          typedef LocalFunctionSpace<GFSU> LFSU;
          typedef DiscreteGridFunction<GFSU, U> DGF;
          typedef typename GFSU::Traits::FiniteElementMapType FEM;
          typedef InterpolateBackendStandard IB;
          typedef CoeffsToLocalFunctionAdapter<typename U::ElementType,DGF,FEM,Element> CTLFA;
          typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType RangeType;

          //typedef std::map<
          //                 int,
          //                 std::map<
          // 	                        IdType,
          //                          std::vector<typename U::ElementType>
          //		               >
          //		      > MapType;
          typedef std::map<IdType,std::vector<typename U::ElementType> > MapType;

          public:

          /*! @brief The constructor.
           *
           * @param grid_       The grid we want to adapt
           * @param gfsu_       The ansatz space, we need to update it
           * @param adaptation_ The Adaptation object deciding what to refine/coarsen
           * @param projection_ The Projection used when Elems vanish
           */
          GridAdaptor(Grid& grid_, GFSU& gfsu_,
              Adaptation& adaptation_, Projection& projection_)
            : grid(grid_), gfsu(gfsu_),
            adaptation(adaptation_), projection(projection_) {}

          /* @brief Complete one cycle of adaptation, i.e.
           *        @arg mark the Elems on the grid
           *        @arg get a backup of the solution
           *        @arg adapt the grid
           *        @arg update the function space
           *        @arg move the solution back onto the grid
           *
           * @param[in] u The current solution
           */
          void adapt (U& u)
          {
            // mark elems that are to be refined / coarsened
            markGrid(u);

            // prepare the grid for refinement
            grid.preAdapt();

            // save u
            MapType transferMap;
            backupData(u,transferMap);

            // adapt the grid
            grid.adapt();

            // update the function spaces
            gfsu.update();

            // reset u
            u = U(gfsu,0.0);
            replayData(u,transferMap);

            // clean up
            grid.postAdapt();
          }

          private:

          /* @brief @todo
           *
           * @param[in] u The current solution
           */
          void markGrid(const U& u)
          {
            adaptation.prepare(u);

            // iterate over all elems
            LeafGridView leafView = grid.leafView();
            for (LeafIterator it = leafView.template begin<0,Dune::Interior_Partition>();
                it!=leafView.template end<0,Dune::Interior_Partition>(); ++it)
            {
              adaptation.mark(*it,u);
            }
          }

          /* @brief @todo
           *
           * @param[in]  u           The solution that will be saved
           * @param[out] transferMap The map containing the solution during adaptation
           */
          void backupData(U& u, MapType& transferMap)
          {
            const IdSet& idset = grid.globalIdSet();
            LFSU lfsu(gfsu);
            DGF dgf(gfsu,u);
            const FEM& fem = gfsu.finiteElementMap();
            IB ib = IB();
            std::vector<typename U::ElementType> ul;

            // iterate over all elems
            LeafGridView leafView = grid.leafView();
            for (LeafIterator it = leafView.template begin<0,Dune::Interior_Partition>();
                it!=leafView.template end<0,Dune::Interior_Partition>(); ++it)
            {
              const Element& e = *it;

              // save local coeffs in map
              lfsu.bind(e);
              //lfsu.vread(u,transferMap[e.level()][idset.id(e)]);
              lfsu.vread(u,transferMap[idset.id(e)]);

              // save local coeffs of father in map
              if (e.mightVanish())
              {
                ElementPointer father = e.father();

                projection.template inverseMassMatrix<CTLFA,FEM>(e,fem,0);
                {
                  const int localSize = (fem.find(*father)).localBasis().size();
                  std::vector<typename U::ElementType> uCoarse(localSize,0.);
                  std::vector<typename U::ElementType> coarseBasis(localSize,0.);
                  const typename Element::HierarchicIterator& hbegin = (*father).hbegin(grid.maxLevel());
                  const typename Element::HierarchicIterator& hend   = (*father).hend(grid.maxLevel());

                  for (typename Element::HierarchicIterator hit = hbegin; hit != hend; ++hit)
                  {
                    // only evaluate on entities with data
                    if ((*hit).isLeaf())
                    {
                      typedef GridFunctionToLocalFunctionAdapter<DGF> GFTLFA;
                      GFTLFA gftlfa(dgf,*hit);
                      CTLFA  ctlfa(coarseBasis,fem,*father,*hit);

                      // iterate over canonical basis vectors
                      for (unsigned int i = 0; i < coarseBasis.size(); ++i)
                      {
                        coarseBasis = projection.template inverseMassMatrix<CTLFA,FEM>(e,fem,i);
                        projection.template apply<GFTLFA,CTLFA>(*father,*hit,gftlfa,ctlfa,uCoarse[i]);
                      }
                    }
                  }
                  //transferMap[(*father).level()][idset.id(*father)] = uCoarse;
                  transferMap[idset.id(*father)] = uCoarse;
                }

                // conforming grids may need this
                while ((*father).mightVanish())
                {
                  father = (*father).father();
                  {
                    const int localSize = (fem.find(*father)).localBasis().size();
                    std::vector<typename U::ElementType> uCoarse(localSize,0.);
                    std::vector<typename U::ElementType> coarseBasis(localSize,0.);
                    const typename Element::HierarchicIterator& hbegin = (*father).hbegin(grid.maxLevel());
                    const typename Element::HierarchicIterator& hend   = (*father).hend(grid.maxLevel());

                    for (typename Element::HierarchicIterator hit = hbegin; hit != hend; ++hit)
                    {
                      // only evaluate on entities with data
                      if ((*hit).isLeaf())
                      {
                        typedef GridFunctionToLocalFunctionAdapter<DGF> GFTLFA;
                        GFTLFA gftlfa(dgf,*hit);
                        CTLFA  ctlfa(coarseBasis,fem,*father,*hit);

                        // iterate over canonical basis vectors
                        for (unsigned int i = 0; i < coarseBasis.size(); ++i)
                        {
                          coarseBasis = projection.template inverseMassMatrix<CTLFA,FEM>(e,fem,i);
                          projection.template apply<GFTLFA,CTLFA>(*father,*hit,gftlfa,ctlfa,uCoarse[i]);
                        }
                      }
                    }
                    //transferMap[(*father).level()][idset.id(*father)] = uCoarse;
                    transferMap[idset.id(*father)] = uCoarse;
                  }
                }
              }
            }
          }

          /* @brief @todo
           *
           * @param[out] u           The solution after adaptation
           * @param[in]  transferMap The map that contains the information for the rebuild of u
           */
          void replayData(U& u, MapType& transferMap)
          {
            const IdSet& idset = grid.globalIdSet();
            LFSU lfsu(gfsu);
            const FEM& fem = gfsu.finiteElementMap();
            IB ib = IB();
            std::vector<typename U::ElementType> ul;
            std::vector<typename U::ElementType> ulc;
            U uc(gfsu,0.0);
            const IndexSet& indexset = grid.leafIndexSet();
            std::vector<typename U::ElementType> ug(indexset.size(IndexSet::dimension),0.);
            std::vector<typename U::ElementType> ugc(indexset.size(IndexSet::dimension),0.);
            Dune::LeafMultipleCodimMultipleGeomTypeMapper<Grid,Dune::MCMGElementLayout> mapper(grid);

            // iterate over all elems
            LeafGridView leafView = grid.leafView();
            for (LeafIterator it = leafView.template begin<0,Dune::Interior_Partition>();
                it!=leafView.template end<0,Dune::Interior_Partition>(); ++it)
            {
              const Element& e = *it;
              lfsu.bind(e);
              const IdType& id = idset.id(e);
              const int level = e.level();

              if (e.isNew()) // id is not in map, we have to interpolate
              {
                // find ancestor with data
                int levelAncestor = level - 1;
                ElementPointer pAncestor = e.father();
                //while (transferMap[levelAncestor].find(idset.id(*pAncestor)) == transferMap[levelAncestor].end())
                while (transferMap.find(idset.id(*pAncestor)) == transferMap.end())
                {
                  //this ancestor does not have data, check next one
                  if (levelAncestor == 0)
                    DUNE_THROW(Dune::Exception,
                        "transferMap of GridAdaptor didn't contain ancestor of element with id " << id);
                  pAncestor = (*pAncestor).father();
                  --levelAncestor;
                }

                // this one has the data
                const Element& Ancestor = *pAncestor;
                //CTLFA ctlfa(transferMap[levelAncestor][idset.id(Ancestor)],fem,Ancestor,e);
                CTLFA ctlfa(transferMap[idset.id(Ancestor)],fem,Ancestor,e);

                ul.clear();
                ib.interpolate(fem.find(e),ctlfa,ul);

                lfsu.vadd(ul,u);
              }
              else // this entity is not new and should have data
              {
                //lfsu.vadd(transferMap[level][id],u);
                lfsu.vadd(transferMap[id],u);
              }

              ulc = std::vector<typename U::ElementType>(lfsu.size(),1.0);
              lfsu.vadd(ulc,uc);
            }

            typedef Dune::PDELab::AddDataHandle<GFSU,U> Handle;
            Handle addHandle1(gfsu,u);
            leafView.communicate (addHandle1,
                Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
            Handle addHandle2(gfsu,uc);
            leafView.communicate (addHandle2,
                Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);

            for (LeafIterator it = leafView.template begin<0,Dune::Interior_Partition>();
                it!=leafView.template end<0,Dune::Interior_Partition>(); ++it)
            {
              const Element& e = *it;
              lfsu.bind(e);

              ul = std::vector<typename U::ElementType>(lfsu.size(),0.0);
              ulc = std::vector<typename U::ElementType>(lfsu.size(),0.0);
              lfsu.vread(u,ul);
              lfsu.vread(uc,ulc);

              for (unsigned int i = 0; i<ul.size();++i)
              {
                if(ulc[i]>1)
                {
                  ul[i] /= ulc[i];
                  ulc[i] = 1.;
                }
              }
              lfsu.vwrite(ul,u);
              lfsu.vwrite(ulc,uc);
            }
          }

          Grid&  grid;
          GFSU&  gfsu;
          Adaptation& adaptation;
          Projection& projection;
        };

  } // namespace PDELab
} // namespace Dune

#endif
