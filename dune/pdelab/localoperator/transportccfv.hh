// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_TRANSPORTCCFV_HH
#define DUNE_PDELAB_TRANSPORTCCFV_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/common/genericreferenceelements.hh>
#include<dune/localfunctions/raviartthomas/raviartthomas0q.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>

#include"../common/geometrywrapper.hh"
#include"../common/function.hh"
#include"../gridoperatorspace/gridoperatorspace.hh"
#include"pattern.hh"
#include"flags.hh"
#include"idefault.hh"

namespace Dune {
  namespace PDELab {

	//! traits class for two phase parameter class
	template<typename GV, typename RF>
	struct TransportParameterTraits
	{
	  //! \brief the grid view
	  typedef GV GridViewType;

	  //! \brief Enum for domain dimension
	  enum { 
		//! \brief dimension of the domain
		dimDomain = GV::dimension
	  }; 

	  //! \brief Export type for domain field
	  typedef typename GV::Grid::ctype DomainFieldType;

	  //! \brief domain type
	  typedef Dune::FieldVector<DomainFieldType,dimDomain> DomainType;

	  //! \brief domain type
	  typedef Dune::FieldVector<DomainFieldType,dimDomain-1> IntersectionDomainType;

	  //! \brief Export type for range field
	  typedef RF RangeFieldType;

	  //! \brief range type
	  typedef Dune::FieldVector<RF,GV::dimensionworld> RangeType;

	  //! grid types
	  typedef typename GV::Traits::template Codim<0>::Entity ElementType;
	  typedef typename GV::Intersection IntersectionType;
	};

    //! base class for parameter class
	template<class T, class Imp>
	class TransportSpatialParameterInterface
	{
	public:
	  typedef T Traits;

	  //! velocity field 
	  typename Traits::RangeType
	  v (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return asImp().v(e,x);
	  }

	  //! scalar diffusion coefficient
	  typename Traits::RangeFieldType 
	  D (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return asImp().D(e,x);
	  }

	  //! source term
	  typename Traits::RangeFieldType 
	  q (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return asImp().q(e,x);
	  }

	  //! boundary condition type function
      /** 
       * 0 means Neumann
       * 1 means Dirichlet
       * 2 means Outflow (zero diffusive flux, velocity points outward)
       */
	  int
	  bc (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
	  {
		return asImp().bc(is,x);
	  }

	  //! Dirichlet boundary condition on inflow
	  typename Traits::RangeFieldType 
	  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return asImp().g(e,x);
	  }

	  //! Neumann boundary condition
	  typename Traits::RangeFieldType 
	  j (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return asImp().j(e,x);
	  }

	private:
	  Imp& asImp () {return static_cast<Imp &> (*this);}
	  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
	};


    /*! Adapter that extracts boundary condition type function from parameter class

      \tparam T  model of ConvectionDiffusionParameterInterface
     */
    template<typename T>
    class BoundaryConditionType_Transport
      : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::
                                                      BoundaryGridFunctionTraits<typename T::Traits::GridViewType,int,1,
                                                                                 Dune::FieldVector<int,1> >,
                                                      BoundaryConditionType_Transport<T> >
    {
      const typename T::Traits::GridViewType& gv;
      const T& t;

    public:
      typedef Dune::PDELab::BoundaryGridFunctionTraits<typename T::Traits::GridViewType,int,1,
                                                       Dune::FieldVector<int,1> > Traits;
      typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,BoundaryConditionType_Transport<T> > BaseT;

      BoundaryConditionType_Transport (const typename T::Traits::GridViewType& gv_, const T& t_) : gv(gv_), t(t_) {}

      template<typename I>
      inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig, 
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {  
        y = t.bc(ig.intersection(),x);
      }

      //! get a reference to the GridView
      inline const typename T::Traits::GridViewType& getGridView ()
      {
        return gv;
      }
    };


    /*! Adapter that extracts Dirichlet boundary conditions from parameter class

      \tparam T  model of ConvectionDiffusionParameterInterface
     */
	template<typename T>
	class DirichletBoundaryCondition_Transport 
       : public GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                                                  typename T::Traits::RangeFieldType,
                                                                  1,Dune::FieldVector<typename T::Traits::RangeFieldType,1> >
                                                                  ,DirichletBoundaryCondition_Transport<T> >
	{
	public:
	  typedef Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                               typename T::Traits::RangeFieldType,
                                               1,Dune::FieldVector<typename T::Traits::RangeFieldType,1> > Traits;

      //! constructor 
	  DirichletBoundaryCondition_Transport (const typename Traits::GridViewType& g_, const T& t_) : g(g_), t(t_) {}

      //! \copydoc GridFunctionBase::evaluate()
	  inline void evaluate (const typename Traits::ElementType& e, 
							const typename Traits::DomainType& x, 
							typename Traits::RangeType& y) const
	  {  
        y = t.g(e,x);
	  }

	  inline const typename Traits::GridViewType& getGridView () const
	  {
		return g;
	  }
  
	private:
	  const typename Traits::GridViewType& g;
      const T& t;
	};

	/** a local operator for a cell-centered finite folume scheme for
        the transport equation
	  
		\nabla \cdot \{v u - D \nabla u \} = q in \Omega
                                         u = g on \Gamma_D
           \{v u - D \nabla u \} \cdot \nu = j on \Gamma_N
                                       outflow on \Gamma_O        

        Can be used for stationary and time-dependent computations

		\tparam TP  parameter class implementing ComponentTransportParameterInterface
	*/
    template<typename TP>
	class CCFVSpatialTransportOperator : 
      public NumericalJacobianApplySkeleton<CCFVSpatialTransportOperator<TP> >,
      public NumericalJacobianSkeleton<CCFVSpatialTransportOperator<TP> >,
      public NumericalJacobianApplyBoundary<CCFVSpatialTransportOperator<TP> >,
      public NumericalJacobianBoundary<CCFVSpatialTransportOperator<TP> >,
      public NumericalJacobianApplyVolumePostSkeleton<CCFVSpatialTransportOperator<TP> >,
      public NumericalJacobianVolumePostSkeleton<CCFVSpatialTransportOperator<TP> >,
      public FullSkeletonPattern,
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
	{
      enum { dim = TP::Traits::GridViewType::dimension };

	public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

	  // residual assembly flags
      enum { doAlphaVolume  = true };
      enum { doAlphaSkeleton  = true };
	  enum { doAlphaVolumePostSkeleton = true };
      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume    = true };

      enum { doSkeletonTwoSided = true }; // need to see face from both sides for CFL calculation

      CCFVSpatialTransportOperator (TP& tp_) 
		: tp(tp_)
	  {
	  }

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
        cellinflux = 0.0; // prepare dt computation
	  }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
                            LocalMatrix<R>& mat) const
      {
        // do nothing; alpha_volume only needed for dt computations
      }

	  // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_skeleton (const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                           R& r_s, R& r_n) const
	  {
		// domain and range field type
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
		typedef typename LFSU::Traits::GridFunctionSpaceType::Traits::BackendType B;

        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>& 
          face_local = Dune::GenericReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();

        // face center in element coordinates
        Dune::FieldVector<DF,IG::dimension> face_center_in_element = ig.geometryInInside().global(face_local);

        // evaluate velocity
        typename TP::Traits::RangeType v(tp.v(*(ig.inside()),face_center_in_element));

        // the normal velocity
        RF vn = v*ig.centerUnitOuterNormal();

        // convective flux
        RF u_upwind=0;
        if (vn>=0) u_upwind = x_s[0]; else u_upwind = x_n[0];
        r_s[0] += (u_upwind*vn)*face_volume;
        if (vn>=0)
          cellinflux += vn*face_volume; // dt computation

        // evaluate diffusion coefficients
        const Dune::FieldVector<DF,IG::dimension>&
          inside_local = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,IG::dimension>& 
          outside_local = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.outside()->type()).position(0,0);
        typename TP::Traits::RangeFieldType D_inside = tp.D(*(ig.inside()),inside_local);
        typename TP::Traits::RangeFieldType D_outside = tp.D(*(ig.outside()),outside_local);
        typename TP::Traits::RangeFieldType D_avg = 2.0/(1.0/(D_inside+1E-30) + 1.0/(D_outside+1E-30));

        // distance between cell centers in global coordinates
        Dune::FieldVector<DF,IG::dimension> 
          inside_global = ig.inside()->geometry().center();
        Dune::FieldVector<DF,IG::dimension> 
          outside_global = ig.outside()->geometry().center();
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();
 
        // diffusive flux
        // note: we do only one-sided evaluation here
        r_s[0] -= (D_avg*(x_n[0]-x_s[0])/distance)*face_volume;
      }

      // post skeleton: compute time step allowable for cell; to be done later
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume_post_skeleton(const EG& eg, const LFSU& lfsu, const X& x,
                                      const LFSV& lfsv, R& r) const
      {
		// domain and range field type
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
        const int dim = EG::Geometry::dimension;

        if (!first_stage) return; // time step calculation is only done in first stage

        // cell center
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);
        
        // compute optimal dt for this cell
        typename TP::Traits::RangeFieldType cellcapacity = tp.c(eg.entity(),inside_local)*eg.geometry().volume();
        typename TP::Traits::RangeFieldType celldt = cellcapacity/(cellinflux+1E-30);
        dtmin = std::min(dtmin,celldt);
      }

	  // skeleton integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_boundary (const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
	  {
		// domain and range field type
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
		typedef typename LFSU::Traits::GridFunctionSpaceType::Traits::BackendType B;
    
        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>& 
          face_local = Dune::GenericReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();
        Dune::FieldVector<DF,dim> face_center_in_element = ig.geometryInInside().global(face_local);

        // evaluate boundary condition type
        int bc = tp.bc(ig.intersection(),face_local);

        // do things depending on boundary condition type
        if (bc==0) // Neumann boundary
          {
            typename TP::Traits::RangeFieldType j = tp.j(*(ig.inside()),face_center_in_element);
            r_s[0] += j*face_volume;
            return;
         }

        // evaluate velocity
        typename TP::Traits::RangeType v(tp.v(*(ig.inside()),face_center_in_element));

        // the normal velocity
        RF vn = v*ig.centerUnitOuterNormal();
        if (vn>=0)
          cellinflux += vn*face_volume; // dt computation

        if (bc==2) // Outflow boundary
          {
            r_s[0] += vn*x_s[0]*face_volume;
            return;
          }

        if (bc==1) // Dirichlet boundary
          {
            typename TP::Traits::RangeFieldType g;
            if (vn>=0) g=x_s[0]; else g=tp.g(*(ig.inside()),face_center_in_element);
            const Dune::FieldVector<DF,IG::dimension>&
              inside_local = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
            typename TP::Traits::RangeFieldType D_inside = tp.D(*(ig.inside()),inside_local);
            Dune::FieldVector<DF,IG::dimension> 
              inside_global = ig.inside()->geometry().center();
            Dune::FieldVector<DF,IG::dimension> 
              outside_global = ig.geometry().center();
            inside_global -= outside_global;
            RF distance = inside_global.two_norm();
            r_s[0] += (g*vn - D_inside*(g-x_s[0])/distance)*face_volume;
            return;
          }
      }

 	  // volume integral depending only on test functions
	  template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
		// domain and range field type
		typedef typename LFSV::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSV::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
        const int dim = EG::Geometry::dimension;

        // cell center
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // evaluate source term
        typename TP::Traits::RangeFieldType q = tp.q(eg.entity(),inside_local);

        r[0] -= q*eg.geometry().volume();
      }

      //! set time in parameter class
      void setTime (typename TP::Traits::RangeFieldType t)
      {
      }

      //! to be called once before each time step
      void preStep (typename TP::Traits::RangeFieldType time, typename TP::Traits::RangeFieldType dt,
                    int stages)
      {
      }
      
      //! to be called once before each stage
      void preStage (typename TP::Traits::RangeFieldType time, int r)
      {
        if (r==1)
          {
            first_stage = true;
            dtmin = 1E100;
          }
        else first_stage = false;
      }

      //! to be called once at the end of each stage
      void postStage ()
      {
      }

      //! to be called once before each stage
      typename TP::Traits::RangeFieldType suggestTimestep (typename TP::Traits::RangeFieldType dt) const
      {
        return dtmin;
      }
      
	private:
	  TP& tp;
      bool first_stage;
      mutable typename TP::Traits::RangeFieldType dtmin; // accumulate minimum dt here
      mutable typename TP::Traits::RangeFieldType cellinflux;
	};


    //! base class for parameter class
	template<class T, class Imp>
	class TransportTemporalParameterInterface
	{
	public:
	  typedef T Traits;

	  //! source term
	  typename Traits::RangeFieldType 
	  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return asImp().c(e,x);
	  }

	private:
	  Imp& asImp () {return static_cast<Imp &> (*this);}
	  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
	};

    /** a local operator for the storage operator
     *
     * \f{align*}{
     \int_\Omega c(x) uv dx
     * \f}
     */
    template<class TP>
	class CCFVTemporalOperator : public NumericalJacobianApplyVolume<CCFVTemporalOperator<TP> >,
                                 public FullVolumePattern,
                                 public LocalOperatorDefaultFlags,
                                 public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
	{
	public:
      // pattern assembly flags
      enum { doPatternVolume = true };

	  // residual assembly flags
      enum { doAlphaVolume = true };

      CCFVTemporalOperator (TP& tp_) 
		: tp(tp_)
	  {
	  }

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
		// domain and range field type
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // cell center
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);
        
        // evaluate capacity
        typename TP::Traits::RangeFieldType c = tp.c(eg.entity(),inside_local);
        
        // residual contribution
        r[0] += c*x[0]*eg.geometry().volume();
	  }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
                            LocalMatrix<R>& mat) const
      {
		// domain and range field type
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType JacobianType;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // cell center
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);
        
        // evaluate capacity
        typename TP::Traits::RangeFieldType c = tp.c(eg.entity(),inside_local);
        
        // residual contribution
        mat(0,0) += c*eg.geometry().volume();
      }

	private:
	  TP& tp;
	};

  }
}

#endif
