// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_INSTATIONARYGRIDOPERATORSPACE_HH
#define DUNE_PDELAB_INSTATIONARYGRIDOPERATORSPACE_HH

#include<map>

#include<dune/common/exceptions.hh>
#include<dune/common/geometrytype.hh>
#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>

#include<dune/istl/io.hh>

#include"../common/geometrywrapper.hh"
#include"../gridfunctionspace/gridfunctionspace.hh"
#include"../gridfunctionspace/constraints.hh"
#include"../gridfunctionspace/interpolate.hh"
#include"localmatrix.hh"
#include"gridoperatorspaceutilities.hh"

namespace Dune {
  namespace PDELab {

    //! Base parameter class for time stepping scheme parameters
    /**
     * \tparam R C++ type of the floating point parameters
     */
    template<class R> 
    class TimeSteppingParameterInterface
    {
    public:
      typedef R RealType;
      
      /*! \brief Return true if method is implicit
      */
      virtual bool implicit () const = 0;
      
      /*! \brief Return number of stages of the method
      */
      virtual unsigned s () const = 0;
      
      /*! \brief Return entries of the A matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R a (int r, int i) const = 0;
      
      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R b (int r, int i) const = 0;
      
      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,r
      */
      virtual R d (int r) const = 0;
      
      /*! \brief Return name of the scheme
      */
      virtual std::string name () const = 0;
      
      //! every abstract base class has a virtual destructor
      virtual ~TimeSteppingParameterInterface () {}
    };


    //! Parameters specifying implicit euler
    /**
     * \tparam R C++ type of the floating point parameters
     */
    template<class R> 
    class ImplicitEulerParameter : public TimeSteppingParameterInterface<R>
    {
    public:
      
      ImplicitEulerParameter ()
      {
        D[0] = 0.0;  D[1] = 1.0;
        A[0][0] = -1.0; A[0][1] = 1.0;
        B[0][0] = 0.0;  B[0][1] = 1.0;
      }

      /*! \brief Return true if method is implicit
      */
      virtual bool implicit () const
      {
        return true;
      }
      
      /*! \brief Return number of stages s of the method
      */
      virtual unsigned s () const
      {
        return 1;
      }
      
      /*! \brief Return entries of the A matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R a (int r, int i) const
      {
        return A[r-1][i];
      }
      
      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R b (int r, int i) const
      {
        return B[r-1][i];
      }
      
      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R d (int i) const
      {
        return D[i];
      }
      
      /*! \brief Return name of the scheme
      */
      virtual std::string name () const
      {
        return std::string("implicit Euler");
      }

    private:
      Dune::FieldVector<R,2> D;
      Dune::FieldMatrix<R,1,2> A;
      Dune::FieldMatrix<R,1,2> B;
    };

	//================================================
	// The operator
	//================================================
	//! The generic assembler for time-dependent problems
    /**
     * \tparam TReal type to represent time values (and coefficients of time-stepping schemes)
     * \tparam R type that stores a residual vector
     * \tparam GFSU GridFunctionSpace for ansatz functions
     * \tparam GFSV GridFunctionSpace for test functions
     * \tparam LA   local operator assembler for spatial derivatives
     * \tparam LM   local operator assembler for temporal derivative
     * \tparam CU   assembled constraints for the space U
     * \tparam CV   assembled constraints for the space V
     * \tparam B    linear algebra backend
     * \tparam nonoverlapping_mode switch to assemble for nonoverlapping grids
     */
	template<typename TReal,
             typename R,
             typename GFSU, 
             typename GFSV, 
             typename LA, 
             typename LM,
			 typename CU=EmptyTransformation,
			 typename CV=EmptyTransformation,
			 typename B=StdVectorFlatMatrixBackend, 
             bool nonoverlapping_mode=false>
	class InstationaryGridOperatorSpace : public GridOperatorBase<GFSU,GFSV,CU,CV,B>
	{
	  // extract useful types
	  typedef typename GFSU::Traits::GridViewType GV;
      typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
	  typedef typename GV::Traits::template Codim<0>::Entity Element;
	  typedef typename GV::IntersectionIterator IntersectionIterator;
	  typedef typename IntersectionIterator::Intersection Intersection;

	public:
      typedef GridOperatorBase<GFSU,GFSV,CU,CV,B> Base;
	  typedef typename Base::Traits Traits;

      template<typename E>
      struct MatrixContainer
      {
        //! \brief define Type as the Type of a Matrix of E's
        typedef typename B::template Matrix<InstationaryGridOperatorSpace,E> Type;
      private:
        MatrixContainer () {}
      };

      //! construct 
	  InstationaryGridOperatorSpace (const TimeSteppingParameterInterface<TReal>& method_, 
                                     const GFSU& gfsu_, const GFSV& gfsv_, LA& la_, LM& lm_) 
		: Base(gfsu_,gfsv_), la(la_), lm(lm_), method(&method_), r0(gfsv,0.0)
	  {}

      //! construct using default time stepper
	  InstationaryGridOperatorSpace (const GFSU& gfsu_, const GFSV& gfsv_, LA& la_, LM& lm_) 
		: Base(gfsu_,gfsv_), la(la_), lm(lm_), method(&defaultmethod), r0(gfsv,0.0)
	  {}

      //! construct, with constraints
	  InstationaryGridOperatorSpace (const TimeSteppingParameterInterface<TReal>& method_, 
                                     const GFSU& gfsu_, const CU& cu,
                                     const GFSV& gfsv_, const CV& cv,
                                     LA& la_, LM& lm_) 
        : Base(gfsu_,cu,gfsv_,cv), la(la_), lm(lm_), method(&method_), r0(gfsv,0.0)
	  {}

      //! construct, with constraints and default time stepper
	  InstationaryGridOperatorSpace (const GFSU& gfsu_, const CU& cu,
                                     const GFSV& gfsv_, const CV& cv,
                                     LA& la_, LM& lm_) 
        : Base(gfsu_,cu,gfsv_,cv), la(la_), lm(lm_), method(&defaultmethod), r0(gfsv,0.0)
	  {}

      //! get dimension of space u
	  typename GFSU::Traits::SizeType globalSizeU () const
	  {
		return gfsu.globalSize();
	  }

      //! get dimension of space v
	  typename GFSV::Traits::SizeType globalSizeV () const
	  {
		return gfsv.globalSize();
	  }

      //! get the trial grid function space
      const GFSU& trialGridFunctionSpace() const
      {
        return gfsu;
      }

      //! get the test grid function space
      const GFSV& testGridFunctionSpace() const
      {
        return gfsv;
      }


      /**\brief Construct global sparsity pattern from local description

         This function can be called by the Matrix to get the sparsity pattern. 
         Assumes that the pattern is initially empty.
      */
      template<typename P>
      void fill_pattern (P& globalpattern) const
      {
 		// make local function spaces
		typedef typename GFSU::LocalFunctionSpace LFSU;
		LFSU lfsu(gfsu);
		typedef typename GFSV::LocalFunctionSpace LFSV;
		LFSV lfsv(gfsv);

        for (ElementIterator it = gfsu.gridview().template begin<0>();
             it!=gfsu.gridview().template end<0>(); ++it)
          {
 			// bind local function spaces to element
			lfsu.bind(*it);
			lfsv.bind(*it);

            LocalSparsityPattern localpattern;

            // get local pattern of spatial operator
            if (method->implicit())
              LocalAssemblerCallSwitch<LA,LA::doPatternVolume>::
                pattern_volume(la,lfsu,lfsv,localpattern);
            // add pattern of temporal operator
            LocalAssemblerCallSwitch<LM,LM::doPatternVolume>::
              pattern_volume(lm,lfsu,lfsv,localpattern);

            // skeleton and boundary pattern
            if((method->implicit() && (LA::doPatternSkeleton ||
                                       LA::doPatternBoundary)) ||
               LM::doPatternSkeleton || LM::doPatternBoundary)
            {

              // local function spaces in neighbor
              LFSU lfsun(gfsu);
              LFSV lfsvn(gfsv);

              IntersectionIterator endit = gfsu.gridview().iend(*it);
              for(IntersectionIterator iit = gfsu.gridview().ibegin(*it);
                  iit!=endit; ++iit)
              {
                // skeleton term
                if (iit->neighbor() &&
                    ((method->implicit() && LA::doPatternSkeleton) ||
                     LM::doPatternSkeleton) )
                {
                  // bind local function spaces to neighbor element
                  lfsun.bind(*(iit->outside()));
                  lfsvn.bind(*(iit->outside()));

                  LocalSparsityPattern localpattern_sn, localpattern_ns;
                  // spatial part
                  if(method->implicit())
                    LocalAssemblerCallSwitch<LA,LA::doPatternSkeleton>::
                      pattern_skeleton(la,lfsu,lfsv,lfsun,lfsvn,
                                       localpattern_sn, localpattern_ns);

                  // temporal part
                  LocalAssemblerCallSwitch<LM,LM::doPatternSkeleton>::
                    pattern_skeleton(lm,lfsu,lfsv,lfsun,lfsvn,
                                     localpattern_sn, localpattern_ns);

                  // translate local to global indices and add to global
                  // pattern
                  for (size_t k=0; k<localpattern_sn.size(); ++k)
                    add_entry(globalpattern,
                              lfsv.globalIndex(localpattern_sn[k].i()),
                              lfsun.globalIndex(localpattern_sn[k].j())
                              );

                  for (size_t k=0; k<localpattern_ns.size(); ++k)
                    add_entry(globalpattern,
                              lfsvn.globalIndex(localpattern_ns[k].i()),
                              lfsu.globalIndex(localpattern_ns[k].j())
                              );
                }
                // boundary term
                if (iit->boundary()) {
                  // spatial part
                  if(method->implicit())
                    LocalAssemblerCallSwitch<LA,LA::doPatternBoundary>::
                      pattern_boundary(la,lfsu,lfsv,localpattern);

                  // temporal part
                  LocalAssemblerCallSwitch<LM,LM::doPatternBoundary>::
                    pattern_boundary(lm,lfsu,lfsv,localpattern);
                }
			  }
            } // if((method->implicit() && LA::doPatternSkeleton) ||
              //    LM::doPatternSkeleton)

            // get local pattern of spatial operator
            if (method->implicit())
              LocalAssemblerCallSwitch<LA,LA::doPatternVolumePostSkeleton>::
                pattern_volume_post_skeleton(la,lfsu,lfsv,localpattern);
            // add pattern of temporal operator
            LocalAssemblerCallSwitch<LM,LM::doPatternVolumePostSkeleton>::
              pattern_volume_post_skeleton(lm,lfsu,lfsv,localpattern);

            // translate local to global indices and add to global pattern
            for (size_t k=0; k<localpattern.size(); ++k)
              add_entry(globalpattern,
                        lfsv.globalIndex(localpattern[k].i()),
                        lfsu.globalIndex(localpattern[k].j())
                        );
          } // Element loop
      }

      //! parametrize assembler with a time-stepping method
      void setMethod (const TimeSteppingParameterInterface<TReal>& method_)
      {
        method = &method_;
      }

      //! parametrize assembler with a time-stepping method
      /**
       * Invokes preStep(start_time, dt, nstages) on each local operator.
       */
      void preStep (const TimeSteppingParameterInterface<TReal>& method_, TReal time_, TReal dt_)
      {
        setMethod(method_);
        preStep(time_, dt_);
      }

      //! parametrize assembler with a time-stepping method
      /**
       * Invokes preStep(start_time, dt, nstages) on each local operator.
       */
      void preStep (TReal time_, TReal dt_)
      {
        time = time_;
        dt = dt_;
        la.preStep(time,dt,method->s());
        lm.preStep(time,dt,method->s());
      }

      //! to be called after step is completed
      /**
       * Invokes postStep() on the temporal local operator only.
       */
      void postStep ()
      {
        lm.postStep();
      }

      //! to be called after stage is completed
      /**
       * Invokes postStage() on the local operators.
       */
      void postStage ()
      {
        la.postStage();
        lm.postStage();
      }

      //! to be called once before each stage
      TReal suggestTimestep (TReal dt) const
      {
        TReal suggested_dt = la.suggestTimestep(dt);
        if (gfsu.gridview().comm().size()>1)
          suggested_dt =  gfsu.gridview().comm().min(suggested_dt);
        return suggested_dt;
      }

      //! Interpolate constrained values
      /**
       * \param stage Stage number in which to evaluate f.
       * \param xold  Vector with the old values.  Used to obtain the
       *              non-constrained values.
       * \param f     Function to evaluate to obtain the contrained values.
       * \param x     Where to store the combination of xold and the
       *              interpolated values.
       *
       * \note xold and x may not refer to the same object.
       *
       * Invokes setTime(time_of_stage) on f.
       */
	  template<typename F, typename X> 
      void interpolate (unsigned stage, const X& xold, F& f, X& x) const
      {
        // set time in boundary value function
        f.setTime(time+method->d(stage)*dt);

        // make x obey the boundary values
        Dune::PDELab::interpolate(f,gfsu,x);

        // copy non-constrained dofs from old time step
        Dune::PDELab::copy_nonconstrained_dofs(*pconstraintsv,xold,x);
      }

      //! set stage number to do next; assemble constant part of residual
      /**
       * Must be called before evaluating the residual for a certain stage.
       * Calls preStage() on the local operators.  Calls setTime() as
       * necessary on the local operators.
       */
	  template<typename X> 
      void preStage (unsigned stage_, const std::vector<X*>& x)
      {
        // process arguments
        stage = stage_;
        if (x.size()!=stage)
          DUNE_THROW(Exception,"wrong number of solutions in InstationaryGridOperatorSpace");
        if (stage<1 || stage>method->s())
          DUNE_THROW(Exception,"invalid stage number in InstationaryGridOperatorSpace");
 
        // visit each face only once
        const int chunk=1<<28;
        int offset = 0;
        const typename GV::IndexSet& is=gfsu.gridview().indexSet();
        std::map<Dune::GeometryType,int> gtoffset;

		// make local function spaces
		typedef typename GFSU::LocalFunctionSpace LFSU;
		LFSU lfsu(gfsu);
		typedef typename GFSV::LocalFunctionSpace LFSV;
		LFSV lfsv(gfsv);

        // extract coefficients of time stepping scheme
        std::vector<TReal> a(stage);
        for (size_t i=0; i<stage; ++i) a[i] = method->a(stage,i);
        std::vector<TReal> b(stage);
        for (size_t i=0; i<stage; ++i) b[i] = method->b(stage,i);
        std::vector<TReal> d(stage);
        for (size_t i=0; i<stage; ++i) d[i] = method->d(i);

        bool needsSkeleton = LA::doAlphaSkeleton||LA::doAlphaBoundary||LA::doLambdaSkeleton||LA::doLambdaBoundary;

        //std::cout << "preStage: stage " << stage << std::endl;

        // clear constant part residual before assembling
        r0 = 0.0;

        // prepare local operators for stage
        la.preStage(time+method->d(stage)*dt,stage);
        lm.preStage(time+method->d(stage)*dt,stage);

		// traverse grid view
		for (ElementIterator it = gfsu.gridview().template begin<0>();
			 it!=gfsu.gridview().template end<0>(); ++it)
		  {
            // assign offset for geometry type;
            if (gtoffset.find(it->type())==gtoffset.end())
              {
                gtoffset[it->type()] = offset;
                offset += chunk;
              }

            // compute unique id
            int id = is.index(*it)+gtoffset[it->type()];

            // skip ghost and overlap
            if (nonoverlapping_mode && it->partitionType()!=Dune::InteriorEntity)
              continue; 

			// bind local function spaces to element
			lfsu.bind(*it);
			lfsv.bind(*it);

            // loop over all previous time steps
            for (unsigned i=0; i<stage; ++i)
              {
                // set time in local operators for evaluation
                la.setTime(time+d[i]*dt);
                lm.setTime(time+d[i]*dt);

                // allocate local data container
                std::vector<typename X::ElementType> xl(lfsu.size());
                std::vector<typename R::ElementType> rl_a(lfsv.size(),0.0);
                std::vector<typename R::ElementType> rl_m(lfsv.size(),0.0);

                // read coefficents
                lfsu.vread(*x[i],xl);
                bool doM = a[i]>1E-6 || a[i]<-1E-6;
                bool doA = b[i]>1E-6 || b[i]<-1E-6;

                //std::cout << "R0 " << "stage=" << i << " time=" << time << " d_i*dt=" << d[i]*dt 
                //          << " doM=" << doM << " doA=" << doA << " skel=" << needsSkeleton << std::endl;

                // volume evaluation
                if (doA)
                  {
                    LocalAssemblerCallSwitch<LA,LA::doAlphaVolume>::
                      alpha_volume(la,ElementGeometry<Element>(*it),lfsu,xl,lfsv,rl_a);
                    LocalAssemblerCallSwitch<LA,LA::doLambdaVolume>::
                      lambda_volume(la,ElementGeometry<Element>(*it),lfsv,rl_a);
                  }
                if (doM)
                  {
                    LocalAssemblerCallSwitch<LM,LM::doAlphaVolume>::
                      alpha_volume(lm,ElementGeometry<Element>(*it),lfsu,xl,lfsv,rl_m);
                  }

                // skip if no intersection iterator is needed
                // note: LM has no skeleton and boundary terms !
                if (doA && needsSkeleton)
                  {
                    // local function spaces in neighbor
                    LFSU lfsun(gfsu);
                    LFSV lfsvn(gfsv);

                    // traverse intersections
                    unsigned int intersection_index = 0;
                    IntersectionIterator endit = gfsu.gridview().iend(*it);
                    for (IntersectionIterator iit = gfsu.gridview().ibegin(*it); 
                         iit!=endit; ++iit, ++intersection_index)
                      {
                        // skeleton term
                        if (iit->neighbor() && (LA::doAlphaSkeleton||LA::doLambdaSkeleton) )
                          {
                            // assign offset for geometry type;
                            Dune::GeometryType gtn = iit->outside()->type();
                            if (gtoffset.find(gtn)==gtoffset.end())
                              {
                                gtoffset[gtn] = offset;
                                offset += chunk;
                              }
                        
                            // compute unique id for neighbor
                            int idn = is.index(*(iit->outside()))+gtoffset[gtn];
                          
                            // unique vist of intersection
                            if (LA::doSkeletonTwoSided || id>idn || 
                                (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity) )
                              {
                                // bind local function spaces to neighbor element
                                lfsun.bind(*(iit->outside()));
                                lfsvn.bind(*(iit->outside()));
                            
                                // allocate local data container
                                std::vector<typename X::ElementType> xn(lfsun.size());
                                std::vector<typename R::ElementType> rn(lfsvn.size(),0.0);
                            
                                // read coefficents
                                lfsun.vread(*x[i],xn);
                            
                                // skeleton evaluation
                                LocalAssemblerCallSwitch<LA,LA::doAlphaSkeleton>::
                                  alpha_skeleton(la,IntersectionGeometry<Intersection>(*iit,intersection_index),lfsu,xl,lfsv,lfsun,xn,lfsvn,rl_a,rn);
                                LocalAssemblerCallSwitch<LA,LA::doLambdaSkeleton>::
                                  lambda_skeleton
                                  (la,IntersectionGeometry<Intersection>
                                             (*iit,intersection_index),
                                   lfsv,lfsvn,rl_a,rn);
                            
                                // accumulate result (note: r needs to be cleared outside)
                                for (size_t k=0; k<rn.size(); ++k) rn[k] *= b[i]*dt;
                                lfsvn.vadd(rn,r0);
                              }
                          }
               
                        // boundary term
                        if (iit->boundary())
                          {
                            LocalAssemblerCallSwitch<LA,LA::doAlphaBoundary>::
                              alpha_boundary(la,IntersectionGeometry<Intersection>(*iit,intersection_index),lfsu,xl,lfsv,rl_a);
                            LocalAssemblerCallSwitch<LA,LA::doLambdaBoundary>::
                              lambda_boundary(la,IntersectionGeometry<Intersection>(*iit,intersection_index),lfsv,rl_a);
                          }
                      }
                  }

                if (doA)
                  {
                    LocalAssemblerCallSwitch<LA,LA::doAlphaVolumePostSkeleton>::
                      alpha_volume_post_skeleton(la,ElementGeometry<Element>(*it),lfsu,xl,lfsv,rl_a);
                    LocalAssemblerCallSwitch<LA,LA::doLambdaVolumePostSkeleton>::
                      lambda_volume_post_skeleton(la,ElementGeometry<Element>(*it),lfsv,rl_a);
                    
                    // accumulate result (note: r needs to be cleared outside)
                    for (size_t k=0; k<rl_a.size(); ++k) rl_a[k] *= b[i]*dt; 
                    lfsv.vadd(rl_a,r0);
                  }
                if (doM)
                  {
                    for (size_t k=0; k<rl_m.size(); ++k) rl_m[k] *= a[i]; 
                    lfsv.vadd(rl_m,r0);
                  }
              }
          }

        //Dune::printvector(std::cout,r.base(),"const residual","row",4,9,1);
      }


      //! set stage number to do next; assemble constant part of residual
      /**
       * This is essentially a combination of preStage() and residual() for
       * the case of an explicit jacobian.  It is mainly used to determine the
       * time step size from the matrix.
       *
       * in explicit mode we assume that 
       *  A) the problem is linear in the d_t term
       *  B) the jacobian is block diagonal
       * this means that the system can always be solved by one step of a Jacobi preconditioner
       * without even checking for the residual.
       * From B also follows that the time local operator has only alpha_volume
       * \param[in] stage_ the stage we are in
       * \param[in] vector of pointers to the solutions in previous stages
       * \param[out] mat the block diagonal Jacobian to be assembled; we assume it is zero on entry!
       * \param[out] alpha temporal part of the residual; we assume it is zero on entry!
       * \param[out] beta spatial part of residual; we assume it is zero on entry!
       *
       * Calls preStage() on the local operators, and setTime() as
       * apropriate.  Assumes that preStep() has been called before.
       */
	  template<typename X, typename A> 
      void explicit_jacobian_residual (unsigned stage_, const std::vector<X*>& x, A& mat, R& alpha, R& beta)
      {
        // process arguments
        stage = stage_;
        if (x.size()!=stage+1)
          DUNE_THROW(Exception,"wrong number of solutions in InstationaryGridOperatorSpace");
        if (stage<1 || stage>method->s())
          DUNE_THROW(Exception,"invalid stage number in InstationaryGridOperatorSpace");
        if (method->implicit())
          DUNE_THROW(Exception,"explicit mode called with implicit scheme");
 
        // visit each face only once
        const int chunk=1<<28;
        int offset = 0;
        const typename GV::IndexSet& is=gfsu.gridview().indexSet();
        std::map<Dune::GeometryType,int> gtoffset;

		// make local function spaces
		typedef typename GFSU::LocalFunctionSpace LFSU;
		LFSU lfsu(gfsu);
		typedef typename GFSV::LocalFunctionSpace LFSV;
		LFSV lfsv(gfsv);

        // extract coefficients of time stepping scheme
        std::vector<TReal> a(stage);
        for (size_t i=0; i<stage; ++i) a[i] = method->a(stage,i);
        std::vector<TReal> b(stage);
        for (size_t i=0; i<stage; ++i) b[i] = method->b(stage,i);
        std::vector<TReal> d(stage);
        for (size_t i=0; i<stage; ++i) d[i] = method->d(i);
        TReal d_r = method->d(stage);

        bool needsSkeleton = LA::doAlphaSkeleton||LA::doAlphaBoundary||LA::doLambdaSkeleton||LA::doLambdaBoundary;

        // prepare local operators for stage
        la.preStage(time+method->d(stage)*dt,stage);
        lm.preStage(time+method->d(stage)*dt,stage);

		// traverse grid view
		for (ElementIterator it = gfsu.gridview().template begin<0>();
			 it!=gfsu.gridview().template end<0>(); ++it)
		  {
            // assign offset for geometry type;
            if (gtoffset.find(it->type())==gtoffset.end())
              {
                gtoffset[it->type()] = offset;
                offset += chunk;
              }

            // compute unique id
            int id = is.index(*it)+gtoffset[it->type()];

            // skip ghost and overlap
            if (nonoverlapping_mode && it->partitionType()!=Dune::InteriorEntity)
              continue; 

			// bind local function spaces to element
			lfsu.bind(*it);
			lfsv.bind(*it);

            // residual part
            // loop over all previous time steps (stages)
            for(unsigned i=0; i<stage; ++i)
              {
                // set time in local operators for evaluation
                la.setTime(time+d[i]*dt);
                lm.setTime(time+d[i]*dt);

                // allocate local data container
                std::vector<typename X::ElementType> xl(lfsu.size());
                std::vector<typename R::ElementType> rl_a(lfsv.size(),0.0);
                std::vector<typename R::ElementType> rl_m(lfsv.size(),0.0);

                // read coefficents
                lfsu.vread(*x[i],xl);
                bool doM = a[i]>1E-6 || a[i]<-1E-6;
                bool doA = b[i]>1E-6 || b[i]<-1E-6;

                //std::cout << "R0 " << "stage=" << i << " time=" << time << " d_i*dt=" << d[i]*dt 
                //          << " doM=" << doM << " doA=" << doA << " skel=" << needsSkeleton << std::endl;

                // volume evaluation
                if (doA)
                  {
                    LocalAssemblerCallSwitch<LA,LA::doAlphaVolume>::
                      alpha_volume(la,ElementGeometry<Element>(*it),lfsu,xl,lfsv,rl_a);
                    LocalAssemblerCallSwitch<LA,LA::doLambdaVolume>::
                      lambda_volume(la,ElementGeometry<Element>(*it),lfsv,rl_a);
                  }
                if (doM)
                  {
                    LocalAssemblerCallSwitch<LM,LM::doAlphaVolume>::
                      alpha_volume(lm,ElementGeometry<Element>(*it),lfsu,xl,lfsv,rl_m);
                  }

                // skip if no intersection iterator is needed
                // note: LM has no skeleton and boundary terms !
                if (doA && needsSkeleton)
                  {
                    // local function spaces in neighbor
                    LFSU lfsun(gfsu);
                    LFSV lfsvn(gfsv);

                    // traverse intersections
                    unsigned int intersection_index = 0;
                    IntersectionIterator endit = gfsu.gridview().iend(*it);
                    for (IntersectionIterator iit = gfsu.gridview().ibegin(*it); 
                         iit!=endit; ++iit, ++intersection_index)
                      {
                        // skeleton term
                        if (iit->neighbor() && (LA::doAlphaSkeleton||LA::doLambdaSkeleton) )
                          {
                            // assign offset for geometry type;
                            Dune::GeometryType gtn = iit->outside()->type();
                            if (gtoffset.find(gtn)==gtoffset.end())
                              {
                                gtoffset[gtn] = offset;
                                offset += chunk;
                              }
                        
                            // compute unique id for neighbor
                            int idn = is.index(*(iit->outside()))+gtoffset[gtn];
                          
                            // unique vist of intersection
                            if (LA::doSkeletonTwoSided || id>idn || 
                                (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity) )
                              {
                                // bind local function spaces to neighbor element
                                lfsun.bind(*(iit->outside()));
                                lfsvn.bind(*(iit->outside()));
                            
                                // allocate local data container
                                std::vector<typename X::ElementType> xn(lfsun.size());
                                std::vector<typename R::ElementType> rn(lfsvn.size(),0.0);
                            
                                // read coefficents
                                lfsun.vread(*x[i],xn);
                            
                                // skeleton evaluation
                                LocalAssemblerCallSwitch<LA,LA::doAlphaSkeleton>::
                                  alpha_skeleton(la,IntersectionGeometry<Intersection>(*iit,intersection_index),lfsu,xl,lfsv,lfsun,xn,lfsvn,rl_a,rn);
                                LocalAssemblerCallSwitch<LA,LA::doLambdaSkeleton>::
                                  lambda_skeleton
                                  (la,IntersectionGeometry<Intersection>
                                           (*iit,intersection_index),
                                   lfsv,lfsvn,rl_a,rn);
                            
                                // accumulate result (note: r needs to be cleared outside)
                                for (size_t k=0; k<rn.size(); ++k) rn[k] *= -b[i];
                                lfsvn.vadd(rn,beta);
                              }
                          }
               
                        // boundary term
                        if (iit->boundary())
                          {
                            LocalAssemblerCallSwitch<LA,LA::doAlphaBoundary>::
                              alpha_boundary(la,IntersectionGeometry<Intersection>(*iit,intersection_index),lfsu,xl,lfsv,rl_a);
                            LocalAssemblerCallSwitch<LA,LA::doLambdaBoundary>::
                              lambda_boundary(la,IntersectionGeometry<Intersection>(*iit,intersection_index),lfsv,rl_a);
                          }
                      }
                  }

                if (doA)
                  {
                    LocalAssemblerCallSwitch<LA,LA::doAlphaVolumePostSkeleton>::
                      alpha_volume_post_skeleton(la,ElementGeometry<Element>(*it),lfsu,xl,lfsv,rl_a);
                    LocalAssemblerCallSwitch<LA,LA::doLambdaVolumePostSkeleton>::
                      lambda_volume_post_skeleton(la,ElementGeometry<Element>(*it),lfsv,rl_a);
                    
                    // accumulate result (note: beta needs to be cleared outside)
                    for (size_t k=0; k<rl_a.size(); ++k) rl_a[k] *= -b[i]; 
                    lfsv.vadd(rl_a,beta);
                  }
                if (doM)
                  {
                    for (size_t k=0; k<rl_m.size(); ++k) rl_m[k] *= -a[i]; 
                    lfsv.vadd(rl_m,alpha);
                  }
              }

            // Jacobian part
            // Note: 
            // - we are explicit; there is no spatial part here
            // - temporal part has only alpha_volume

			// allocate local data container
			std::vector<typename X::ElementType> xl(lfsu.size());
			LocalMatrix<typename A::ElementType> ml(lfsv.size(),lfsu.size(),0.0);

            // set time in local operator for evaluation
            lm.setTime(time+d_r*dt);

            // read coefficents; this is only a dummy since Jacobian should not depend on solution !
            // but of course it is required to give this parameter
            lfsu.vread(*x[stage],xl);

            // compute local jacobian
            LocalAssemblerCallSwitch<LM,LM::doAlphaVolume>::
              jacobian_volume(lm,ElementGeometry<Element>(*it),lfsu,xl,lfsv,ml);

			// accumulate to global matrix
            etadd(lfsv,lfsu,ml,mat); // scheme is normalized 
          }

        // set trivial conditions for constrained degrees of freedom
        typedef typename CV::const_iterator global_row_iterator;	  
        for (global_row_iterator cit=pconstraintsv->begin(); cit!=pconstraintsv->end(); ++cit)
          set_trivial_row(cit->first,cit->second,mat);
        
		// set residual to zero on constrained dofs of spatial part (which is scaled by dt)
		Dune::PDELab::constrain_residual(*pconstraintsv,beta);

		// copy solution on constrained dofs from solution of stage to temporal part (which is not scaled)
        // this makes the boundary conditions appear in the solution !
		Dune::PDELab::copy_constrained_dofs(*pconstraintsu,*x[stage],alpha);

        //Dune::printvector(std::cout,r.base(),"const residual","row",4,9,1);
      }

	  //! generic evaluation of residual
      /**
       * \param r residual (needs to be cleared before this method is called)
       *
       * Invokes setTime(time_of_current_stage) on the local operators.
       * preStage() must have been called before this method to assemble the
       * constant part of the residual and to set the current stage number.
       */
	  template<typename X> 
	  void residual (const X& x, R& r) const
	  {
        //Dune::printvector(std::cout,x.base(),"x on entry to residual","row",4,9,1);

        // copy constant part of residual
        r = r0; // assumes assignment operator on vectors.
        //Dune::printvector(std::cout,r.base(),"r after copy in residual","row",4,9,1);

        // visit each face only once
        const int chunk=1<<28;
        int offset = 0;
        const typename GV::IndexSet& is=gfsu.gridview().indexSet();
        std::map<Dune::GeometryType,int> gtoffset;

		// make local function spaces
		typedef typename GFSU::LocalFunctionSpace LFSU;
		LFSU lfsu(gfsu);
		typedef typename GFSV::LocalFunctionSpace LFSV;
		LFSV lfsv(gfsv);

        // extract coefficients of time stepping scheme
        TReal b_rr = method->b(stage,stage);
        TReal d_r = method->d(stage);
        bool implicit = method->implicit();

        // set time in local operators for evaluation
        la.setTime(time+d_r*dt);
        lm.setTime(time+d_r*dt);

		// traverse grid view
		for (ElementIterator it = gfsu.gridview().template begin<0>();
			 it!=gfsu.gridview().template end<0>(); ++it)
		  {
            // assign offset for geometry type;
            if (gtoffset.find(it->type())==gtoffset.end())
              {
                gtoffset[it->type()] = offset;
                offset += chunk;
              }

            // compute unique id
            int id = is.index(*it)+gtoffset[it->type()];

            // skip ghost and overlap
            if (nonoverlapping_mode && it->partitionType()!=Dune::InteriorEntity)
              continue; 

			// bind local function spaces to element
			lfsu.bind(*it);
			lfsv.bind(*it);

			// allocate local data container
			std::vector<typename X::ElementType> xl(lfsu.size());
			std::vector<typename R::ElementType> rl_a(lfsv.size(),0.0);
			std::vector<typename R::ElementType> rl_m(lfsv.size(),0.0);

			// read coefficents
			lfsu.vread(x,xl);

			// volume evaluation
            if (implicit)
              {
                LocalAssemblerCallSwitch<LA,LA::doAlphaVolume>::
                  alpha_volume(la,ElementGeometry<Element>(*it),lfsu,xl,lfsv,rl_a);
                LocalAssemblerCallSwitch<LA,LA::doLambdaVolume>::
                  lambda_volume(la,ElementGeometry<Element>(*it),lfsv,rl_a);
              }
			LocalAssemblerCallSwitch<LM,LM::doAlphaVolume>::
              alpha_volume(lm,ElementGeometry<Element>(*it),lfsu,xl,lfsv,rl_m);

            //std::cout << "residual " << "stage=" << stage << " time=" << time << " d_i*dt=" << d_r*dt 
            //          << " b_rr=" << b_rr << " implicit=" << implicit << std::endl;

			// skip if no intersection iterator is needed
 			if (implicit&&(LA::doAlphaSkeleton||LA::doAlphaBoundary||LA::doLambdaSkeleton||LA::doLambdaBoundary))
              {
                // local function spaces in neighbor
                LFSU lfsun(gfsu);
                LFSV lfsvn(gfsv);

                // traverse intersections
                unsigned int intersection_index = 0;
                IntersectionIterator endit = gfsu.gridview().iend(*it);
                for (IntersectionIterator iit = gfsu.gridview().ibegin(*it); 
                     iit!=endit; ++iit, ++intersection_index)
                  {
                    // skeleton term
                    if (iit->neighbor() && (LA::doAlphaSkeleton||LA::doLambdaSkeleton) )
                      {
                        // assign offset for geometry type;
                        Dune::GeometryType gtn = iit->outside()->type();
                        if (gtoffset.find(gtn)==gtoffset.end())
                          {
                            gtoffset[gtn] = offset;
                            offset += chunk;
                          }
                        
                        // compute unique id for neighbor
                        int idn = is.index(*(iit->outside()))+gtoffset[gtn];
                          
                        // unique vist of intersection
                        if (LA::doSkeletonTwoSided || id>idn || 
                            (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity) )
                          {
                            // bind local function spaces to neighbor element
                            lfsun.bind(*(iit->outside()));
                            lfsvn.bind(*(iit->outside()));
                            
                            // allocate local data container
                            std::vector<typename X::ElementType> xn(lfsun.size());
                            std::vector<typename R::ElementType> rn(lfsvn.size(),0.0);
                            
                            // read coefficents
                            lfsun.vread(x,xn);
                            
                            // skeleton evaluation
                            LocalAssemblerCallSwitch<LA,LA::doAlphaSkeleton>::
                              alpha_skeleton(la,IntersectionGeometry<Intersection>(*iit,intersection_index),lfsu,xl,lfsv,lfsun,xn,lfsvn,rl_a,rn);
                            LocalAssemblerCallSwitch<LA,LA::doLambdaSkeleton>::
                              lambda_skeleton
                              (la,IntersectionGeometry<Intersection>
                                       (*iit,intersection_index),
                               lfsv,lfsvn,rl_a,rn);
                            
                            // accumulate result (note: r needs to be cleared outside)
                            for (size_t i=0; i<rn.size(); ++i) rn[i] *= b_rr*dt;
                            lfsvn.vadd(rn,r);
                          }
                      }
               
                    // boundary term
                    if (iit->boundary())
                      {
                        LocalAssemblerCallSwitch<LA,LA::doAlphaBoundary>::
                          alpha_boundary(la,IntersectionGeometry<Intersection>(*iit,intersection_index),lfsu,xl,lfsv,rl_a);
                        LocalAssemblerCallSwitch<LA,LA::doLambdaBoundary>::
                          lambda_boundary(la,IntersectionGeometry<Intersection>(*iit,intersection_index),lfsv,rl_a);
                      }
                  }
              }

            if (implicit)
              {
                LocalAssemblerCallSwitch<LA,LA::doAlphaVolumePostSkeleton>::
                  alpha_volume_post_skeleton(la,ElementGeometry<Element>(*it),lfsu,xl,lfsv,rl_a);
                LocalAssemblerCallSwitch<LA,LA::doLambdaVolumePostSkeleton>::
                  lambda_volume_post_skeleton(la,ElementGeometry<Element>(*it),lfsv,rl_a);

                // accumulate result (note: r needs to be cleared outside)
                for (size_t i=0; i<rl_a.size(); ++i) rl_a[i] *= b_rr*dt; 
                lfsv.vadd(rl_a,r);
              }
 
			lfsv.vadd(rl_m,r); // scheme is normalized !
		  }

		// set residual to zero on constrained dofs
		Dune::PDELab::constrain_residual(*pconstraintsv,r);

        //Dune::printvector(std::cout,r.base(),"residual","row",4,9,1);
	  }

	  //! generic application of Jacobian
      /**
       * Invokes setTime(time_of_current_stage) on the local operators.
       * preStage() must have been called before this method to set the
       * current stage number.
       */
	  template<typename X, typename Y> 
	  void jacobian_apply (X& x, Y& y) const
	  {
        // visit each face only once
        const int chunk=1<<28;
        int offset = 0;
        const typename GV::IndexSet& is=gfsu.gridview().indexSet();
        std::map<Dune::GeometryType,int> gtoffset;

		// make local function spaces
		typedef typename GFSU::LocalFunctionSpace LFSU;
		LFSU lfsu(gfsu);
		typedef typename GFSV::LocalFunctionSpace LFSV;
		LFSV lfsv(gfsv);

        // extract coefficients of time stepping scheme
        TReal b_rr = method->b(stage,stage);
        TReal d_r = method->d(stage);
        bool implicit = method->implicit();

        // set time in local operators for evaluation
        la.setTime(time+d_r*dt);
        lm.setTime(time+d_r*dt);

		// traverse grid view
		for (ElementIterator it = gfsu.gridview().template begin<0>();
			 it!=gfsu.gridview().template end<0>(); ++it)
		  {
            // assign offset for geometry type;
            if (gtoffset.find(it->type())==gtoffset.end())
              {
                gtoffset[it->type()] = offset;
                offset += chunk;
              }

            // compute unique id
            int id = is.index(*it)+gtoffset[it->type()];

            // skip ghost and overlap
            if (nonoverlapping_mode && it->partitionType()!=Dune::InteriorEntity)
              continue; 

			// bind local function spaces to element
			lfsu.bind(*it);
			lfsv.bind(*it);

			// allocate local data container
			std::vector<typename X::ElementType> xl(lfsu.size());
			std::vector<typename Y::ElementType> yl_a(lfsv.size(),0.0);
			std::vector<typename Y::ElementType> yl_m(lfsv.size(),0.0);

			// read coefficents
			lfsu.vread(x,xl);

			// volume evaluation
            if (implicit)
              {
                LocalAssemblerCallSwitch<LA,LA::doAlphaVolume>::
                  jacobian_apply_volume(la,ElementGeometry<Element>(*it),lfsu,xl,lfsv,yl_a);
              }
            LocalAssemblerCallSwitch<LM,LM::doAlphaVolume>::
              jacobian_apply_volume(lm,ElementGeometry<Element>(*it),lfsu,xl,lfsv,yl_m);

			// skeleton and boundary evaluation
			if (implicit&&(LA::doAlphaSkeleton||LA::doAlphaBoundary))
			  {
                // local function spaces in neighbor
                LFSU lfsun(gfsu);
                LFSV lfsvn(gfsv);

                unsigned int intersection_index = 0;
				IntersectionIterator endit = gfsu.gridview().iend(*it);
				for (IntersectionIterator iit = gfsu.gridview().ibegin(*it); 
					 iit!=endit; ++iit, ++intersection_index)
				  {
                    // skeleton term
                    if (iit->neighbor() && LA::doAlphaSkeleton )
                      {
                        // assign offset for geometry type;
                        Dune::GeometryType gtn = iit->outside()->type();
                        if (gtoffset.find(gtn)==gtoffset.end())
                          {
                            gtoffset[gtn] = offset;
                            offset += chunk;
                          }
                        
                        // compute unique id for neighbor
                        int idn = is.index(*(iit->outside()))+gtoffset[gtn];
                          
                        // unique vist of intersection
                        if (LA::doSkeletonTwoSided || id>idn ||
                            (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity))
                          {
                            // bind local function spaces to neighbor element
                            lfsun.bind(*(iit->outside()));
                            lfsvn.bind(*(iit->outside()));
                            
                            // allocate local data container
                            std::vector<typename X::ElementType> xn(lfsun.size());
                            std::vector<typename Y::ElementType> yn(lfsvn.size(),0.0);
                            
                            // read coefficents
                            lfsun.vread(x,xn);
                            
                            // skeleton evaluation
                            LocalAssemblerCallSwitch<LA,LA::doAlphaSkeleton>::
                              jacobian_apply_skeleton(la,IntersectionGeometry<Intersection>(*iit,intersection_index),lfsu,xl,lfsv,lfsun,xn,lfsvn,yl_a,yn);

                            // accumulate result (note: r needs to be cleared outside)
                            for (size_t i=0; i<yn.size(); ++i) yn[i] *= b_rr*dt; 
                            lfsvn.vadd(yn,y);
                          }
                      }

                    // boundary term
                    if (iit->boundary())
                      {
                        LocalAssemblerCallSwitch<LA,LA::doAlphaBoundary>::
                          jacobian_apply_boundary(la,IntersectionGeometry<Intersection>(*iit,intersection_index),lfsu,xl,lfsv,yl_a);
                      }
				  }
			  }

            if (implicit)
              {
                LocalAssemblerCallSwitch<LA,LA::doAlphaVolumePostSkeleton>::
                  jacobian_apply_volume_post_skeleton(la,ElementGeometry<Element>(*it),lfsu,xl,lfsv,yl_a);
                for (size_t i=0; i<yl_a.size(); ++i) yl_a[i] *= b_rr*dt; 
                lfsv.vadd(yl_a,y);
              }

			// accumulate result (note: r needs to be cleared outside)
			lfsv.vadd(yl_m,y); // scheme is normalized!
		  }

		// set residual to zero on constrained dofs
		Dune::PDELab::copy_constrained_dofs(*pconstraintsu,x,y);
	  }

	  //! generic assembly of Jacobian
      /**
       * \param x Where (in the space spanned by the dofs) to evaluate the Jacobian
       * \param a Jacobian (needs to be cleared before passed to this method)
       *
       * Invokes setTime(time_of_current_stage) on the local operators.
       * preStage() must have been called before this method to set the
       * current stage number.
       */
	  template<typename X, typename A> 
	  void jacobian (const X& x, A& a) const
	  {
        // visit each face only once
        const int chunk=1<<28;
        int offset = 0;
        const typename GV::IndexSet& is=gfsu.gridview().indexSet();
        std::map<Dune::GeometryType,int> gtoffset;

		// make local function spaces
		typedef typename GFSU::LocalFunctionSpace LFSU;
		LFSU lfsu(gfsu);
		typedef typename GFSV::LocalFunctionSpace LFSV;
		LFSV lfsv(gfsv);

        // extract coefficients of time stepping scheme
        TReal b_rr = method->b(stage,stage);
        TReal d_r = method->d(stage);
        bool implicit = method->implicit();

        // set time in local operators for evaluation
        la.setTime(time+d_r*dt);
        lm.setTime(time+d_r*dt);

		// traverse grid view
		for (ElementIterator it = gfsu.gridview().template begin<0>();
			 it!=gfsu.gridview().template end<0>(); ++it)
		  {
            // assign offset for geometry type;
            if (gtoffset.find(it->type())==gtoffset.end())
              {
                gtoffset[it->type()] = offset;
                offset += chunk;
              }

            // compute unique id
            const typename GV::IndexSet::IndexType id = is.index(*it)+gtoffset[it->type()];
            //            std::cout << "[" << gfsu.gridview().comm().rank() << "] " << " element: " << id << std::endl;

            // skip ghost and overlap
            if (nonoverlapping_mode && it->partitionType()!=Dune::InteriorEntity)
              continue; 
 
			// bind local function spaces to element
			lfsu.bind(*it);
			lfsv.bind(*it);

			// allocate local data container
			std::vector<typename X::ElementType> xl(lfsu.size());
			LocalMatrix<typename A::ElementType> al(lfsv.size(),lfsu.size(),0.0);
			LocalMatrix<typename A::ElementType> ml(lfsv.size(),lfsu.size(),0.0);

			// read coefficents
			lfsu.vread(x,xl);

			// volume evaluation
            if (implicit)
              {
                LocalAssemblerCallSwitch<LA,LA::doAlphaVolume>::
                  jacobian_volume(la,ElementGeometry<Element>(*it),lfsu,xl,lfsv,al);
              }
            LocalAssemblerCallSwitch<LM,LM::doAlphaVolume>::
              jacobian_volume(lm,ElementGeometry<Element>(*it),lfsu,xl,lfsv,ml);

			// skeleton and boundary evaluation
			if (implicit&&(LA::doAlphaSkeleton||LA::doAlphaBoundary))
			  {
                // local function spaces in neighbor
                LFSU lfsun(gfsu);
                LFSV lfsvn(gfsv);
                
                unsigned int intersection_index = 0;
				IntersectionIterator endit = gfsu.gridview().iend(*it);
				for (IntersectionIterator iit = gfsu.gridview().ibegin(*it); 
					 iit!=endit; ++iit, ++intersection_index)
				  {
                    // skeleton term
                    if (iit->neighbor() && LA::doAlphaSkeleton )
                      {
                        // assign offset for geometry type;
                        Dune::GeometryType gtn = iit->outside()->type();
                        if (gtoffset.find(gtn)==gtoffset.end())
                          {
                            gtoffset[gtn] = offset;
                            offset += chunk;
                          }
                        
                        // compute unique id for neighbor
                        const typename GV::IndexSet::IndexType idn = is.index(*(iit->outside()))+gtoffset[gtn];
                          
                        // unique vist of intersection
                        if (LA::doSkeletonTwoSided || id>idn ||
                            (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity) )
                          {
                            // bind local function spaces to neighbor element
                            lfsun.bind(*(iit->outside()));
                            lfsvn.bind(*(iit->outside()));
                            
                            // allocate local data container
                            std::vector<typename X::ElementType> xn(lfsun.size());
                            LocalMatrix<typename A::ElementType> al_sn(lfsv.size() ,lfsun.size(),0.0);
                            LocalMatrix<typename A::ElementType> al_ns(lfsvn.size(),lfsu.size() ,0.0);
                            LocalMatrix<typename A::ElementType> al_nn(lfsvn.size(),lfsun.size(),0.0);
                            
                            // read coefficents
                            lfsun.vread(x,xn);
                            
                            // skeleton evaluation
                            LocalAssemblerCallSwitch<LA,LA::doAlphaSkeleton>::
                              jacobian_skeleton(la,IntersectionGeometry<Intersection>(*iit,intersection_index),
                                                lfsu,xl,lfsv,lfsun,xn,lfsvn,al,al_sn,al_ns,al_nn);

                            // accumulate result
                            al_sn *= b_rr*dt; etadd(lfsv,lfsun,al_sn,a);
                            al_ns *= b_rr*dt; etadd(lfsvn,lfsu,al_ns,a);
                            al_nn *= b_rr*dt; etadd(lfsvn,lfsun,al_nn,a);
                          }
                      }

                    // boundary term
                    if (iit->boundary())
                      {
                        LocalAssemblerCallSwitch<LA,LA::doAlphaBoundary>::
                          jacobian_boundary(la,IntersectionGeometry<Intersection>(*iit,intersection_index),lfsu,xl,lfsv,al);
                      }
				  }
			  }

            if (implicit)
              {
                LocalAssemblerCallSwitch<LA,LA::doAlphaVolumePostSkeleton>::
                  jacobian_volume_post_skeleton(la,ElementGeometry<Element>(*it),lfsu,xl,lfsv,al);
                al *= b_rr*dt;
                etadd(lfsv,lfsu,al,a);
              }

			// accumulate result (note: a needs to be cleared outside)
            etadd(lfsv,lfsu,ml,a); // scheme is normalized 
		  }

         typedef typename CV::const_iterator global_row_iterator;	  
         for (global_row_iterator cit=pconstraintsv->begin(); cit!=pconstraintsv->end(); ++cit)
           set_trivial_row(cit->first,cit->second,a);

         //printmatrix(std::cout,a.base(),"global stiffness matrix","row",9,1);
 	  }

      using Base::gfsu;
      using Base::gfsv;
      using Base::pconstraintsu;
      using Base::pconstraintsv;
	  LA& la;
      LM& lm;
      const TimeSteppingParameterInterface<TReal> *method;
      TReal time, dt;
      unsigned stage;
      R r0;
      ImplicitEulerParameter<TReal> defaultmethod;
	};

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
