// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDOPERATORSPACE_HH
#define DUNE_PDELAB_GRIDOPERATORSPACE_HH

#include<map>

#include<dune/common/deprecated.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/geometrytype.hh>

#include "../common/geometrywrapper.hh"
#include "../constraints/constraints.hh"
#include "../gridfunctionspace/gridfunctionspace.hh"
#include "../gridfunctionspace/localvector.hh"
#include "localmatrix.hh"

#include "gridoperatorspaceutilities.hh"


namespace Dune {
  namespace PDELab {

	//================================================
	// The operator
	//================================================

    /**
       \brief The generic assembler

     * \tparam GFSU GridFunctionSpace for ansatz functions
     * \tparam GFSV GridFunctionSpace for test functions
     * \tparam LA   local operator assembler (provided by user)
     * \tparam CU   Constraints maps for the individual dofs (trial space)
     * \tparam CV   Constraints maps for the individual dofs (test space)
     * \tparam B The vector backend used for the coefficient vector
     * \tparam nonoverlapping_mode Indicates whether assembling is done for overlap cells
     */
    /**
       \deprecated This class is deprecated, please use Dune::PDELab::GridOperator
    */
	template<typename GFSU, typename GFSV, typename LA,
			 typename CU=EmptyTransformation,
			 typename CV=EmptyTransformation,
			 typename B=StdVectorFlatMatrixBackend, 
             bool nonoverlapping_mode=false,
             typename ST=NoSubTriangulation<typename GFSU::Traits::GridViewType> >
	class GridOperatorSpace : public GridOperatorBase<GFSU,GFSV,CU,CV,B>
	{
	  // extract useful types
	  typedef typename GFSU::Traits::GridViewType GV;
      typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
	  typedef typename GV::Traits::template Codim<0>::Entity Element;
	  typedef typename GV::IntersectionIterator IntersectionIterator;
	  typedef typename IntersectionIterator::Intersection Intersection;
      typedef typename ST::SubIntersectionIterator SubIntersectionIterator;
      typedef typename ST::SubEntityIterator SubEntityIterator;
      typedef ST SubTriangulation;

	public:
      typedef GridOperatorBase<GFSU,GFSV,CU,CV,B> Base;
	  typedef typename Base::Traits Traits;

      //! construct GridOperatorSpace
      /**
         \deprecated This class is deprecated, please use Dune::PDELab::GridOperator
       */
	  GridOperatorSpace (const GFSU& gfsu_, const GFSV& gfsv_, 
                         const LA & la_, const ST & st_) DUNE_DEPRECATED
		: Base(gfsu_,gfsv_), la(la_), sub_triangulation(st_)
	  { }

      //! construct GridOperatorSpace
      /**
         \deprecated This class is deprecated, please use Dune::PDELab::GridOperator
       */
	  GridOperatorSpace (const GFSU& gfsu_, const GFSV& gfsv_, 
                         const LA & la_) DUNE_DEPRECATED
		: Base(gfsu_,gfsv_), la(la_), sub_triangulation(ST(gfsu_.gridview(),NoSubTriangulationImp()))
	  { }

      //! construct GridOperatorSpace, with constraints
      /**
         \deprecated This class is deprecated, please use Dune::PDELab::GridOperator
       */
	  GridOperatorSpace (const GFSU& gfsu_, const CU& cu,
						 const GFSV& gfsv_, const CV& cv, 
                         const LA & la_, const ST & st_) DUNE_DEPRECATED
		: Base(gfsu_,cu,gfsv_,cv), la(la_), sub_triangulation(st_)
	  { }

      //! construct GridOperatorSpace, with constraints
      /**
         \deprecated This class is deprecated, please use Dune::PDELab::GridOperator
       */
	  GridOperatorSpace (const GFSU& gfsu_, const CU& cu,
						 const GFSV& gfsv_, const CV& cv, 
                         const LA & la_) DUNE_DEPRECATED
		: Base(gfsu_,cu,gfsv_,cv), la(la_), sub_triangulation(ST(gfsu_.gridview(),NoSubTriangulationImp()))
	  { }


      template<typename E>
      struct MatrixContainer
      {
        //! \brief define Type as the Type of a Matrix of E's
        typedef typename B::template Matrix<E> Type;
      private:
        MatrixContainer () {}
      };

      /**\brief Construct global sparsity pattern from local description

         This function can be called by the Matrix to get the sparsity pattern. 
         Assumes that the pattern is initially empty.
      */
      template<typename P>
      void fill_pattern (P& globalpattern) const
      {
        // map each cell to unique id
        MultiGeomUniqueIDMapper<GV> cell_mapper(gfsu.gridview());

        const bool has_subtriangulation = ST::hasSubTriangulation;
        
        for (ElementIterator it = gfsu.gridview().template begin<0>();
             it!=gfsu.gridview().template end<0>(); ++it)
          {
 			// bind local function spaces to element
			lfsu.bind(*it);
			lfsv.bind(*it);

            // compute unique id
            typename GV::IndexSet::IndexType id = cell_mapper.map(*it);

            // get local pattern of operator
            LocalSparsityPattern localpattern;
            LocalAssemblerCallSwitch<LA,LA::doPatternVolume>::
              pattern_volume(la,lfsu,lfsv,localpattern);


            if(has_subtriangulation){
              // translate local to global indices and add to global pattern
              for (size_t k=0; k<localpattern.size(); ++k)
                add_entry(globalpattern,
                          lfsv.globalIndex(localpattern[k].i()),
                          lfsu.globalIndex(localpattern[k].j())
                          );

              // reset local sparsity pattern
              localpattern = LocalSparsityPattern();
            }

            // skeleton and boundary pattern
            if(LA::doPatternSkeleton || LA::doPatternBoundary) {

              sub_triangulation.create(*it);
              SubIntersectionIterator iit = sub_triangulation.ibegin();
              SubIntersectionIterator eiit = sub_triangulation.iend();
              for(; iit!=eiit; ++iit)
              {
                // skeleton term
                if(iit->neighbor() && LA::doPatternSkeleton) {

                  // Recompute inside id
                  typename GV::IndexSet::IndexType ids = has_subtriangulation ? 
                    cell_mapper.map(*(iit->insideHostEntity())) : id;

                  // compute unique id
                  typename GV::IndexSet::IndexType idn = cell_mapper.map(*(iit->outside()));

                  // Visit face if id is bigger
                  bool visit_face = ids>idn;
                  // or interior is a ghost
                  visit_face |= (nonoverlapping_mode && (iit->inside())->partitionType()!=Dune::InteriorEntity);
                  // or local operator wishes to visit both sides
                  visit_face |= LA::doSkeletonTwoSided;
                  // or this is a domain interface
                  visit_face |= iit->insideDomainIndex() != iit->outsideDomainIndex();

                  if (visit_face){

                    // bind local function spaces to self element
                    if(has_subtriangulation){
                      lfsu.bind(*(iit->insideHostEntity()));
                      lfsv.bind(*(iit->insideHostEntity()));
                    }

                    // bind local function spaces to neighbor element
                    lfsun.bind(*(iit->outside()));
                    lfsvn.bind(*(iit->outside()));

                    // get pattern
                    LocalSparsityPattern localpattern_sn, localpattern_ns;
                    LocalAssemblerCallSwitch<LA,LA::doPatternSkeleton>::
                      pattern_skeleton(la,lfsu,lfsv,lfsun,lfsvn,
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
                }

                // boundary term
                if(iit->boundary()){

                  if(has_subtriangulation){
                    // bind local function spaces to self element
                    lfsu.bind(*(iit->insideHostEntity()));
                    lfsv.bind(*(iit->insideHostEntity()));
                  }

                  LocalAssemblerCallSwitch<LA,LA::doPatternBoundary>::
                    pattern_boundary(la,lfsu,lfsv,localpattern);

                  if(has_subtriangulation){
                    // translate local to global indices and add to global pattern
                    for (size_t k=0; k<localpattern.size(); ++k)
                      add_entry(globalpattern,
                                lfsv.globalIndex(localpattern[k].i()),
                                lfsu.globalIndex(localpattern[k].j())
                                );

                    // reset local sparsity pattern
                    localpattern = LocalSparsityPattern();
                  }

                } // boundary
			  } // iit
            } // if(LA::doPatternSkeleton || LA::doPatternBoundary)

            LocalAssemblerCallSwitch<LA,LA::doPatternVolumePostSkeleton>::
              pattern_volume_post_skeleton(la,lfsu,lfsv,localpattern);

            // translate local to global indices and add to global pattern
            for (size_t k=0; k<localpattern.size(); ++k)
              add_entry(globalpattern,
                        lfsv.globalIndex(localpattern[k].i()),
                        lfsu.globalIndex(localpattern[k].j())
                        );
          }
      }


	  //! generic evaluation of residual
      /**
       * \param r residual (needs to be cleared before this method is called)
       */
	  template<typename X, typename R> 
	  void residual (const X& x, R& r) const
	  {
        const bool has_subtriangulation = ST::hasSubTriangulation;

        // map each cell to unique id
        MultiGeomUniqueIDMapper<GV> cell_mapper(gfsu.gridview());

        // allocate local data container
        LocalVector<typename X::ElementType, TrialSpaceTag> xl;
        LocalVector<typename R::ElementType, TestSpaceTag> rl;
        LocalVector<typename X::ElementType, TrialSpaceTag> xn;
        LocalVector<typename R::ElementType, TestSpaceTag> rn;

		// traverse grid view
		for (ElementIterator it = gfsu.gridview().template begin<0>();
			 it!=gfsu.gridview().template end<0>(); ++it)
		  {
            // compute unique id
            typename GV::IndexSet::IndexType id = cell_mapper.map(*it);

            // skip ghost and overlap
            if (nonoverlapping_mode && it->partitionType()!=Dune::InteriorEntity)
              continue; 

            sub_triangulation.create(*it);
            SubEntityIterator sit = sub_triangulation.begin();
            SubEntityIterator esit = sub_triangulation.end();
            for(; sit != esit; ++sit){
        
              // bind local function spaces to element
              lfsu.bind( (sit->hostEntity()) );
              lfsv.bind( (sit->hostEntity()) );

              // allocate local data container
              xl.resize(lfsu.size());
              rl.assign(lfsv.size(),0.0);

              // read coefficents
              lfsu.vread(x,xl);

              // Rebind finite elements to the sub entity
              typename ST::BindEntityPartVisitor rebind_visitor(*sit);
              TypeTree::applyToTree(lfsu,rebind_visitor);
              TypeTree::applyToTree(lfsv,rebind_visitor);

              // volume evaluation
              LocalAssemblerCallSwitch<LA,LA::doAlphaVolume>::
                alpha_volume(la,*sit,lfsu,xl,lfsv,rl);
              LocalAssemblerCallSwitch<LA,LA::doLambdaVolume>::
                lambda_volume(la,*sit,lfsv,rl);

              // Accumulate local residuals for each sub entity if we
              // have a sub triangulation.
              if(has_subtriangulation)
                lfsv.vadd(rl,r);
            }

            // skip if no intersection iterator is needed
            if (LA::doAlphaSkeleton||LA::doAlphaBoundary||LA::doLambdaSkeleton||LA::doLambdaBoundary)
              {
                // traverse intersections
                SubIntersectionIterator iit = sub_triangulation.ibegin();
                SubIntersectionIterator eiit = sub_triangulation.iend();
                for (; iit!=eiit; ++iit)
                  {

                    typename ST::BindInsideIntersectionPartVisitor in_rebind_visitor(*iit);
                    typename ST::BindOutsideIntersectionPartVisitor out_rebind_visitor(*iit);
                    typename ST::BindIntersectionPartVisitor rebind_intersection_visitor(*iit);

                    if(has_subtriangulation){
                      // Rebind local function spaces and read coefficents
                      lfsu.bind(*( iit->insideHostEntity() ));
                      lfsv.bind(*( iit->insideHostEntity() ));

                      // allocate local data container
                      xl.resize(lfsu.size());
                      rl.assign(lfsv.size(),0.0);

                      // read coefficents
                      lfsu.vread(x,xl);
                    }

                    // skeleton term
                    if (iit->neighbor() && (LA::doAlphaSkeleton||LA::doLambdaSkeleton) )
                      {
                        // compute unique id for neighbor
                        const typename GV::IndexSet::IndexType ids = has_subtriangulation ? 
                          cell_mapper.map(*(iit->insideHostEntity())) : id;
                        const typename GV::IndexSet::IndexType idn = cell_mapper.map(*(iit->outside()));

                        // Visit face if id is bigger
                        bool visit_face = ids > idn || LA::doSkeletonTwoSided;
                        // or interior is a ghost
                        visit_face |= (nonoverlapping_mode && 
                                       (iit->inside())->partitionType()!=Dune::InteriorEntity);
                        // or this is a domain interface
                        visit_face |= iit->insideDomainIndex() != iit->outsideDomainIndex();
                          
                        // unique vist of intersection
                        if (visit_face)
                          {
                            // bind local function spaces to neighbor element
                            lfsun.bind(*(iit->outside()));
                            lfsvn.bind(*(iit->outside()));
                            
                            // allocate local data container
                            xn.resize(lfsun.size());
                            rn.assign(lfsvn.size(),0.0);
                            
                            // read coefficents
                            lfsun.vread(x,xn);
                            
                            // Prepare the local functions for global evaluation 
                            TypeTree::applyToTree(lfsu,in_rebind_visitor);
                            TypeTree::applyToTree(lfsun,out_rebind_visitor);
                            TypeTree::applyToTree(lfsv,in_rebind_visitor);
                            TypeTree::applyToTree(lfsvn,out_rebind_visitor);

                            // ST::BindInsideSubIntersection::rebind(lfsu,*iit);
                            // ST::BindOutsideSubIntersection::rebind(lfsun,*iit);

                            // ST::BindInsideSubIntersection::rebind(lfsv,*iit);
                            // ST::BindOutsideSubIntersection::rebind(lfsvn,*iit);

                            // skeleton evaluation
                            LocalAssemblerCallSwitch<LA,LA::doAlphaSkeleton>::
                              alpha_skeleton(la,*iit,lfsu,xl,lfsv,lfsun,xn,lfsvn,rl,rn);
                            LocalAssemblerCallSwitch<LA,LA::doLambdaSkeleton>::
                              lambda_skeleton(la,*iit,lfsv,lfsvn,rl,rn);

                            // accumulate result (note: r needs to be cleared outside)
                            lfsvn.vadd(rn,r);
                          }
                      }

                    // If domain boundary
                    if(iit->boundary()){

                      // Prepare the local functions for global evaluation 
                      // ST::BindSubIntersection::rebind(lfsu,*iit);
                      // ST::BindSubIntersection::rebind(lfsv,*iit);
                      TypeTree::applyToTree(lfsu,rebind_intersection_visitor);
                      TypeTree::applyToTree(lfsv,rebind_intersection_visitor);
                              
                      LocalAssemblerCallSwitch<LA,LA::doAlphaBoundary>::
                        alpha_boundary(la,*iit,lfsu,xl,lfsv,rl);
                      LocalAssemblerCallSwitch<LA,LA::doLambdaBoundary>::
                        lambda_boundary(la,*iit,lfsv,rl);
                    }

                    // Accumulate local residuals for each sub entity if we
                    // have a sub triangulation.
                    if(has_subtriangulation)
                      lfsv.vadd(rl,r);

                  } // iit
              } // do skeleton

            
            if(!has_subtriangulation){
              LocalAssemblerCallSwitch<LA,LA::doAlphaVolumePostSkeleton>::
                alpha_volume_post_skeleton(la,ElementGeometry<Element>(*it),lfsu,xl,lfsv,rl);
              LocalAssemblerCallSwitch<LA,LA::doLambdaVolumePostSkeleton>::
                lambda_volume_post_skeleton(la,ElementGeometry<Element>(*it),lfsv,rl);

              // accumulate result (note: r needs to be cleared outside)
              lfsv.vadd(rl,r);
            }

		  } // it

		// set residual to zero on constrained dofs
		Dune::PDELab::constrain_residual(*pconstraintsv,r);
	  }

      //! \brief Assemble constant part of residual for operators with purely
      //!        linear alpha_*() methods
      /**
       * The result of this method is identical to calling residual() with
       * parameter \c x initialized to 0.  However, this method works only
       * when the local operator has purely linear alpha_volume(),
       * alpha_skeleton(), alpha_boundary() and alpha_volume_post_skeleton()
       * methods.  Note that this is a stronger demand than requiring the
       * operator to be affine -- affine operators may still have an affine
       * shift in their alpha_*() methods.  For this method to work any affine
       * shift must be implemented in the lambda_*() methods.
       *
       * \note This method is meaningless in the context of non-linear
       *       operators.
       *
       * \param r residual (needs to be cleared before this method is called)
       */
      template<typename R>
      void zero_residual(R& r) const 
	  {
        const bool has_subtriangulation = ST::hasSubTriangulation;

        // map each cell to unique id
        MultiGeomUniqueIDMapper<GV> cell_mapper(gfsu.gridview());

        // allocate local data container
        LocalVector<typename R::ElementType, TestSpaceTag> rl;
        LocalVector<typename R::ElementType, TestSpaceTag> rn;

		// traverse grid view
		for (ElementIterator it = gfsu.gridview().template begin<0>();
			 it!=gfsu.gridview().template end<0>(); ++it)
		  {
            // compute unique id
            typename GV::IndexSet::IndexType id = cell_mapper.map(*it);


            // skip ghost and overlap
            if (nonoverlapping_mode && it->partitionType()!=Dune::InteriorEntity)
              continue; 

            sub_triangulation.create(*it);
            SubEntityIterator sit = sub_triangulation.begin();
            SubEntityIterator esit = sub_triangulation.end();
            for(; sit != esit; ++sit){
        
              // bind local function spaces to element
              lfsv.bind( (sit->hostEntity()) );

              // allocate local data container
              rl.assign(lfsv.size(),0.0);

              // Rebind finite elements to the sub entity
              typename ST::BindEntityPartVisitor rebind_visitor(*sit);
              TypeTree::applyToTree(lfsv,rebind_visitor);
              // ST::BindSubEntity::rebind(lfsv,*sit);

              // volume evaluation
              LocalAssemblerCallSwitch<LA,LA::doLambdaVolume>::
                lambda_volume(la,*sit,lfsv,rl);

              // Accumulate local residuals for each sub entity if we
              // have a sub triangulation.
              if(has_subtriangulation)
                lfsv.vadd(rl,r);
            }

            // skip if no intersection iterator is needed
            if (LA::doLambdaSkeleton||LA::doLambdaBoundary)
              {
                // traverse intersections
                SubIntersectionIterator iit = sub_triangulation.ibegin();
                SubIntersectionIterator eiit = sub_triangulation.iend();
                for (; iit!=eiit; ++iit)
                  {

                    typename ST::BindInsideIntersectionPartVisitor in_rebind_visitor(*iit);
                    typename ST::BindOutsideIntersectionPartVisitor out_rebind_visitor(*iit);
                    typename ST::BindIntersectionPartVisitor rebind_intersection_visitor(*iit);

                    if(has_subtriangulation){
                      // Rebind local function spaces and read coefficents
                      lfsv.bind(*( iit->insideHostEntity() ));

                      // allocate local data container
                      rl.assign(lfsv.size(),0.0);
                    }

                    // skeleton term
                    if (iit->neighbor() && (LA::doLambdaSkeleton) )
                      {
                        // compute unique id for neighbor
                        const typename GV::IndexSet::IndexType ids = has_subtriangulation ? 
                          cell_mapper.map(*(iit->insideHostEntity())) : id;
                        const typename GV::IndexSet::IndexType idn = cell_mapper.map(*(iit->outside()));

                        // Visit face if id is bigger
                        bool visit_face = ids > idn || LA::doSkeletonTwoSided;
                        // or interior is a ghost
                        visit_face |= (nonoverlapping_mode && 
                                       (iit->inside())->partitionType()!=Dune::InteriorEntity);
                        // or this is a domain interface
                        visit_face |= iit->insideDomainIndex() != iit->outsideDomainIndex();
                          
                        // unique vist of intersection
                        if (visit_face)
                          {
                            // bind local function spaces to neighbor element
                            lfsvn.bind(*(iit->outside()));
                            
                            // allocate local data container
                            rn.assign(lfsvn.size(),0.0);
                            
                            // Prepare the local functions for global evaluation 
                            TypeTree::applyToTree(lfsv,in_rebind_visitor);
                            TypeTree::applyToTree(lfsvn,out_rebind_visitor);
                            // ST::BindInsideSubIntersection::rebind(lfsv,*iit);
                            // ST::BindOutsideSubIntersection::rebind(lfsvn,*iit);

                            // skeleton evaluation
                            LocalAssemblerCallSwitch<LA,LA::doLambdaSkeleton>::
                              lambda_skeleton(la,*iit,lfsv,lfsvn,rl,rn);

                            // accumulate result (note: r needs to be cleared outside)
                            lfsvn.vadd(rn,r);
                          }
                      }

                    // If domain boundary
                    if(iit->boundary()){

                      // Prepare the local functions for global evaluation 
                      //ST::BindSubIntersection::rebind(lfsv,*iit);
                      TypeTree::applyToTree(lfsv,rebind_intersection_visitor);
                              
                      LocalAssemblerCallSwitch<LA,LA::doLambdaBoundary>::
                        lambda_boundary(la,*iit,lfsv,rl);
                    }

                    // Accumulate local residuals for each sub entity if we
                    // have a sub triangulation.
                    if(has_subtriangulation)
                      lfsv.vadd(rl,r);

                  } // iit
              } // do skeleton

            
            if(!has_subtriangulation){
              LocalAssemblerCallSwitch<LA,LA::doLambdaVolumePostSkeleton>::
                lambda_volume_post_skeleton(la,ElementGeometry<Element>(*it),lfsv,rl);

              // accumulate result (note: r needs to be cleared outside)
              lfsv.vadd(rl,r);
            }

		  } // it

		// set residual to zero on constrained dofs
		Dune::PDELab::constrain_residual(*pconstraintsv,r);
	  }

	  //! generic application of Jacobian
	  template<typename X, typename Y> 
	  void jacobian_apply (X& x, Y& y) const
	  {
        const bool has_subtriangulation = ST::hasSubTriangulation;

        // map each cell to unique id
        MultiGeomUniqueIDMapper<GV> cell_mapper(gfsu.gridview());

        // allocate local data container
        LocalVector<typename X::ElementType, TrialSpaceTag> xl;
        LocalVector<typename Y::ElementType, TestSpaceTag> yl;
        LocalVector<typename X::ElementType, TrialSpaceTag> xn;
        LocalVector<typename Y::ElementType, TestSpaceTag> yn;

		// traverse grid view
		for (ElementIterator it = gfsu.gridview().template begin<0>();
			 it!=gfsu.gridview().template end<0>(); ++it)
		  {
            // compute unique id
            typename GV::IndexSet::IndexType id = cell_mapper.map(*it);

            // skip ghost and overlap
            if (nonoverlapping_mode && it->partitionType()!=Dune::InteriorEntity)
              continue; 

            sub_triangulation.create(*it);
            SubEntityIterator sit = sub_triangulation.begin();
            SubEntityIterator esit = sub_triangulation.end();
            for(; sit != esit; ++sit){
        
              // bind local function spaces to element
              lfsu.bind( (sit->hostEntity()) );
              lfsv.bind( (sit->hostEntity()) );

              // allocate local data container
              xl.resize(lfsu.size());
              yl.assign(lfsv.size(),0.0);

              // read coefficents
              lfsu.vread(x,xl);

              // Rebind finite elements to the sub entity
              typename ST::BindEntityPartVisitor rebind_visitor(*sit);
              TypeTree::applyToTree(lfsu,rebind_visitor);
              TypeTree::applyToTree(lfsv,rebind_visitor);
              // ST::BindSubEntity::rebind(lfsu,*sit);
              // ST::BindSubEntity::rebind(lfsv,*sit);

              // volume evaluation
              LocalAssemblerCallSwitch<LA,LA::doAlphaVolume>::
                jacobian_apply_volume(la,*sit,lfsu,xl,lfsv,yl);

              // Accumulate local residuals for each sub entity if we
              // have a sub triangulation.
              if(has_subtriangulation)
                lfsv.vadd(yl,y);
            }

            // skip if no intersection iterator is needed
            if (LA::doAlphaSkeleton||LA::doAlphaBoundary)
              {
                // traverse intersections
                SubIntersectionIterator iit = sub_triangulation.ibegin();
                SubIntersectionIterator eiit = sub_triangulation.iend();
                for (; iit!=eiit; ++iit)
                  {

                    typename ST::BindInsideIntersectionPartVisitor in_rebind_visitor(*iit);
                    typename ST::BindOutsideIntersectionPartVisitor out_rebind_visitor(*iit);
                    typename ST::BindIntersectionPartVisitor rebind_intersection_visitor(*iit);

                    if(has_subtriangulation){
                      // Rebind local function spaces and read coefficents
                      lfsu.bind(*( iit->insideHostEntity() ));
                      lfsv.bind(*( iit->insideHostEntity() ));

                      // allocate local data container
                      xl.resize(lfsu.size());
                      yl.assign(lfsv.size(),0.0);

                      // read coefficents
                      lfsu.vread(x,xl);
                    }

                    // skeleton term
                    if (iit->neighbor() && (LA::doAlphaSkeleton||LA::doLambdaSkeleton) )
                      {
                        // compute unique id for neighbor
                        const typename GV::IndexSet::IndexType ids = has_subtriangulation ? 
                          cell_mapper.map(*(iit->insideHostEntity())) : id;
                        const typename GV::IndexSet::IndexType idn = cell_mapper.map(*(iit->outside()));

                        // Visit face if id is bigger
                        bool visit_face = ids > idn || LA::doSkeletonTwoSided;
                        // or interior is a ghost
                        visit_face |= (nonoverlapping_mode && 
                                       (iit->inside())->partitionType()!=Dune::InteriorEntity);
                        // or this is a domain interface
                        visit_face |= iit->insideDomainIndex() != iit->outsideDomainIndex();
                          
                        // unique vist of intersection
                        if (visit_face)
                          {
                            // bind local function spaces to neighbor element
                            lfsun.bind(*(iit->outside()));
                            lfsvn.bind(*(iit->outside()));
                            
                            // allocate local data container
                            xn.resize(lfsun.size());
                            yn.assign(lfsvn.size(),0.0);
                            
                            // read coefficents
                            lfsun.vread(x,xn);
                            
                            // Prepare the local functions for global evaluation 
                            TypeTree::applyToTree(lfsu,in_rebind_visitor);
                            TypeTree::applyToTree(lfsun,out_rebind_visitor);
                            TypeTree::applyToTree(lfsv,in_rebind_visitor);
                            TypeTree::applyToTree(lfsvn,out_rebind_visitor);
                            // ST::BindInsideSubIntersection::rebind(lfsu,*iit);
                            // ST::BindOutsideSubIntersection::rebind(lfsun,*iit);
                            // ST::BindInsideSubIntersection::rebind(lfsv,*iit);
                            // ST::BindOutsideSubIntersection::rebind(lfsvn,*iit);

                            // skeleton evaluation
                            LocalAssemblerCallSwitch<LA,LA::doAlphaSkeleton>::
                              jacobian_apply_skeleton(la,*iit,lfsu,xl,lfsv,lfsun,xn,lfsvn,yl,yn);

                            // accumulate result (note: r needs to be cleared outside)
                            lfsvn.vadd(yn,y);
                          }
                      }

                    // If domain boundary
                    if(iit->boundary()){

                      // Prepare the local functions for global evaluation 
                      // ST::BindSubIntersection::rebind(lfsu,*iit);
                      // ST::BindSubIntersection::rebind(lfsv,*iit);
                      TypeTree::applyToTree(lfsu,rebind_intersection_visitor);
                      TypeTree::applyToTree(lfsv,rebind_intersection_visitor);

                      LocalAssemblerCallSwitch<LA,LA::doAlphaBoundary>::
                        jacobian_apply_boundary(la,*iit,lfsu,xl,lfsv,yl);                              
                    }

                    // Accumulate local residuals for each sub entity if we
                    // have a sub triangulation.
                    if(has_subtriangulation)
                      lfsv.vadd(yl,y);

                  } // iit
              } // do skeleton

            if(!has_subtriangulation){
              LocalAssemblerCallSwitch<LA,LA::doAlphaVolumePostSkeleton>::
                jacobian_apply_volume_post_skeleton(la,ElementGeometry<Element>(*it),lfsu,xl,lfsv,yl);

              // accumulate result (note: r needs to be cleared outside)
              lfsv.vadd(yl,y);
            }

		  }

		// set residual to zero on constrained dofs
		Dune::PDELab::copy_constrained_dofs(*pconstraintsu,x,y);
	  }

	  //! generic assembly of Jacobian
      /**
       * \param x Where (in the space spanned by the dofs) to evaluate the Jacobian
       * \param a Jacobian (needs to be cleared before passed to this method)
       */
	  template<typename X, typename A> 
	  void jacobian (const X& x, A& a) const
	  {
        const bool has_subtriangulation = ST::hasSubTriangulation;

        // map each cell to unique id
        MultiGeomUniqueIDMapper<GV> cell_mapper(gfsu.gridview());

        // allocate local data container
        LocalVector<typename X::ElementType, TrialSpaceTag> xn;
        LocalVector<typename X::ElementType, TrialSpaceTag> xl;
        LocalMatrix<typename A::ElementType> al;
        LocalMatrix<typename A::ElementType> al_sn;
        LocalMatrix<typename A::ElementType> al_ns;
        LocalMatrix<typename A::ElementType> al_nn;
        
		// traverse grid view
		for (ElementIterator it = gfsu.gridview().template begin<0>();
			 it!=gfsu.gridview().template end<0>(); ++it)
		  {

            // compute unique id
            const typename GV::IndexSet::IndexType id = cell_mapper.map(*it);

            // skip ghost and overlap
            if (nonoverlapping_mode && it->partitionType()!=Dune::InteriorEntity)
              continue; 

            sub_triangulation.create(*it);
            SubEntityIterator sit = sub_triangulation.begin();
            SubEntityIterator esit = sub_triangulation.end();
            for(; sit != esit; ++sit){
        
              // bind local function spaces to element
              lfsu.bind( (sit->hostEntity()) );
              lfsv.bind( (sit->hostEntity()) );

              // allocate local data container
              xl.resize(lfsu.size());
              al.assign(lfsv.size() ,lfsu.size(),0.0);

              // read coefficents
              lfsu.vread(x,xl);

              // Rebind finite elements to the sub entity
              typename ST::BindEntityPartVisitor rebind_visitor(*sit);
              TypeTree::applyToTree(lfsu,rebind_visitor);
              TypeTree::applyToTree(lfsv,rebind_visitor);
              // ST::BindSubEntity::rebind(lfsu,*sit);
              // ST::BindSubEntity::rebind(lfsv,*sit);

              // volume evaluation
              LocalAssemblerCallSwitch<LA,LA::doAlphaVolume>::
                jacobian_volume(la,*sit,lfsu,xl,lfsv,al);

              // Accumulate local residuals for each sub entity if we
              // have a sub triangulation.
              if(has_subtriangulation)
                etadd(lfsv,lfsu,al,a);
            }

            // skip if no intersection iterator is needed
            if (LA::doAlphaSkeleton||LA::doAlphaBoundary)
              {
                // traverse intersections
                SubIntersectionIterator iit = sub_triangulation.ibegin();
                SubIntersectionIterator eiit = sub_triangulation.iend();
                for (; iit!=eiit; ++iit)
                  {

                    typename ST::BindInsideIntersectionPartVisitor in_rebind_visitor(*iit);
                    typename ST::BindOutsideIntersectionPartVisitor out_rebind_visitor(*iit);
                    typename ST::BindIntersectionPartVisitor rebind_intersection_visitor(*iit);

                    if(has_subtriangulation){
                      // Rebind local function spaces and read coefficents
                      lfsu.bind(*( iit->insideHostEntity() ));
                      lfsv.bind(*( iit->insideHostEntity() ));

                      // allocate local data container
                      xl.resize(lfsu.size());
                      al.assign(lfsv.size() ,lfsu.size(),0.0);

                      // read coefficents
                      lfsu.vread(x,xl);
                    }

                    // skeleton term
                    if (iit->neighbor() && (LA::doAlphaSkeleton||LA::doLambdaSkeleton) )
                      {
                        // compute unique id for neighbor
                        const typename GV::IndexSet::IndexType ids = has_subtriangulation ? 
                          cell_mapper.map(*(iit->insideHostEntity())) : id;
                        const typename GV::IndexSet::IndexType idn = cell_mapper.map(*(iit->outside()));

                        // Visit face if id is bigger
                        bool visit_face = ids > idn || LA::doSkeletonTwoSided;
                        // or interior is a ghost
                        visit_face |= (nonoverlapping_mode && 
                                       (iit->inside())->partitionType()!=Dune::InteriorEntity);
                        // or this is a domain interface
                        visit_face |= iit->insideDomainIndex() != iit->outsideDomainIndex();
                          
                        // unique vist of intersection
                        if (visit_face)
                          {
                            // bind local function spaces to neighbor element
                            lfsun.bind(*(iit->outside()));
                            lfsvn.bind(*(iit->outside()));
                            
                            // allocate local data container
                            xn.resize(lfsun.size());
                            al_sn.assign(lfsv.size() ,lfsun.size(),0.0);
                            al_ns.assign(lfsvn.size(),lfsu.size() ,0.0);
                            al_nn.assign(lfsvn.size(),lfsun.size(),0.0);
                            
                            // read coefficents
                            lfsun.vread(x,xn);
                            
                            // Prepare the local functions for global evaluation 
                            TypeTree::applyToTree(lfsu,in_rebind_visitor);
                            TypeTree::applyToTree(lfsun,out_rebind_visitor);
                            TypeTree::applyToTree(lfsv,in_rebind_visitor);
                            TypeTree::applyToTree(lfsvn,out_rebind_visitor);

                            // ST::BindInsideSubIntersection::rebind(lfsu,*iit);
                            // ST::BindOutsideSubIntersection::rebind(lfsun,*iit);

                            // ST::BindInsideSubIntersection::rebind(lfsv,*iit);
                            // ST::BindOutsideSubIntersection::rebind(lfsvn,*iit);

                            // skeleton evaluation
                            LocalAssemblerCallSwitch<LA,LA::doAlphaSkeleton>::
                              jacobian_skeleton(la,*iit,lfsu,xl,lfsv,lfsun,xn,lfsvn,al,al_sn,al_ns,al_nn);

                            // accumulate result
                            etadd(lfsv,lfsun,al_sn,a);
                            etadd(lfsvn,lfsu,al_ns,a);
                            etadd(lfsvn,lfsun,al_nn,a);
                          }
                      }

                    // If domain boundary
                    if(iit->boundary()){

                      // Prepare the local functions for global evaluation 
                      TypeTree::applyToTree(lfsu,rebind_intersection_visitor);
                      TypeTree::applyToTree(lfsv,rebind_intersection_visitor);

                      // ST::BindSubIntersection::rebind(lfsu,*iit);
                      // ST::BindSubIntersection::rebind(lfsv,*iit);

                      LocalAssemblerCallSwitch<LA,LA::doAlphaBoundary>::
                        jacobian_boundary(la,*iit,lfsu,xl,lfsv,al);                              
                    }

                    // Accumulate local residuals for each sub entity if we
                    // have a sub triangulation.
                    if(has_subtriangulation)
                      etadd(lfsv,lfsu,al,a);

                  } // iit
              } // do skeleton

            if(!has_subtriangulation){
              LocalAssemblerCallSwitch<LA,LA::doAlphaVolumePostSkeleton>::
                jacobian_volume_post_skeleton(la,ElementGeometry<Element>(*it),lfsu,xl,lfsv,al);

              // accumulate result (note: a needs to be cleared outside)
              etadd(lfsv,lfsu,al,a);
            }
		  }

		
         typedef typename CV::const_iterator global_row_iterator;	  
         for (global_row_iterator cit=pconstraintsv->begin(); cit!=pconstraintsv->end(); ++cit)
           set_trivial_row(cit->first,cit->second,a);
 	  }

	private:
	  const LA& la;
      const SubTriangulation sub_triangulation;
      using Base::gfsu;
      using Base::gfsv;
      using Base::pconstraintsu;
      using Base::pconstraintsv;
      using Base::lfsu;
      using Base::lfsv;
      using Base::lfsun;
      using Base::lfsvn;
	};

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
