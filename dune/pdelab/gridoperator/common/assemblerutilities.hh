// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_ASSEMBLERUTILITIES_HH
#define DUNE_PDELAB_ASSEMBLERUTILITIES_HH

#include <dune/pdelab/constraints/constraintstransformation.hh>
#include <dune/pdelab/gridoperatorspace/localmatrix.hh>

namespace Dune{
  namespace PDELab{

    /** Traits of the local assembler

        \tparam GO The grid operator

    */
    template<typename GO>
    struct LocalAssemblerTraits
    {

      //! The trial grid function space.
      typedef typename GO::Traits::TrialGridFunctionSpace TrialGridFunctionSpace;

      //! The test grid function space.
      typedef typename GO::Traits::TestGridFunctionSpace TestGridFunctionSpace;


      //! The type of the trial grid function space constraints.
      typedef typename GO::Traits::TrialGridFunctionSpaceConstraints TrialGridFunctionSpaceConstraints;

      //! The type of the test grid function space constraints.
      typedef typename GO::Traits::TestGridFunctionSpaceConstraints TestGridFunctionSpaceConstraints;


      //! The matrix backend of the grid operator.
      typedef typename GO::Traits::MatrixBackend MatrixBackend;


      //! The field type of the domain (solution).
      typedef typename GO::Traits::DomainField DomainField;

      //! The type of the domain (solution).
      typedef typename GO::Traits::Domain Solution;


      //! The field type of the range (residual).
      typedef typename GO::Traits::RangeField RangeField;

      //! The type of the range (residual).
      typedef typename GO::Traits::Range Residual;


      //! The field type of the jacobian.
      typedef typename GO::Traits::JacobianField JacobianField;

      //! The type of the jacobian.
      typedef typename GO::Traits::Jacobian Jacobian;

      //! The matrix pattern
      typedef typename Jacobian::Pattern MatrixPattern;


    };

    //! Translation helper from intersection method return values to intersection type.
    /**
     * This struct can be used to query an intersection for its type in a convenient way.
     * Classification occurs in accordance with the specification for the return values
     * of Intersection::neighbor() and Intersection::boundary() given in the documentation
     * of the Intersection class.
     */
    struct IntersectionType
    {

      enum Type {
        processor = 0, //!< processor boundary intersection (neighbor() == false && boundary() == false)
        skeleton = 1,  //!< skeleton intersection (neighbor() == true && boundary() == false)
        boundary = 2,  //!< domain boundary intersection (neighbor() == false && boundary() == true)
        periodic = 3   //!< periodic boundary intersection (neighbor() == true && boundary() == true)
      };

      //! Returns the type of the intersection.
      template<typename Intersection>
      static Type get(const Intersection& is)
      {
        return static_cast<Type>(1*is.neighbor() + 2*is.boundary());
      }

    };

    /**
       \brief Base class for local assembler

       This class provides some generic behavior required for local
       assembler implementations. This includes the access of the global
       vectors and matrices via local indices and local function spaces
       with regard to the constraint mappings.

       \tparam B    The matrix backend
       \tparam CU   Constraints maps for the individual dofs (trial space)
       \tparam CV   Constraints maps for the individual dofs (test space)
    */
    template<typename B,
             typename CU=EmptyTransformation,
             typename CV=EmptyTransformation>
    class LocalAssemblerBase{
    public:

      typedef typename B::size_type SizeType;

      //! construct AssemblerSpace
      LocalAssemblerBase ()
        : pconstraintsu(&emptyconstraintsu), pconstraintsv(&emptyconstraintsv)
      {}

      //! construct AssemblerSpace, with constraints
      LocalAssemblerBase (const CU& cu, const CV& cv)
        :pconstraintsu(&cu), pconstraintsv(&cv)
      {}

      //! get the constraints on the trial grid function space
      const CU& trialConstraints() const
      {
        return *pconstraintsu;
      }

      //! get the constraints on the test grid function space
      const CV& testConstraints() const
      {
        return *pconstraintsv;
      }


      /** \brief Transforms a vector \f$ \boldsymbol{x} \f$ from \f$
          V\f$ to \f$ V'\f$. If postrestrict == true then
          \f$\boldsymbol{R}^T_{\boldsymbol{\tilde U}', \boldsymbol{U}'}
          \boldsymbol{S}_{\boldsymbol{\tilde V}}\f$ is applied
          instead of the full transformation.  */
      template<typename X>
      void forwardtransform(X & x, const bool postrestrict = false)
      {
        typedef typename CV::const_iterator global_col_iterator;
        for (global_col_iterator cit=pconstraintsv->begin(); cit!=pconstraintsv->end(); ++cit){
          typedef typename global_col_iterator::value_type::first_type GlobalIndex;
          const GlobalIndex & contributor = cit->first;

          typedef typename global_col_iterator::value_type::second_type ContributedMap;
          typedef typename ContributedMap::const_iterator global_row_iterator;
          const ContributedMap & contributed = cit->second;
          global_row_iterator it  = contributed.begin();
          global_row_iterator eit = contributed.end();

          for(;it!=eit;++it)
            x[it->first] += it->second * x[contributor];
        }

        if(postrestrict)
          for (global_col_iterator cit=pconstraintsv->begin(); cit!=pconstraintsv->end(); ++cit)
            x[cit->first]=0.;
      }

      /** \brief Transforms a vector \f$ \boldsymbol{x} \f$ from \f$
          V'\f$ to \f$ V\f$. If prerestrict == true then
          \f$\boldsymbol{S}^T_{\boldsymbol{\tilde U}}\f$ is applied
          instead of the full transformation.  */
      template<typename X>
      void backtransform(X & x, const bool prerestrict = false)
      {
        typedef typename CV::const_iterator global_col_iterator;
        for (global_col_iterator cit=pconstraintsv->begin(); cit!=pconstraintsv->end(); ++cit){
          typedef typename global_col_iterator::value_type::first_type GlobalIndex;
          const GlobalIndex & contributor = cit->first;

          typedef typename global_col_iterator::value_type::second_type ContributedMap;
          typedef typename ContributedMap::const_iterator global_row_iterator;
          const ContributedMap & contributed = cit->second;
          global_row_iterator it  = contributed.begin();
          global_row_iterator eit = contributed.end();

          if(prerestrict)
            x[contributor] = 0.;

          for(;it!=eit;++it)
            x[contributor] += it->second * x[it->first];
        }
      }

    protected:

      /** \brief read local stiffness matrix for entity */
      template<typename GCView, typename T>
      void eread (const GCView& globalcontainer_view, LocalMatrix<T>& localcontainer) const
      {
        for (int i = 0; i < localcontainer.N(); ++i)
          for (int j = 0; j < localcontainer.M(); ++j)
            localcontainer(i,j) = globalcontainer_view(i,j);
      }

      /** \brief write local stiffness matrix for entity */
      template<typename T, typename GCView>
      void ewrite (const LocalMatrix<T>& localcontainer, GCView& globalcontainer_view) const
      {
        for (int i = 0; i < localcontainer.N(); ++i)
          for (int j = 0; j < localcontainer.M(); ++j)
            globalcontainer_view(i,j) = localcontainer(i,j);
      }

      /** \brief write local stiffness matrix for entity */
      template<typename T, typename GCView>
      void eadd (const LocalMatrix<T>& localcontainer, GCView& globalcontainer_view) const
      {
        for (int i = 0; i < localcontainer.N(); ++i)
          for (int j = 0; j < localcontainer.M(); ++j)
            globalcontainer_view.add(i,j,localcontainer(i,j));
      }


      /** \brief Add local matrix to global matrix,
          and apply Dirichlet constraints in a symmetric
          fashion. Apart from that, identical to etadd(). */
      template<typename LFSV, typename LFSU, typename T, typename GC>
      void etadd_symmetric (const LFSV& lfsv, const LFSU& lfsu, LocalMatrix<T>& localcontainer, GC& globalcontainer) const
      {
        const CU & cu = *pconstraintsu;

        typedef typename CU::const_iterator global_ucol_iterator;

        for (size_t j = 0; j < lfsu.size(); ++j)
          {
            global_ucol_iterator cuit = cu.find(lfsu.globalIndex(j));

            // If this column is not constrained or the constraint is not of
            // Dirichlet type, abort
            if (cuit == cu.end() || cuit->second.size() > 0)
              continue;

            // clear out the current column
            for (size_t i = 0; i < lfsv.size(); ++i)
              {
                // we do not need to update the residual, since the solution
                // (i.e. the correction) for the Dirichlet DOF is 0 by definition
                localcontainer(lfsv,i,lfsu,j) = 0.0;
              }
          }

        // hand off to standard etadd() method
        etadd(lfsv,lfsu,localcontainer,globalcontainer);
      }


      /** \brief Add local matrix \f$m\f$ to global Jacobian \f$J\f$
          and apply constraints transformation. Hence we perform: \f$
          \boldsymbol{J} := \boldsymbol{J} + \boldsymbol{S}_{
          \boldsymbol{\tilde V}} m \boldsymbol{S}^T_{
          \boldsymbol{\tilde U}} \f$*/
      template<typename LFSV, typename LFSU, typename T, typename GC>
      void etadd (const LFSV& lfsv, const LFSU& lfsu, const LocalMatrix<T>& localcontainer, GC& globalcontainer) const
      {

        typename B::template Accessor<LFSV,LFSU,T> accessor(globalcontainer,lfsv,lfsu);

        for (size_t i=0; i<lfsv.size(); i++)
          for (size_t j=0; j<lfsu.size(); j++){
            SizeType gi = lfsv.globalIndex(i);
            SizeType gj = lfsu.globalIndex(j);

            // Get global constraints containers for test and ansatz space
            const CV & cv = *pconstraintsv;
            const CU & cu = *pconstraintsu;

            typedef typename CV::const_iterator global_vcol_iterator;
            typedef typename global_vcol_iterator::value_type::second_type global_vrow_type;
            typedef typename global_vrow_type::const_iterator global_vrow_iterator;

            typedef typename CU::const_iterator global_ucol_iterator;
            typedef typename global_ucol_iterator::value_type::second_type global_urow_type;
            typedef typename global_urow_type::const_iterator global_urow_iterator;

            // Check whether the global indices are constrained indices
            global_vcol_iterator gvcit = cv.find(gi);
            global_ucol_iterator gucit = cu.find(gj);

            // Set constrained_v true if gi is constrained dof
            bool constrained_v(false);
            global_vrow_iterator gvrit;
            if(gvcit!=cv.end()){
              gvrit = gvcit->second.begin();
              constrained_v = true;
            }

            T vf = 1;
            do{
              // if gi is index of constrained dof
              if(constrained_v){

                if(gvrit == gvcit->second.end())
                  break;

                // otherwise set gi to an index to a contributed dof
                // and set vf to the contribution weight
                gi = gvrit->first;
                vf = gvrit->second;
              }

              // Set constrained_u true if gj is constrained dof
              bool constrained_u(false);
              global_urow_iterator gurit;
              if(gucit!=cu.end()){
                gurit = gucit->second.begin();
                constrained_u = true;
                if(gurit == gucit->second.end()){
                  T t = localcontainer(lfsv,i,lfsu,j) * vf;
                  if(t != 0.0)                 // entry might not be present in the matrix
                    accessor.add(i,j,t);
                }
              }

              T uf = 1;
              do{
                // if gj is index of constrained dof
                if(constrained_u){

                  if(gurit == gucit->second.end())
                    break;

                  // otherwise set gj to an index to a contributed dof
                  // and set uf to the contribution weight
                  gj = gurit->first;
                  uf = gurit->second;
                }

                // add weighted local entry to global matrix
                T t = localcontainer(lfsv,i,lfsu,j) * uf * vf;
                if (t != 0.0)                 // entry might not be present in the matrix
                  accessor.add(i,j,t);

                if(constrained_u && gurit != gucit->second.end())
                  ++gurit;
                else
                  break;

              }while(true);

              if(constrained_v && gvrit != gvcit->second.end())
                ++gvrit;
              else
                break;

            }while(true);

          }
      }


      /** \brief Adding matrix entry to pattern with respect to the
          constraints contributions. This assembles the entries addressed
          by etadd(..). See the documentation there for more information
          about the matrix pattern. */
      template<typename P, typename LFSVIndices, typename LFSUIndices, typename Index>
      void add_entry(P & globalpattern,
                     const LFSVIndices& lfsv_indices, Index i,
                     const LFSUIndices& lfsu_indices, Index j) const
      {
        typedef typename LFSVIndices::ConstraintsIterator VConstraintsIterator;
        typedef typename LFSUIndices::ConstraintsIterator UConstraintsIterator;

        const bool constrained_v = lfsv_indices.constrained(i);
        const bool constrained_u = lfsu_indices.constrained(j);

        if(lfsv_indices.container_index(i) == lfsu_indices.container_index(j))
          globalpattern.add_link(lfsv_indices.container_index(i),lfsu_indices.container_index(j));

        if(!constrained_v)
          {
            if (!constrained_u || lfsu_indices.dirichlet_constraint(j))
              {
                globalpattern.add_link(lfsv_indices.container_index(i),lfsu_indices.container_index(j));
              }
            else
              {
                for (UConstraintsIterator gurit = lfsu_indices.constraints_begin(j); gurit != lfsu_indices.constraints_end(j); ++gurit)
                  globalpattern.add_link(lfsv_indices.container_index(i),gurit->containerIndex());
              }
          }
        else
          {
            if (lfsv_indices.dirichlet_constraint(i))
              {
                globalpattern.add_link(lfsv_indices.container_index(i),lfsu_indices.container_index(j));
              }
            else
              {
                for(VConstraintsIterator gvrit = lfsv_indices.constraints_begin(i); gvrit != lfsv_indices.constraints_end(i); ++gvrit)
                  {
                    if (!constrained_u || lfsu_indices.dirichlet_constraint(j))
                      {
                        globalpattern.add_link(gvrit->containerIndex(),lfsu_indices.container_index(j));
                      }
                    else
                      {
                        for (UConstraintsIterator gurit = lfsu_indices.constraints_begin(j); gurit != lfsu_indices.constraints_end(j); ++gurit)
                          globalpattern.add_link(gvrit->containerIndex(),gurit->containerIndex());
                      }
                  }
              }
          }
      }

      /** \brief insert dirichlet constraints for row and assemble
          T^T_U in constrained rows
      */
      template<typename GI, typename GC, typename CG>
      void set_trivial_row (GI i, const CG & cv_i, GC& globalcontainer) const
      {
        //std::cout << "clearing row " << i << std::endl;
        // set all entries in row i to zero
        B::clear_row(i,globalcontainer,1);

        // set diagonal element to 1
        // B::access(globalcontainer,i,i) = 1;
      }

      template<typename GC>
      void handle_dirichlet_constraints(GC& globalcontainer) const
      {
        B::flush(globalcontainer);
        typedef typename CV::const_iterator global_row_iterator;
        for (global_row_iterator cit=pconstraintsv->begin(); cit!=pconstraintsv->end(); ++cit)
          set_trivial_row(cit->first,cit->second,globalcontainer);
        B::finalize(globalcontainer);
      }

      /* constraints */
      const CU* pconstraintsu;
      const CV* pconstraintsv;
      static CU emptyconstraintsu;
      static CV emptyconstraintsv;
    };

    template<typename B, typename CU, typename CV>
    CU LocalAssemblerBase<B,CU,CV>::emptyconstraintsu;
    template<typename B, typename CU, typename CV>
    CV LocalAssemblerBase<B,CU,CV>::emptyconstraintsv;

  } // namespace PDELab
} // namespace Dune

#endif //DUNE_PDELAB_ASSEMBLERUTILITIES_HH
