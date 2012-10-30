// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_ASSEMBLERUTILITIES_HH
#define DUNE_PDELAB_ASSEMBLERUTILITIES_HH

#include <dune/pdelab/common/unordered_map.hh>
#include <dune/pdelab/common/unordered_set.hh>

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

      //! Data structure for storing border-border matrix pattern entries in a communication-optimized form
      typedef unordered_map<
        typename GO::Traits::TestGridFunctionSpace::Ordering::Traits::DOFIndex,
        unordered_set<
          typename GO::Traits::TrialGridFunctionSpace::Ordering::Traits::GlobalDOFIndex
          >
        > BorderPattern;

      typedef typename GO::BorderDOFExchanger BorderDOFExchanger;

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
      template<typename T, typename GCView>
      void etadd_symmetric (const LocalMatrix<T>& localcontainer, GCView& globalcontainer_view) const
      {

        typedef typename GCView::RowIndexCache LFSVIndexCache;
        typedef typename GCView::ColIndexCache LFSUIndexCache;

        const LFSVIndexCache& lfsv_indices = globalcontainer_view.rowIndexCache();
        const LFSUIndexCache& lfsu_indices = globalcontainer_view.colIndexCache();

        typedef typename LFSVIndexCache::LocalFunctionSpace LFSV;
        const LFSV& lfsv = lfsv_indices.localFunctionSpace();

        typedef typename LFSUIndexCache::LocalFunctionSpace LFSU;
        const LFSU& lfsu = lfsu_indices.localFunctionSpace();

        for (size_t j = 0; j < lfsu_indices.size(); ++j)
          {
            if (lfsu_indices.constrained(j) && lfsu_indices.dirichlet_constraint(j))
              {
                // clear out the current column
                for (size_t i = 0; i < lfsv_indices.size(); ++i)
                  {
                    // we do not need to update the residual, since the solution
                    // (i.e. the correction) for the Dirichlet DOF is 0 by definition
                    localcontainer(lfsv,i,lfsu,j) = 0.0;
                  }
              }
          }

        // hand off to standard etadd() method
        etadd(localcontainer,globalcontainer_view);
      }


      template<typename T, typename GCView>
      void etadd (const LocalMatrix<T>& localcontainer, GCView& globalcontainer_view) const
      {

        typedef typename GCView::RowIndexCache LFSVIndexCache;
        typedef typename GCView::ColIndexCache LFSUIndexCache;

        const LFSVIndexCache& lfsv_indices = globalcontainer_view.rowIndexCache();
        const LFSUIndexCache& lfsu_indices = globalcontainer_view.colIndexCache();

        typedef typename LFSVIndexCache::LocalFunctionSpace LFSV;
        const LFSV& lfsv = lfsv_indices.localFunctionSpace();

        typedef typename LFSUIndexCache::LocalFunctionSpace LFSU;
        const LFSU& lfsu = lfsu_indices.localFunctionSpace();

        for (size_t i = 0; i<lfsv_indices.size(); ++i)
          for (size_t j = 0; j<lfsu_indices.size(); ++j)
            {

              if (localcontainer(lfsv,i,lfsu,j) == 0.0)
                continue;

              const bool constrained_v = lfsv_indices.constrained(i);
              const bool constrained_u = lfsu_indices.constrained(j);

              typedef typename LFSVIndexCache::ConstraintsIterator VConstraintsIterator;
              typedef typename LFSUIndexCache::ConstraintsIterator UConstraintsIterator;

              if (constrained_v)
                {
                  if (lfsv_indices.dirichlet_constraint(i))
                    continue;

                  for (VConstraintsIterator vcit = lfsv_indices.constraints_begin(i); vcit != lfsv_indices.constraints_end(i); ++vcit)
                    if (constrained_u)
                      {
                        if (lfsu_indices.dirichlet_constraint(j))
                          {
                            T value = localcontainer(lfsv,i,lfsu,j) * vcit->weight();
                            if (value != 0.0)
                              globalcontainer_view.add(vcit->containerIndex(),j,value);
                          }
                        else
                          {
                            for (UConstraintsIterator ucit = lfsu_indices.constraints_begin(j); ucit != lfsu_indices.constraints_end(j); ++ucit)
                              {
                                T value = localcontainer(lfsv,i,lfsu,j) * vcit->weight() * ucit->weight();
                                if (value != 0.0)
                                  globalcontainer_view.add(vcit->containerIndex(),ucit->containerIndex(),value);
                              }
                          }
                      }
                }
              else
                {
                  if (constrained_u)
                    {
                      if (lfsu_indices.dirichlet_constraint(j))
                        {
                          T value = localcontainer(lfsv,i,lfsu,j);
                          if (value != 0.0)
                            globalcontainer_view.add(i,j,value);
                        }
                      else
                        {
                          for (UConstraintsIterator ucit = lfsu_indices.constraints_begin(j); ucit != lfsu_indices.constraints_end(j); ++ucit)
                            {
                              T value = localcontainer(lfsv,i,lfsu,j) * ucit->weight();
                              if (value != 0.0)
                                globalcontainer_view.add(i,ucit->containerIndex(),value);
                            }
                        }
                    }
                  else
                    globalcontainer_view.add(i,j,localcontainer(lfsv,i,lfsu,j));
                }
            }
      }


      template<typename Pattern, typename RI, typename CI>
      typename enable_if<
        is_same<RI,CI>::value
        >::type
      add_diagonal_entry(Pattern& pattern, const RI& ri, const CI& ci) const
      {
        if (ri == ci)
          pattern.add_link(ri,ci);
      }

      template<typename Pattern, typename RI, typename CI>
      typename enable_if<
        !is_same<RI,CI>::value
        >::type
      add_diagonal_entry(Pattern& pattern, const RI& ri, const CI& ci) const
      {}

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

        add_diagonal_entry(globalpattern,
                           lfsv_indices.container_index(i),
                           lfsu_indices.container_index(j)
                           );

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
      template<typename GFSV, typename GC, typename C>
      void set_trivial_rows(const GFSV& gfsv, GC& globalcontainer, const C& c) const
      {
        typedef typename C::const_iterator global_row_iterator;
        for (global_row_iterator cit = c.begin(); cit != c.end(); ++cit)
          if (cit->second.size() == 0)
            globalcontainer.clear_row(cit->first,1);
      }

      template<typename GFSV, typename GC>
      void set_trivial_rows(const GFSV& gfsv, GC& globalcontainer, const EmptyTransformation& c) const
      {
      }

      template<typename GFSV, typename GC>
      void handle_dirichlet_constraints(const GFSV& gfsv, GC& globalcontainer) const
      {
        globalcontainer.flush();
        set_trivial_rows(gfsv,globalcontainer,*pconstraintsv);
        globalcontainer.finalize();
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
