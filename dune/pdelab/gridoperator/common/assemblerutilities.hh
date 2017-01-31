// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDOPERATOR_COMMON_ASSEMBLERUTILITIES_HH
#define DUNE_PDELAB_GRIDOPERATOR_COMMON_ASSEMBLERUTILITIES_HH

#include <algorithm>
#include <tuple>

#include <dune/pdelab/constraints/common/constraintstransformation.hh>
#include <dune/pdelab/gridoperator/common/localmatrix.hh>

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
      typedef typename MatrixBackend::template Pattern<
        Jacobian,
        TestGridFunctionSpace,
        TrialGridFunctionSpace
        > MatrixPattern;

      //! The helper class to exchange data on the processor boundary
      typedef typename GO::BorderDOFExchanger BorderDOFExchanger;

    };

    // ********************************************************************************
    // default local pattern implementation
    // ********************************************************************************

    /**
       \brief Entry in sparsity pattern

       The sparsity pattern of a linear operator is described by by connecting
       degrees of freedom in one element with degrees of freedom in the
       same element (intra) or an intersecting element (inter).

       This numbering is with respect to the depth-first canonical order of the
       degrees of freedom of an entity.

       \nosubgrouping
    */
    class SparsityLink : public std::tuple<int,int>
    {
    public:
      //! \brief Standard constructor for uninitialized local index
      SparsityLink ()
      {}

      //! \brief Initialize all components
      SparsityLink (int i, int j)
        : std::tuple<int,int>(i,j)
      {}

      //! \brief Return first component
      inline int i () const
      {
        return std::get<0>(*this);
      }

      //! \brief Return second component
      inline int j () const
      {
        return std::get<1>(*this);
      }

      //! \brief Set both components
      void set (int i, int j)
      {
        std::get<0>(*this) = i;
        std::get<1>(*this) = j;
      }
    };

    /**
       \brief Layout description for a sparse linear operator
       \see SparsityLink

       \nosubgrouping
    */
    class LocalSparsityPattern
      : public std::vector<SparsityLink>
    {

      // make push_back() inaccessible
      using std::vector<SparsityLink>::push_back;

    public:

      //! Adds a link between DOF i of lfsv and DOF j of lfsu.
      /**
       * This methods adds a link between the DOF i of the test local test space lfsv
       * and the DOF j of the local ansatz space lfsu.
       *
       * \param lfsv  The local test space.
       * \param i     Index of the DOF in the test space lfsv.
       * \param lfsu  The local ansatz space.
       * \param j     Index of the DOF in the ansatz space lfsu.
       */
      template<typename LFSV, typename LFSU>
      void addLink(const LFSV& lfsv, std::size_t i, const LFSU& lfsu, std::size_t j)
      {
        std::vector<SparsityLink>::push_back(
          SparsityLink(
            lfsv.localIndex(i),
            lfsu.localIndex(j)
          )
        );
      }

    };



    // ********************************************************************************
    // Assembler base class
    // ********************************************************************************

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
    class LocalAssemblerBase
    {
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
      typename std::enable_if<
        AlwaysTrue<X>::value && !std::is_same<
          CV,
          EmptyTransformation
          >::value
        >::type
      forwardtransform(X & x, const bool postrestrict = false) const
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
            {
              // typename X::block_type block(x[contributor]);
              // block *= it->second;
              // x[it->first] += block;
              x[it->first] += it->second * x[contributor];
            }
        }

        if(postrestrict)
          for (global_col_iterator cit=pconstraintsv->begin(); cit!=pconstraintsv->end(); ++cit)
            x[cit->first]=0.;
      }


      // Disable forwardtransform for EmptyTransformation
      template<typename X>
      typename std::enable_if<
        AlwaysTrue<X>::value && std::is_same<
          CV,
          EmptyTransformation
          >::value
        >::type
      forwardtransform(X & x, const bool postrestrict = false) const
      {}


      /** \brief Transforms a vector \f$ \boldsymbol{x} \f$ from \f$
          V'\f$ to \f$ V\f$. If prerestrict == true then
          \f$\boldsymbol{S}^T_{\boldsymbol{\tilde U}}\f$ is applied
          instead of the full transformation.  */
      template<typename X>
      typename std::enable_if<
        AlwaysTrue<X>::value && !std::is_same<
          CV,
          EmptyTransformation
          >::value
        >::type
      backtransform(X & x, const bool prerestrict = false) const
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
            {
              // typename X::block_type block(x[it->first]);
              // block *= it->second;
              // x[contributor] += block;
              x[contributor] += it->second * x[it->first]; // PB: 27 Sep 12 this was the old version
            }
        }
      }

      // disable backtransform for empty transformation
      template<typename X>
      typename std::enable_if<
        AlwaysTrue<X>::value && std::is_same<
          CV,
          EmptyTransformation
          >::value
        >::type
      backtransform(X & x, const bool prerestrict = false) const
      {}


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

      //! Scatter local jacobian to global container.
      template<typename M, typename GCView>
      typename std::enable_if<
        AlwaysTrue<M>::value && !std::is_same<
          CV,
          EmptyTransformation
          >::value
      >::type
      scatter_jacobian(M& local_container, GCView& global_container_view, bool symmetric_mode) const
      {
        typedef typename GCView::RowIndexCache LFSVIndexCache;
        typedef typename GCView::ColIndexCache LFSUIndexCache;

        const LFSVIndexCache& lfsv_indices = global_container_view.rowIndexCache();
        const LFSUIndexCache& lfsu_indices = global_container_view.colIndexCache();

        if (lfsv_indices.constraintsCachingEnabled() && lfsu_indices.constraintsCachingEnabled())
          if (symmetric_mode)
            etadd_symmetric(local_container,global_container_view);
          else
            etadd(local_container,global_container_view);
        else
          {

            typedef typename LFSVIndexCache::LocalFunctionSpace LFSV;
            const LFSV& lfsv = lfsv_indices.localFunctionSpace();

            typedef typename LFSUIndexCache::LocalFunctionSpace LFSU;
            const LFSU& lfsu = lfsu_indices.localFunctionSpace();

            // optionally clear out columns that belong to Dirichlet-constrained DOFs to keep matrix symmetric
            if (symmetric_mode)
              {
                typedef typename LFSUIndexCache::ContainerIndex CI;

                for (size_t j = 0; j < lfsu_indices.size(); ++j)
                  {
                    const CI& container_index = lfsu_indices.containerIndex(j);
                    const typename CU::const_iterator cit = pconstraintsu->find(container_index);
                    if (cit != pconstraintsu->end())
                      {
                        // make sure we only have Dirichlet constraints
                        assert(cit->second.empty());
                        // clear out the current column
                        for (size_t i = 0; i < lfsv_indices.size(); ++i)
                          {
                            // we do not need to update the residual, since the solution
                            // (i.e. the correction) for the Dirichlet DOF is 0 by definition
                            local_container(lfsv,i,lfsu,j) = 0.0;
                          }
                      }
                  }
              }

            // write entries without considering constraints.
            // Dirichlet-constrained rows will be fixed in a postprocessing step.
            for (auto it = local_container.begin(); it != local_container.end(); ++it)
              {
                // skip 0 entries because they might not be present in the pattern
                if (*it == 0.0)
                  continue;
                global_container_view.add(it.row(),it.col(),*it);
              }
          }
      }

      // specialization for empty constraints container
      template<typename M, typename GCView>
      typename std::enable_if<
        AlwaysTrue<M>::value && std::is_same<
          CV,
          EmptyTransformation
          >::value
      >::type
      scatter_jacobian(M& local_container, GCView& global_container_view, bool symmetric_mode) const
      {
        // write entries without considering constraints.
        // Dirichlet-constrained rows will be fixed in a postprocessing step.
        for (auto it = local_container.begin(); it != local_container.end(); ++it)
          {
            // skip 0 entries because they might not be present in the pattern
            if (*it == 0.0)
              continue;
            global_container_view.add(it.row(),it.col(),*it);
          }
      }

      /** \brief Add local matrix to global matrix,
          and apply Dirichlet constraints in a symmetric
          fashion. Apart from that, identical to etadd(). */
      template<typename M, typename GCView>
      void etadd_symmetric (M& localcontainer, GCView& globalcontainer_view) const
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
            if (lfsu_indices.isConstrained(j) && lfsu_indices.isDirichletConstraint(j))
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


      template<typename M, typename GCView>
      void etadd (const M& localcontainer, GCView& globalcontainer_view) const
      {

        typedef typename M::value_type T;

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

              const bool constrained_v = lfsv_indices.isConstrained(i);
              const bool constrained_u = lfsu_indices.isConstrained(j);

              typedef typename LFSVIndexCache::ConstraintsIterator VConstraintsIterator;
              typedef typename LFSUIndexCache::ConstraintsIterator UConstraintsIterator;

              if (constrained_v)
                {
                  if (lfsv_indices.isDirichletConstraint(i))
                    continue;

                  for (VConstraintsIterator vcit = lfsv_indices.constraintsBegin(i); vcit != lfsv_indices.constraintsEnd(i); ++vcit)
                    if (constrained_u)
                      {
                        if (lfsu_indices.isDirichletConstraint(j))
                          {
                            T value = localcontainer(lfsv,i,lfsu,j) * vcit->weight();
                            if (value != 0.0)
                              globalcontainer_view.add(vcit->containerIndex(),j,value);
                          }
                        else
                          {
                            for (UConstraintsIterator ucit = lfsu_indices.constraintsBegin(j); ucit != lfsu_indices.constraintsEnd(j); ++ucit)
                              {
                                T value = localcontainer(lfsv,i,lfsu,j) * vcit->weight() * ucit->weight();
                                if (value != 0.0)
                                  globalcontainer_view.add(vcit->containerIndex(),ucit->containerIndex(),value);
                              }
                          }
                      }
                    else
                      {
                        T value = localcontainer(lfsv,i,lfsu,j) * vcit->weight();
                        if (value != 0.0)
                          globalcontainer_view.add(vcit->containerIndex(),j,value);
                      }
                }
              else
                {
                  if (constrained_u)
                    {
                      if (lfsu_indices.isDirichletConstraint(j))
                        {
                          T value = localcontainer(lfsv,i,lfsu,j);
                          if (value != 0.0)
                            globalcontainer_view.add(i,j,value);
                        }
                      else
                        {
                          for (UConstraintsIterator ucit = lfsu_indices.constraintsBegin(j); ucit != lfsu_indices.constraintsEnd(j); ++ucit)
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
      typename std::enable_if<
        std::is_same<RI,CI>::value
        >::type
      add_diagonal_entry(Pattern& pattern, const RI& ri, const CI& ci) const
      {
        if (ri == ci)
          pattern.add_link(ri,ci);
      }

      template<typename Pattern, typename RI, typename CI>
      typename std::enable_if<
        !std::is_same<RI,CI>::value
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

        const bool constrained_v = lfsv_indices.isConstrained(i);
        const bool constrained_u = lfsu_indices.isConstrained(j);

        add_diagonal_entry(globalpattern,
                           lfsv_indices.containerIndex(i),
                           lfsu_indices.containerIndex(j)
                           );

        if(!constrained_v)
          {
            if (!constrained_u || lfsu_indices.isDirichletConstraint(j))
              {
                globalpattern.add_link(lfsv_indices.containerIndex(i),lfsu_indices.containerIndex(j));
              }
            else
              {
                for (UConstraintsIterator gurit = lfsu_indices.constraintsBegin(j); gurit != lfsu_indices.constraintsEnd(j); ++gurit)
                  globalpattern.add_link(lfsv_indices.containerIndex(i),gurit->containerIndex());
              }
          }
        else
          {
            if (lfsv_indices.isDirichletConstraint(i))
              {
                globalpattern.add_link(lfsv_indices.containerIndex(i),lfsu_indices.containerIndex(j));
              }
            else
              {
                for(VConstraintsIterator gvrit = lfsv_indices.constraintsBegin(i); gvrit != lfsv_indices.constraintsEnd(i); ++gvrit)
                  {
                    if (!constrained_u || lfsu_indices.isDirichletConstraint(j))
                      {
                        globalpattern.add_link(gvrit->containerIndex(),lfsu_indices.containerIndex(j));
                      }
                    else
                      {
                        for (UConstraintsIterator gurit = lfsu_indices.constraintsBegin(j); gurit != lfsu_indices.constraintsEnd(j); ++gurit)
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

#endif // DUNE_PDELAB_GRIDOPERATOR_COMMON_ASSEMBLERUTILITIES_HH
