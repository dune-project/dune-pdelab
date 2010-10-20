// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_MULTISTEP_CACHE_HH
#define DUNE_PDELAB_MULTISTEP_CACHE_HH

#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup MultiStepMethods Multi-Step Methods
    //! \ingroup PDELab
    //! \{

    //! Base exception for cache related errors
    class CacheError : public Exception {};

    //! Exception thrown when a requested item is not in the cache
    class NotInCache : public CacheError {};

    //! Exception thrown when a stored item is already in the cache
    class AlreadyInCache : public CacheError {};

    //! Policy class for the MultiStepCache
    /**
     * \tparam Step Type used for step number.  May be signed or unsigned
     *              integral.  If it is unsigned, no negative step number may
     *              occur, ever.
     * \tparam Time Type used for temporal values.  Should be floating-point.
     *
     * This is a virtual base class.  However, in addition to being a base
     * class it can also be used as-is, and will provide safe behaviour in
     * this case: cache everything, but don't assume any affine behaviour and
     * cache-items will be evicted once they are older than the number of
     * steps given in the contructor.
     *
     * \nosubgrouping
     */
    template<class Step = int, class Time = double>
    class MultiStepCachePolicy {
    protected:
      //! The number of steps in the scheme, set by preStep().
      /**
       * This means to compute value \f$u_n\f$, the scheme will need the next
       * \c stepsOfScheme old values, back to (and including)
       * \f$u_{n-\text{\tt stepsOfScheme}}\f$.
       */
      Step stepsOfScheme;
      //! The number of the step set via preStep().
      /**
       * After preStep(), this denotes the number of the step currently being
       * computed.  After postStep(), this denotes the number of the step that
       * has just been computed.
       */
      Step currentStep;

    public:
      //! any virtual base needs a virtual destructor
      virtual ~MultiStepCachePolicy() {}

      //! \name methods to determine what to cache
      //! \{

      //! whether to cache residual evaluations
      /** Returns \c true by default */
      virtual bool cacheResidualValue(std::size_t order, Step step) const
      { return true; }
      //! whether to cache jacobians (only relevant for affine operators)
      /** Returns \c true by default */
      virtual bool cacheJacobian(std::size_t order, Step step) const
      { return true; }
      //! whether to cache zero-residuals (only relevant for affine operators)
      /** Returns \c true by default */
      virtual bool cacheZeroResidual(std::size_t order, Step step) const
      { return true; }
      //! whether to cache composed jacobians (only relevant for affine operators)
      /** Returns \c true by default */
      virtual bool cacheComposedJacobian(Step step) const
      { return true; }

      //! \}

      //! \name methods to determine properties of the operators
      //! \{

      //! Whether the operator is affine
      /** Returns \c false by default */
      virtual bool isAffine(std::size_t order, Step step) const
      { return false; }
      //! Whether the composed operator is affine
      /** Returns \c false by default */
      virtual bool isComposedAffine(Step step) const
      { return false; }
      //! Whether the operator has purely linear alpha_*() methods
      /** Returns \c false by default */
      virtual bool hasPureLinearAlpha(std::size_t order, Step step) const
      { return false; }
      //! Whether a jacobian can be reused
      /** Returns \c false by default */
      virtual bool canReuseJacobian(std::size_t order,
                                    Step requested, Step available) const
      { return false; }
      //! Whether a zero-residual can be reused
      /** Returns \c false by default */
      virtual bool canReuseZeroResidual(std::size_t order,
                                        Step requested, Step available) const
      { return false; }
      //! Whether a composed jacobian can be reused
      /** Returns \c false by default */
      virtual bool canReuseComposedJacobian(Step requested,
                                            Step available) const
      { return false; }

      //! \}

      //! \name methods for pre-/post-processing
      //! \{

      //! method called before starting to compute a step
      /**
       * \param step           Number of the step that will be computed,
       *                       i.e. when \f$u_n\f$ is computed, then
       *                       \f$n=\f$\c step.
       * \param stepsOfScheme_ Number of steps in the scheme: the oldest
       *                       required value is \f$u_{\text{\tt
       *                       step}-\text{\tt stepsOfScheme}}\f$.
       * \param endTime        \f$t_n\f$, the time at step \c step.
       * \param dt             Time step size, i.e. \f$t_n-t_{n-1}\f$.
       *
       * The default implementation saves the parameters \c step and \c
       * stepsOfScheme_ in the protected member variables \c currentStep and
       * \c stepsOfScheme, respectively.
       */
      virtual void preStep(Step step, Step stepsOfScheme_,
                           Time endTime, Time dt)
      { currentStep = step; stepsOfScheme = stepsOfScheme_; }
      //! called after the new values for a time step have been computed
      /** does nothing by default. */
      virtual void postStep() { }

      //! \}

      //! \name methods to determine stale items
      //! \{

      //! determine whether a residual value can be removed
      virtual bool canEvictResidualValue(std::size_t order,
                                         Step step) const
      { return step + stepsOfScheme < currentStep; }
      //! determine whether a Jacobian can be removed
      virtual bool canEvictJacobian(std::size_t order, Step step) const
      { return step + stepsOfScheme < currentStep; }
      //! determine whether a zero-residual can be removed
      virtual bool canEvictZeroResidual(std::size_t order,
                                        Step step) const
      { return step + stepsOfScheme < currentStep; }
      //! determine whether a vector of unknowns can be removed
      virtual bool canEvictUnknowns(Step step) const
      { return step + stepsOfScheme < currentStep; }
      //! determine whether a composed Jacobian can be removed
      virtual bool canEvictComposedJacobian(Step step) const
      { return step + stepsOfScheme < currentStep; }
    };

    //! Cache for the CachedMultiStepGridOperatorSpace
    /**
     * \tparam VectorU Type of vectors for the unknowns.
     * \tparam VectorV Type of vectors for the residuals.
     * \tparam Matrix  Type of the Jacobians.
     * \tparam Step    Type of the step counters (should be integral).
     * \tparam Time    Type of the temporal values (should be floating-point).
     *
     * The MultiStepCache caches four kinds of data:
     * \li Residual values for a given time step and vector of unknowns:
     *     \f$r_j(t_n,u_n)\f$.  This is useful even for non-linear problems.
     * \li Jacobian matrices of the local operators \f$J(r_j|_{t_n})\f$.  This
     *     is mostly useful for affine operators when the Jacobians can be
     *     reused.
     * \li Zero-residuals of the local operators \f$r_j(t_n,0)\f$.  This is
     *     mostly useful for affine operators when the zero-residuals can be
     *     reused.
     * \li Jacobian matrices of the composed system
     *     \f$\sum_{j=0}^p\frac{\alpha{0j}}{(\Delta t)^j}J(r_j|_{t^n})\f$.
     *     These are the matrices which are solved by the solver in the end.
     *     This is mostly useful for affine operators when the Jacobians can be
     *     reused.
     *
     * In addition, there is one item which is store in the cache just like
     * the cached items, but is present mostly because it is a convenient way
     * for the user code to provide an manage that information:
     * \li The old values of the unknown vectors \f$u_n\f$.  This is not
     *     really a cache but a way for the user code to provide those values
     *     to the GridOperatorSpace.
     *
     * It is always valid for the cache implementation to silently refuse to
     * store a value (except for the vectors of unknowns).  The grid operator
     * space must never expect to be able to store a value and immediately be
     * able to extract it again.  The grid operator space must always be
     * prepared to recompute a value that cannot be extracted from the cache,
     * and should try to store that value in the cache afterwards.
     *
     * \note The cache keeps pointers to the values it stores.  The user code
     *       must make sure that any value stored in the cache is not later
     *       modified, any such modification results in undefined behaviour.
     */
    template<class VectorU, class VectorV, class Matrix,
             class Step = int, class Time = double>
    struct MultiStepCache {
      typedef MultiStepCachePolicy<Step, Time> Policy;

    private:
      typedef std::map<Step, shared_ptr<const Matrix> > MatrixMap;
      typedef typename MatrixMap::const_iterator MatrixIterator;
      typedef std::map<Step, shared_ptr<const VectorV> > ResidualMap;
      typedef typename ResidualMap::const_iterator ResidualIterator;
      typedef std::map<Step, shared_ptr<const VectorU> > UnknownMap;
      typedef typename UnknownMap::const_iterator UnknownIterator;

      // non-linear caching across time steps
      std::vector<ResidualMap> residualValues;

      // affine operators
      std::vector<MatrixMap> jacobians;
      std::vector<ResidualMap> zeroResiduals;

      // composed jacobians
      MatrixMap composedJacobians;

      // old values of the unknowns
      UnknownMap unknowns;

      // policy object
      shared_ptr<Policy> policy;

    public:

      //! \name construction and policy management
      //! \{

      MultiStepCache(const shared_ptr<Policy> &policy_ =
                           shared_ptr<Policy>(new Policy)) :
        policy(policy_)
      {
        if(!policy)
          DUNE_THROW(CacheError,
                     "MultiStepCache constructed with policy == NULL");
      }

      shared_ptr<Policy> getPolicy()
      { return policy; }
      shared_ptr<const Policy> getPolicy() const
      { return policy; }
      void setPolicy(const shared_ptr<Policy>& policy_)
      {
        if(!policy_)
          DUNE_THROW(CacheError, "MultiStepCache::setPolicy(): attempt to set "
                     "policy == NULL");
        policy = policy_;
      }

      //! \}

      //! \name methods for non

      //! \name non-linear caching across time steps
      //! \{

      //! get residual value from the cache
      /**
       * \param order Extract the residual value of the operator for the
       *              order'th temporal derivative.
       * \param step  Step for which to extract the residual value.
       *
       * \returns A shared pointer to the residual value vector.
       *
       * \throw NotInCache if the requested residual value is not in the
       *                   cache.
       */
      shared_ptr<const VectorV>
      getResidualValue(std::size_t order, Step step) const {
        if(order < residualValues.size()) {
          ResidualIterator it = residualValues[order].find(step);
          if(it != residualValues[order].end())
            return it->second;
        }
        DUNE_THROW(NotInCache, "MultiStepCache::getResidualValue(): The "
                   "requested residual value "
                   "r_" << order << "(t_" << step << ",u_" << step << ") is "
                   "not in the cache");
      }
      //! store a residual value in the cache
      /**
       * \param order         Store the residual value of the operator for the
       *                      order'th temporal derivative.
       * \param step          Step for which to store the residual value.
       * \param residualValue Pointer to the residual value to store.
       *
       * \throw AlreadyInCache if the cache already contains a residual value
       *                       for the given order and step.
       */
      void setResidualValue(std::size_t order, Step step,
                            const shared_ptr<const VectorV> &residualValue)
      {
        if(!policy->cacheResidualValue(order, step))
          return;
        if(order >= residualValues.size())
          residualValues.resize(order+1);
        if(!residualValues[order].insert(std::make_pair(step, residualValue))
           .second)
          DUNE_THROW(Exception, "Residual r_" << order << "(" << step << ", "
                     "0) is already in the cache!");
      }

      //! \}

      //! \name methods for the Jacobians of affine operators
      //! \{

      //! get a Jacobian from the cache
      /**
       * \param order Extract the Jacobian matrix of the operator for the
       *              order'th temporal derivative.
       * \param step  Step for which to extract the jacobian.
       *
       * \returns A shared pointer to the Jacobian matrix.
       *
       * \throw NotInCache if the requested Jacobian is not in the cache.
       *
       * \note This function will attempt to copy the Jacobian from other
       *       Jacobians of the same order.  As a consequence, this method
       *       only works on the mutable cache.
       */
      shared_ptr<const Matrix>
      getJacobian(std::size_t order, Step step) {
        if(order < jacobians.size()) {
          MatrixIterator it = jacobians[order].find(step);
          const MatrixIterator &end = jacobians[order].end();
          if(it != end)
            return it->second;

          // try to copy from another step
          for(it = jacobians[order].begin(); it != end; ++it)
            if(policy->canReuseJacobian(order, step, it->first))
              // assign and return value
              return jacobians[order][step] = it->second;
        }
        DUNE_THROW(NotInCache, "MultiStepCache::getJacobian(): The requested "
                   "Jacobian J(r_" << order << "|_t_" << step << ") is not in "
                   "the cache");
      }
      //! store a Jacobian in the cache
      /**
       * \param order    Store the Jacobian of the operator for the order'th
       *                 temporal derivative.
       * \param step     Step for which to store the jacobian.
       * \param jacobian Pointer to the Jacobian to store.
       *
       * \throw AlreadyInCache if the cache already contains a Jacobian for
       *                       the given order and step.
       */
      void setJacobian(std::size_t order, Step step,
                       const shared_ptr<const Matrix> &jacobian)
      {
        if(!policy->cacheJacobian(order, step))
          return;
        if(order >= jacobians.size())
          jacobians.resize(order+1);
        if(!jacobians[order].insert(std::make_pair(step, jacobian)).second)
          DUNE_THROW(AlreadyInCache, "MultiStepCache::setJacobian(): Jacobian "
                     "J(r_" << order << "|_t_" << step << ") is already in "
                     "the cache!");
      }

      //! \}

      //! \name methods for the zero-residual of affine operators
      //! \{

      //! get zero-residual from the cache
      /**
       * \param order Extract the zero-residual of the operator for the
       *              order'th temporal derivative.
       * \param step  Step for which to extract the zero-residual.
       *
       * \returns A shared pointer to the zero-residual vector.
       *
       * \throw NotInCache if the requested zero-residual is not in the cache.
       *
       * \note This function will attempt to copy the zero-residual from other
       *       zero-residuals of the same order.  As a consequence, this
       *       method only works on the mutable cache.
       */
      shared_ptr<const VectorV>
      getZeroResidual(std::size_t order, Step step) const {
        if(order < zeroResiduals.size()) {
          ResidualIterator it = zeroResiduals[order].find(step);
          const ResidualIterator &end = zeroResiduals[order].end();
          if(it != end)
            return it->second;

          // try to copy from another step
          for(it = zeroResiduals[order].begin(); it != end; ++it)
            if(policy->canReuseZeroResidual(order, step, it->first))
              // assign and return value
              return zeroResiduals[order][step] = it->second;
        }
        DUNE_THROW(NotInCache, "MultiStepCache::getZeroResidual(): The "
                   "requested zero-residual "
                   "r_" << order << "(t_" << step << ",0) is not in the "
                   "cache");
      }
      //! store a zero-residual in the cache
      /**
       * \param order        Store the zero-residual of the operator for the
       *                     order'th temporal derivative.
       * \param step         Step for which to store the zero-residual.
       * \param zeroResidual Pointer to the zero-residual to store.
       *
       * \throw AlreadyInCache if the cache already contains a zero-residual
       *                       for the given order and step.
       */
      void setZeroResidual(std::size_t order, Step step,
                           const shared_ptr<const VectorV> &zeroResidual)
      {
        if(!policy->cacheZeroResidual(order, step))
          return;
        if(order >= zeroResiduals.size())
          zeroResiduals.resize(order+1);
        if(!zeroResiduals[order].insert(std::make_pair(step, zeroResidual))
           .second)
          DUNE_THROW(Exception, "Residual r_" << order << "(" << step << ", "
                     "0) is already in the cache!");
      }

      //! \}

      //! \name methods for the composed Jacobian of affine operators
      //! \{

      //! get a composed Jacobian from the cache
      /**
       * \param step Step for which to extract the jacobian.
       *
       * \returns A shared pointer to the Jacobian matrix.
       *
       * \throw NotInCache if the requested composed Jacobian is not in the
       *                   cache.
       *
       * \note This function will attempt to copy the composed Jacobian from
       *       composed Jacobians of other times steps.  As a consequence,
       *       this method only works on the mutable cache.
       */
      shared_ptr<const Matrix>
      getComposedJacobian(Step step) {
        MatrixIterator it = composedJacobians.find(step);
        const MatrixIterator &end = composedJacobians.end();
        if(it != end)
          return it->second;

        // try to copy from another step
        for(it = composedJacobians.begin(); it != end; ++it)
          if(policy->canReuseComposedJacobian(step, it->first))
            // assign and return value
            return composedJacobians[step] = it->second;

        DUNE_THROW(NotInCache, "MultiStepCache::getComposedJacobian(): The "
                   "requested composed Jacobian for step " << step << " is "
                   "not in the cache");
      }
      //! store a composed Jacobian in the cache
      /**
       * \param step     Step for which to store the composed Jacobian.
       * \param jacobian Pointer to the composed Jacobian to store.
       *
       * \throw AlreadyInCache if the cache already contains a composed
       *                       Jacobian for the given step.
       */
      void setComposedJacobian(Step step,
                               const shared_ptr<const Matrix> &jacobian)
      {
        if(!policy->cacheComposedJacobian(step))
          return;
        if(!composedJacobians.insert(std::make_pair(step, jacobian)).second)
          DUNE_THROW(AlreadyInCache, "MultiStepCache::setComposedJacobian(): "
                     "Composed Jacobian for time step " << step << " is "
                     "already in the cache!");
      }

      //! \}

      //! \name methods to access old values of the unknowns
      //! \{

      //! get vector of unknowns from the cache
      /**
       * \param step Step for which to extract the vector of unknowns.
       *
       * \returns A shared pointer to the vector of unknowns.
       *
       * \throw NotInCache if the cache does not contain a vector of unknowns
       *                   for the given step.
       */
      shared_ptr<const VectorU>
      getUnknowns(Step step) const {
        UnknownIterator it = unknowns.find(step);
        if(it != unknowns.end())
          return it->second;
        DUNE_THROW(Exception, "Unknowns for time step " << step << " missing "
                   "in the cache!");
      }
      //! store a vector of unknowns  in the cache
      /**
       * \param step      Step for which to store the vector of unknowns.
       * \param unknowns_ Pointer to the unknowns to store.
       *
       * \throw AlreadyInCache if the cache already contains a vector of
       *                       unknowns for the given order and step.
       */
      void setUnknowns(Step step, const shared_ptr<const VectorU> &unknowns_) {
        if(unknowns[step])
          DUNE_THROW(Exception, "Unknowns u_" << step << " are already in the "
                     "cache!");
        unknowns[step] = unknowns_;
      }

      //! \}

      //! \name methods to flush the cache
      //! \{

      //! Flush all cached values
      /**
       * This is useful for instance after adaption.  It is equivalent to
       * recreating the cache.
       */
      void flushAll() {
        residualValues.clear();
        jacobians.clear();
        zeroResiduals.clear();
        composedJacobians.clear();
        unknowns.clear();
      }

      //! \}

      //! \name methods for pre-/post-processing

      //! Do some housekeeping before computing a time-step
      /**
       * The method first calls preStep on the policy object, passing it the
       * parameters:
       *
       * \param step          Number of the step that will be computed,
       *                      i.e. when \f$u_n\f$ is computed, then
       *                      \f$n=\f$\c step.
       * \param stepsOfScheme Number of steps in the scheme: the oldest
       *                      required value is \f$u_{\text{\tt
       *                      step}-\text{\tt stepsOfScheme}}\f$.
       * \param endTime       \f$t_n\f$, the time at step \c step.
       * \param dt            Time step size, i.e. \f$t_n-t_{n-1}\f$.
       *
       * It the iterates over each cache item and remove it, if the policy
       * allows that.
       */
      void preStep(Step step, Step stepsOfScheme, Time endTime, Time dt) {
        policy->preStep(step, stepsOfScheme, endTime, dt);

        // residual values
        for(std::size_t order = 0; order < residualValues.size(); ++order) {
          const ResidualIterator end = residualValues[order].end();
          ResidualIterator it = residualValues[order].begin();
          while(it != end)
            if(policy->canEvictResidualValue(order, it->first))
              residualValues[order].erase(it++);
            else
              ++it;
        }

        // Jacobians
        for(std::size_t order = 0; order < jacobians.size(); ++order) {
          const MatrixIterator end = jacobians[order].end();
          MatrixIterator it = jacobians[order].begin();
          while(it != end)
            if(policy->canEvictJacobian(order, it->first))
              jacobians[order].erase(it++);
            else
              ++it;
        }

        // zero-residual
        for(std::size_t order = 0; order < zeroResiduals.size(); ++order) {
          const ResidualIterator end = zeroResiduals[order].end();
          ResidualIterator it = zeroResiduals[order].begin();
          while(it != end)
            if(policy->canEvictZeroResidual(order, it->first))
              zeroResiduals[order].erase(it++);
            else
              ++it;
        }

        // composed Jacobians
        {
          const MatrixIterator end = composedJacobians.end();
          MatrixIterator it = composedJacobians.begin();
          while(it != end)
            if(policy->canEvictComposedJacobian(it->first))
              composedJacobians.erase(it++);
            else
              ++it;
        }

        // unknowns
        {
          const UnknownIterator end = unknowns.end();
          UnknownIterator it = unknowns.begin();
          while(it != end)
            if(policy->canEvictUnknowns(it->first))
              unknowns.erase(it++);
            else
              ++it;
        }
      }

      //! Do some housekeeping after computing a time-step
      /**
       * This just calls postStep() on the policy object.
       */
      void postStep() { policy->postStep(); }

      //! \}

    };

    //! \} group MultiStepMethods
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTISTEP_CACHE_HH
