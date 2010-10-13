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

    //! Cache for the CachingMultiStepGridOperatorSpace
    /**
     * \tparam VectorU Type of vectors for the unknowns.
     * \tparam VectorV Type of vectors for the residuals.
     * \tparam Matrix  Type of the Jacobians.
     * \tparam Step    Type of the step counters (should be integral).
     *
     * The MultiStepCache caches for kinds of data:
     * \li Residual values for a given time step and vector of unknowns:
     *     \f$r_j(t_n,u_n)\f$.  This is useful even for non-linear problems.
     * \li Old values of the unknown vectors \f$u_n\f$.  This is not really a
     *     cache but a way for the user code to provide those values to the
     *     GridOperatorSpace.
     * \li Jacobian matrices of the local operators \f$J(r_j|_{t_n})\f$.  This
     *     is mostly useful fo affine operators when the Jacobians can be
     *     reused.
     * \li Zero-residuals of the local operators \f$r_j(t_n,0)\f$.  This is
     *     mostly useful fo affine operators when the zero-residuals can be
     *     reused.
     *
     * It is always valid for the cache implementation to silently refuse to
     * store a value (except for the vectors of unknowns).  The grid operator
     * space must never expect to be able to store a value and immediately be
     * able to extract it again.  The grid operator space must always be
     * prepared to recompute a value that cannot be extracted from the cache,
     * and should store that value in the cache afterwards.
     *
     * \note The cache keeps pointers to the values it stores.  The user code
     *       must make sure that any value stored in the cache is not later
     *       modified, any such modification results in undefined behaviour.
     */
    template<class VectorU, class VectorV, class Matrix, class Step = int>
    class MultiStepCache {
      typedef std::map<Step, shared_ptr<const Matrix> > MatrixMap;
      typedef typename MatrixMap::const_iterator MatrixIterator;
      typedef std::map<Step, shared_ptr<const VectorV> > ResidualMap;
      typedef typename ResidualMap::const_iterator ResidualIterator;
      typedef std::map<Step, shared_ptr<const VectorU> > UnknownMap;
      typedef typename UnknownMap::const_iterator UnknownIterator;

      // affine operators
      std::vector<MatrixMap> jacobians;
      std::vector<ResidualMap> zeroResiduals;

      // non-linear caching across time steps
      std::vector<ResidualMap> residualValues;

      // old values of the unknowns
      UnknownMap unknowns;

    public:
      ////////////////////////////////////////////////////////////////////////
      //
      // affine operators
      //

      //! get a Jacobian from the cache
      /**
       * \param order Extract the Jacobian matrix of the operator for the
       *              order'th temporal derivative.
       * \param step  Step for which to extract the jacobian.
       *
       * \returns A shared pointer to the Jacobian matrix.
       *
       * \throw NotInCache if the requested Jacobian is not in the cache.
       */
      shared_ptr<const Matrix>
      getJacobian(std::size_t order, Step step) const {
        if(order < jacobians.size()) {
          MatrixIterator it = jacobians[order].find(step);
          if(it != jacobians[order].end())
            return it->second;
        }
        DUNE_THROW(NotInCache, "MultiStepCache::getJacobian(): The requested "
                   "Jacobian J(r_" << order << "|_t_" << step << ") is not in "
                   "the cache");
      }
      //! get a the latest Jacobian from the cache
      /**
       * \param order Extract the Jacobian of the operator for the order'th
       *              temporal derivative.
       * \param step  Step of the extracted jacobian.
       *
       * \returns A shared pointer to the latest Jacobian matrix for the given
       *          order.
       *
       * \throw NotInCache if there is no Jacobian for the requested order in
       *                   the cache.
       */
      shared_ptr<const Matrix>
      getLatestJacobian(std::size_t order, Step &step) const {
        if(order < jacobians.size() && !jacobians[order].empty()) {
          step = jacobians[order].rbegin()->first;
          return jacobians[order].rbegin()->second;
        }
        DUNE_THROW(NotInCache, "MultiStepCache::getLatestJacobian(): No "
                   "Jacobian J(r_" << order << ") in the cache");
      }
      //! store a Jacobian in the cache
      /**
       * \param order Store the Jacobian of the operator for the order'th
       *              temporal derivative.
       * \param step  Step for which to store the jacobian.
       *
       * \throw AlreadyInCache if the cache already contains a Jacobian for
       *                       the given order and step.
       */
      void setJacobian(Step step, std::size_t order,
                       const shared_ptr<const Matrix> &jacobian)
      {
        if(order >= jacobians.size())
          jacobians.resize(order+1);
        if(!jacobians.insert(std::make_pair(step, jacobian)).second)
          DUNE_THROW(AlreadyInCache, "MultiStepCache::setJacobian(): Jacobian "
                     "J(r_" << order << "|_t_" << step << ") is already in "
                     "the cache!");
      }
      //! get zero-residual from the cache
      /**
       * \param order Extract the zero-residual of the operator for the
       *              order'th temporal derivative.
       * \param step  Step for which to extract the zero-residual.
       *
       * \returns A shared pointer to the zero-residual vector.
       *
       * \throw NotInCache if the requested zero-residual is not in the cache.
       */
      shared_ptr<const VectorV>
      getZeroResidual(std::size_t order, Step step) const {
        if(order < zeroResiduals.size()) {
          ResidualIterator it = zeroResiduals[order].find(step);
          if(it != zeroResiduals[order].end())
            return it->second;
        }
        DUNE_THROW(NotInCache, "MultiStepCache::getZeroResidual(): The "
                   "requested zero-residual "
                   "r_" << order << "(t_" << step << ",0) is not in the "
                   "cache");
      }
      //! get the latest zero-residual vector from the cache
      /**
       * \param order Extract the zero-residual of the operator for the
       *              order'th temporal derivative.
       * \param step  Step of the extracted zero-residual.
       *
       * \returns A shared pointer to the zero-residual vector.
       *
       * \throw NotInCache if there is no zero-residual for the requested
       *                   order in the cache.
       */
      shared_ptr<const VectorV>
      getLatestZeroResidual(std::size_t order, Step &step) const {
        if(order < zeroResiduals.size() && !zeroResiduals[order].empty()) {
          step = zeroResiduals[order].rbegin()->first;
          return zeroResiduals[order].rbegin()->second;
        }
        DUNE_THROW(NotInCache, "MultiStepCache::getZeroResidual(): No "
                   "zero-residual r_" << order << "(t,0) in the cache");
      }
      //! store a zero-residual in the cache
      /**
       * \param order Store the zero-residual of the operator for the order'th
       *              temporal derivative.
       * \param step  Step for which to store the zero-residual.
       *
       * \throw AlreadyInCache if the cache already contains a zero-residual
       *                       for the given order and step.
       */
      void setZeroResidual(Step step, std::size_t order,
                           const shared_ptr<const VectorV> &zeroResidual)
      {
        if(order >= zeroResiduals.size())
          zeroResiduals.resize(order+1);
        if(!zeroResiduals.insert(std::make_pair(step, zeroResidual)).second)
          DUNE_THROW(Exception, "Residual r_" << order << "(" << step << ", "
                     "0) is already in the cache!");
      }

      ////////////////////////////////////////////////////////////////////////
      //
      // non-linear caching across time steps
      //

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
       * \param order Store the residual value of the operator for the
       *              order'th temporal derivative.
       * \param step  Step for which to store the residual value.
       *
       * \throw AlreadyInCache if the cache already contains a residual value
       *                       for the given order and step.
       */
      void setResidualValue(Step step, std::size_t order,
                            const shared_ptr<const VectorV> &residualValue)
      {
        if(order >= residualValues.size())
          residualValues.resize(order+1);
        if(!residualValues.insert(std::make_pair(step, residualValue)).second)
          DUNE_THROW(Exception, "Residual r_" << order << "(" << step << ", "
                     "0) is already in the cache!");
      }

      ////////////////////////////////////////////////////////////////////////
      //
      // old values for the unknowns
      //

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
        if(it != residualValues.end())
          return it->second;
        DUNE_THROW(Exception, "Unknowns for time step " << step << " missing "
                   "in the cache!");
      }
      //! store a vector of unknowns  in the cache
      /**
       * \param step Step for which to store the vector of unknowns.
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

      ////////////////////////////////////////////////////////////////////////
      //
      // flushing
      //

      //! Flush old stuff but try to keep one Jacobian and one zero-residual
      /**
       * \param stepToKeep Step number of the oldest values to keep.
       *
       * Keeping at least one Jacobian and zero-residual is importand when the
       * multi-step methods are used to implement a one-step scheme: With only
       * one step, we would always flush the values we just stored in the
       * cache.  This is no problem for the residual values and the vector of
       * unknowns, since we can't reuse thos anyway.  However, flushing the
       * only value for the Jacobians and the zero-residuals inhibits caching
       * across time steps.
       */
      void flushOldKeepLatestAffine(Step stepToKeep) {
        for(std::size_t o = 0; o < jacobians.size(); ++o) {
          MatrixIterator it = jacobians[o].lower_bound(stepToKeep);
          if(it == jacobians[o].end() && it != jacobians[o].begin())
            --it;
          jacobians[o].erase(jacobians[o].begin(), it);
        }
        for(std::size_t o = 0; o < zeroResiduals.size(); ++o) {
          ResidualIterator it = zeroResiduals[o].lower_bound(stepToKeep);
          if(it == zeroResiduals[o].end() && it != zeroResiduals[o].begin())
            --it;
          zeroResiduals[o].erase(zeroResiduals[o].begin(), it);
        }
        for(std::size_t o = 0; o < residualValues.size(); ++o)
          residualValues[o].erase(residualValues[o].begin(),
                                  residualValues[o].lower_bound(stepToKeep));
        unknowns.erase(unknowns.begin(),
                       unknowns.lower_bound(stepToKeep));
      }

      //! Flush all cached values
      /**
       * This is useful for instance after adaption.  It is equivalent to
       * recreating the cache.
       */
      void flushAll() {
        jacobians.clear();
        zeroResiduals.clear();
        residualValues.clear();
        unknowns.clear();
      }

    };

    //! \} group MultiStepMethods
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_MULTISTEP_CACHE_HH
