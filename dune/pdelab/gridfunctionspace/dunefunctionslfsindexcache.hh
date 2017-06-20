// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_DUNEFUNCTIONSLFSINDEXCACHE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_DUNEFUNCTIONSLFSINDEXCACHE_HH

#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/dunefunctionsgridfunctionspace.hh>

namespace Dune {
  namespace PDELab {

    template<typename LFS, typename C>
    class LFSIndexCacheBase<LFS,C,Experimental::DuneFunctionsCacheTag, false>
    {

      enum DOFFlags
        {
          DOF_NONCONSTRAINED = 0,
          DOF_CONSTRAINED = 1<<0,
          DOF_DIRICHLET = 1<<1
        };

    public:

      using LocalFunctionSpace = LFS;
      using GFS                = typename LFS::Traits::GridFunctionSpace;
      using Ordering           = typename GFS::Ordering;
      using DOFIndex           = typename Ordering::Traits::DOFIndex;
      using DI                 = DOFIndex;
      using ContainerIndex     = typename Ordering::Traits::ContainerIndex;
      using CI                 = ContainerIndex;
      using size_type          = std::size_t;

      typedef std::vector<CI> CIVector;
      typedef std::unordered_map<DI,CI> CIMap;

      typedef std::unordered_map<const CI*,std::pair<size_type,bool> > InverseMap;

      struct ConstraintsEntry
        : public std::pair<const CI*,typename C::mapped_type::mapped_type>
      {
        typedef CI ContainerIndex;
        typedef typename C::mapped_type::mapped_type Weight;

        const ContainerIndex& containerIndex() const
        {
          return *(this->first);
        }

        const Weight& weight() const
        {
          return this->second;
        }
      };

      //typedef std::pair<CI,typename C::mapped_type::mapped_type> ConstraintsEntry;

      typedef std::vector<ConstraintsEntry> ConstraintsVector;
      typedef typename ConstraintsVector::const_iterator ConstraintsIterator;

      LFSIndexCacheBase(const LFS& lfs, const C& constraints, bool enable_constraints_caching)
        : _lfs(lfs)
        , _enable_constraints_caching(enable_constraints_caching)
        , _container_indices(lfs.maxSize())
        , _dof_flags(lfs.maxSize(),0)
        , _constraints_iterators(lfs.maxSize())
        , _inverse_cache_built(false)
        , _gfs_constraints(constraints)
      {
      }

      void update()
      {}

      DI dofIndex(size_type i) const
      {
        return _lfs.dofIndex(i);
      }

      CI containerIndex(size_type i) const
      {
        return _lfs.dofIndex(i);
      }

      CI containerIndex(const DI& i) const
      {
        return i;
      }

      bool isConstrained(size_type i) const
      {
        return _dof_flags[i] & DOF_CONSTRAINED;
      }

      bool isDirichletConstraint(size_type i) const
      {
        return _dof_flags[i] & DOF_DIRICHLET;
      }

      ConstraintsIterator constraintsBegin(size_type i) const
      {
        assert(isConstrained(i));
        return _constraints_iterators[i].first;
      }

      ConstraintsIterator constraintsEnd(size_type i) const
      {
        assert(isConstrained(i));
        return _constraints_iterators[i].second;
      }

      const LocalFunctionSpace& localFunctionSpace() const
      {
        return _lfs;
      }

      size_type size() const
      {
        return _lfs.size();
      }

      bool constraintsCachingEnabled() const
      {
        return _enable_constraints_caching;
      }

    private:

      const LFS& _lfs;
      const bool _enable_constraints_caching;
      CIVector _container_indices;
      std::vector<unsigned char> _dof_flags;
      std::vector<std::pair<ConstraintsIterator,ConstraintsIterator> > _constraints_iterators;
      mutable CIMap _container_index_map;
      ConstraintsVector _constraints;
      mutable bool _inverse_cache_built;
      mutable InverseMap _inverse_map;

      const C& _gfs_constraints;

    };


    template<typename LFS>
    class LFSIndexCacheBase<LFS,EmptyTransformation,Experimental::DuneFunctionsCacheTag, false>
    {

    public:

      using LocalFunctionSpace = LFS;
      using GFS                = typename LFS::Traits::GridFunctionSpace;
      using Ordering           = typename GFS::Ordering;
      using DOFIndex           = typename Ordering::Traits::DOFIndex;
      using DI                 = DOFIndex;
      using ContainerIndex     = typename Ordering::Traits::ContainerIndex;
      using CI                 = ContainerIndex;
      using size_type          = std::size_t;

      struct ConstraintsEntry
        : public std::pair<const CI*,double>
      {
        typedef CI ContainerIndex;
        typedef double Weight;

        const ContainerIndex& containerIndex() const
        {
          return *(this->first);
        }

        const Weight& weight() const
        {
          return this->second;
        }
      };

      typedef std::vector<ConstraintsEntry> ConstraintsVector;
      typedef typename ConstraintsVector::const_iterator ConstraintsIterator;

      typedef std::vector<CI> CIVector;
      typedef std::unordered_map<DI,CI> CIMap;

      explicit LFSIndexCacheBase(const LFS& lfs)
        : _lfs(lfs)
      {}

      template<typename C>
      LFSIndexCacheBase(const LFS& lfs, const C& c, bool enable_constraints_caching)
        : _lfs(lfs)
      {}

      DI dofIndex(size_type i) const
      {
        return _lfs.dofIndex(i);
      }

      CI containerIndex(size_type i) const
      {
        return _lfs.dofIndex(i);
      }

      CI containerIndex(const DI& i) const
      {
        return i;
      }

      bool isConstrained(size_type i) const
      {
        return false;
      }

      bool isDirichletConstraint(size_type i) const
      {
        return false;
      }

      ConstraintsIterator constraintsBegin(size_type i) const
      {
        return _constraints.begin();
      }

      ConstraintsIterator constraintsEnd(size_type i) const
      {
        return _constraints.end();
      }

      const LocalFunctionSpace& localFunctionSpace() const
      {
        return _lfs;
      }

      size_type size() const
      {
        return _lfs.size();
      }

      bool constraintsCachingEnabled() const
      {
        return false;
      }

      void update()
      {}

    private:

      const LFS& _lfs;
      CIVector _container_indices;
      mutable CIMap _container_index_map;
      const ConstraintsVector _constraints;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_DUNEFUNCTIONSLFSINDEXCACHE_HH
