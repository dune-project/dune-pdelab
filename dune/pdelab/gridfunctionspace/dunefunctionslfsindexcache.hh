// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_DUNEFUNCTIONSLFSINDEXCACHE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_DUNEFUNCTIONSLFSINDEXCACHE_HH

#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/dunefunctionsgridfunctionspace.hh>

namespace Dune {
  namespace PDELab {

    template<typename LFS, typename C>
    class LFSIndexCacheBase<LFS,C,Experimental::DuneFunctionsCacheTag>
    {

      enum DOFFlags
        {
          DOF_NONCONSTRAINED = 0,
          DOF_CONSTRAINED = 1<<0,
          DOF_DIRICHLET = 1<<1
        };

    public:

      using LocalFunctionSpace = LFS;
      using GFS                = typename LFS::Traits::GridFunctionSpace;;
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
      {
        /*
        if (_enable_constraints_caching)
          {
            _constraints.resize(0);
            std::vector<std::pair<size_type,typename C::const_iterator> > non_dirichlet_constrained_dofs;
            size_type constraint_entry_count = 0;
            for (size_type i = 0; i < _lfs.size(); ++i)
              {
                const CI& container_index = _container_indices[i];
                const typename C::const_iterator cit = _gfs_constraints.find(container_index);
                if (cit == _gfs_constraints.end())
                  {
                    _dof_flags[i] = DOF_NONCONSTRAINED;
                    continue;
                  }

                if (cit->second.size() == 0)
                  {
                    _dof_flags[i] = DOF_CONSTRAINED | DOF_DIRICHLET;
                    _constraints_iterators[i] = make_pair(_constraints.end(),_constraints.end());
                  }
                else
                  {
                    _dof_flags[i] = DOF_CONSTRAINED;
                    constraint_entry_count += cit->second.size();
                    non_dirichlet_constrained_dofs.push_back(make_pair(i,cit));
                  }
              }

            if (constraint_entry_count > 0)
              {
                _constraints.resize(constraint_entry_count);
                typename ConstraintsVector::iterator eit = _constraints.begin();
                for (typename std::vector<std::pair<size_type,typename C::const_iterator> >::const_iterator it = non_dirichlet_constrained_dofs.begin();
                     it != non_dirichlet_constrained_dofs.end();
                     ++it)
                  {
                    _constraints_iterators[it->first].first = eit;
                    for (typename C::mapped_type::const_iterator cit = it->second->second.begin(); cit != it->second->second.end(); ++cit, ++eit)
                      {
                        eit->first = &(cit->first);
                        eit->second = cit->second;
                      }
                    _constraints_iterators[it->first].second = eit;
                  }
              }
          }
        */
      }

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

      /*
      std::pair<size_type,bool> localIndex(const ContainerIndex& ci) const
      {
        if (!_inverse_cache_built)
          build_inverse_cache();
        return _inverse_map[ci];
      }

      size_type offset(size_type i) const
      {
        if (!_inverse_cache_built)
          build_inverse_cache();
        return _offsets[i];
      }

      size_type extendedOffset(size_type i) const
      {
        if (!_inverse_cache_built)
          build_inverse_cache();
        return _extended_offsets[i];
      }
      */

      bool constraintsCachingEnabled() const
      {
        return _enable_constraints_caching;
      }

    private:

      /*
      struct sort_container_indices
      {
        template<typename T>
        bool operator()(const T* a, const T* b) const
        {
          return std::lexicographical_compare(reversed_iterator(a->end()),reversed_iterator(a->begin()),
                                              reversed_iterator(b->end()),reversed_iterator(b->begin())
                                              );
        }
      };
      */


      /*
      void build_inverse_cache()
      {
        size_type i = 0;
        size_type child = 0;
        _offsets[0] = 0;
        for (typename CIVector::const_iterator it = _container_indices.begin(),
               endit = _container_indices.end();
             it != endit;
             ++it, ++i
             )
          {
            _inverse_map.insert(std::make_pair(&(*it),std::make_pair(i,false)));
            if (it->back() != child)
              {
                _offsets[child+1] = i;
                ++child;
              }
          }

        std::vector<const ContainerIndex*> extended_cis;
        extended_cis.reserve(_constraints.size());

        for (typename ConstraintsVector::const_iterator it = _constraints.begin(),
               endit = _constraints.end();
             it != endit;
             ++it
             )
          {
            if (_inverse_map.count(it->first) == 0)
              extended_cis.push_back(it->first);
          }

        std::sort(extended_cis.begin(),extended_cis.end(),sort_container_indices());

        typename std::vector<const ContainerIndex*>::const_iterator endit = std::unique(extended_cis.begin(),extended_cis.end());

        i = 0;
        child = 0;
        for (typename std::vector<const ContainerIndex*>::const_iterator it = extended_cis.begin(); it != endit; ++it, ++i)
          {
            _inverse_map.insert(std::make_pair(&(*it),std::make_pair(i,true)));
            if (it->back() != child)
              {
                _extended_offsets[child+1] = i;
                ++child;
              }
          }

        _inverse_cache_built = true;

      }
      */

      const LFS& _lfs;
      const bool _enable_constraints_caching;
      CIVector _container_indices;
      std::vector<unsigned char> _dof_flags;
      std::vector<std::pair<ConstraintsIterator,ConstraintsIterator> > _constraints_iterators;
      mutable CIMap _container_index_map;
      ConstraintsVector _constraints;
      //mutable std::array<size_type,TypeTree::StaticDegree<LFS>::value> _offsets;
      //mutable std::array<size_type,TypeTree::StaticDegree<LFS>::value> _extended_offsets;
      mutable bool _inverse_cache_built;
      mutable InverseMap _inverse_map;

      const C& _gfs_constraints;

    };


    template<typename LFS>
    class LFSIndexCacheBase<LFS,EmptyTransformation,Experimental::DuneFunctionsCacheTag>
    {

    public:

      using LocalFunctionSpace = LFS;
      //using LocalFunctionSpace = Experimental::LocalFunctionSpace<GFS_>;
      // using LFS                = typename LFS::Traits::GridFunctionSpace;
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

    /*

    template<typename LFS, typename C>
    class LFSIndexCacheBase<LFS,C,SimpleLFSCacheTag>
    {

      enum DOFFlags
        {
          DOF_NONCONSTRAINED = 0,
          DOF_CONSTRAINED = 1<<0,
          DOF_DIRICHLET = 1<<1
        };

    public:

      typedef LFS LocalFunctionSpace;

      typedef typename LFS::Traits::GridFunctionSpace GFS;
      typedef typename GFS::Ordering Ordering;
      typedef typename Ordering::Traits::ContainerIndex CI;
      typedef typename Ordering::Traits::DOFIndex DI;
      typedef std::size_t size_type;

      typedef std::vector<CI> CIVector;
      typedef std::unordered_map<DI,CI> CIMap;

      struct ConstraintsEntry
        : public std::pair<CI,typename C::mapped_type::mapped_type>
      {
        typedef CI ContainerIndex;
        typedef typename C::mapped_type::mapped_type Weight;

        const ContainerIndex& containerIndex() const
        {
          return this->first;
        }

        const Weight& weight() const
        {
          return this->second;
        }
      };

      typedef std::vector<ConstraintsEntry> ConstraintsVector;
      typedef typename ConstraintsVector::const_iterator ConstraintsIterator;

      LFSIndexCacheBase(const LFS& lfs, const C& constraints)
        : _lfs(lfs)
        , _dof_flags(lfs.maxSize())
        , _constraints_iterators(lfs.maxSize())
        , _gfs_constraints(constraints)
      {
      }


      void update()
      {
        _constraints.resize(0);
        std::vector<std::pair<size_type,typename C::const_iterator> > non_dirichlet_constrained_dofs;
        size_type constraint_entry_count = 0;
        for (size_type i = 0; i < _lfs.size(); ++i)
          {
            const DI& dof_index = _lfs.dofIndex(i);
            const typename C::const_iterator cit = _gfs_constraints.find(dof_index);
            if (cit == _gfs_constraints.end())
              {
                _dof_flags[i] = DOF_NONCONSTRAINED;
                continue;
              }

            if (cit->second.size() == 0)
              {
                _dof_flags[i] = DOF_CONSTRAINED | DOF_DIRICHLET;
                _constraints_iterators[i] = make_pair(_constraints.end(),_constraints.end());
              }
            else
              {
                _dof_flags[i] = DOF_CONSTRAINED;
                constraint_entry_count += cit->second.size();
                non_dirichlet_constrained_dofs.push_back(make_pair(i,cit));
              }
          }

        if (constraint_entry_count > 0)
          {
            _constraints.resize(constraint_entry_count);
            typename ConstraintsVector::iterator eit = _constraints.begin();
            for (typename std::vector<std::pair<size_type,typename C::const_iterator> >::const_iterator it = non_dirichlet_constrained_dofs.begin();
                 it != non_dirichlet_constrained_dofs.end();
                 ++it)
              {
                _constraints_iterators[it->first].first = eit;
                for (typename C::mapped_type::const_iterator cit = it->second->second.begin(); cit != it->second->second.end(); ++cit, ++eit)
                  {
                    eit->first = cit->first;
                    eit->second = cit->second;
                  }
                _constraints_iterators[it->first].second = eit;
              }
          }

      }

      const DI& dofIndex(size_type i) const
      {
        return _lfs.dofIndex(i);
      }

      CI containerIndex(size_type i) const
      {
        return CI(_lfs.dofIndex(i)[0]);
      }

      const CI& containerIndex(const DI& i) const
      {
        return CI(i[0]);
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

    private:

      const LFS& _lfs;
      CIVector _container_indices;
      std::vector<unsigned char> _dof_flags;
      std::vector<std::pair<ConstraintsIterator,ConstraintsIterator> > _constraints_iterators;
      mutable CIMap _container_index_map;
      ConstraintsVector _constraints;

      const C& _gfs_constraints;

    };


    template<typename LFS>
    class LFSIndexCacheBase<LFS,EmptyTransformation,SimpleLFSCacheTag>
    {

    public:

      typedef LFS LocalFunctionSpace;
      typedef typename LFS::Traits::GridFunctionSpace GFS;
      typedef typename GFS::Ordering Ordering;
    private:
      typedef typename Ordering::Traits::ContainerIndex CI;
      typedef typename Ordering::Traits::DOFIndex DI;
    public:
      typedef CI ContainerIndex;
      typedef DI DOFIndex;
      typedef std::size_t size_type;

      typedef std::vector<CI> CIVector;
      typedef std::unordered_map<DI,CI> CIMap;

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

      explicit LFSIndexCacheBase(const LFS& lfs)
        : _lfs(lfs)
      {
      }

      template<typename C>
      LFSIndexCacheBase(const LFS& lfs, const C& c)
        : _lfs(lfs)
      {
      }


      void update()
      {
        // there's nothing to do here...
      }

      CI containerIndex(size_type i) const
      {
        return CI(_lfs.dofIndex(i)[0]);
      }

      CI containerIndex(const DI& i) const
      {
        return CI(i[0]);
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

    private:

      const LFS& _lfs;
      mutable CIMap _container_index_map;
      const ConstraintsVector _constraints;

    };


    template<typename LFS, typename C = EmptyTransformation>
    class LFSIndexCache
      : public LFSIndexCacheBase<LFS,C,typename LFS::Traits::GridFunctionSpace::Ordering::CacheTag>
    {

    public:

      template<typename CC>
      LFSIndexCache(const LFS& lfs, const CC& c, bool enable_constraints_caching = !std::is_same<C,EmptyTransformation>::value)
        : LFSIndexCacheBase<LFS,C,typename LFS::Traits::GridFunctionSpace::Ordering::CacheTag>(lfs,c,enable_constraints_caching)
      {
      }

      explicit LFSIndexCache(const LFS& lfs)
        : LFSIndexCacheBase<LFS,C,typename LFS::Traits::GridFunctionSpace::Ordering::CacheTag>(lfs)
      {
      }

    };
    */


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_DUNEFUNCTIONSLFSINDEXCACHE_HH
