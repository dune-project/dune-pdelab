// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_LFSINDEXCACHE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_LFSINDEXCACHE_HH

#include <vector>
#include <stack>
#include <algorithm>
#include <unordered_map>

#include <dune/common/reservedvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/hash.hh>
#include <dune/common/iteratorfacades.hh>

#include <dune/typetree/typetree.hh>

#include <dune/pdelab/constraints/common/constraintstransformation.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>

namespace Dune {
  namespace PDELab {

    template<typename Iterator>
    class DOFIndexViewIterator
      : public RandomAccessIteratorFacade<DOFIndexViewIterator<Iterator>,
                                          const typename std::iterator_traits<Iterator>::value_type::View,
                                          const typename std::iterator_traits<Iterator>::value_type::View
                                          >
    {

      friend class RandomAccessIteratorFacade<
        DOFIndexViewIterator,
        const typename std::iterator_traits<Iterator>::value_type::View,
        const typename std::iterator_traits<Iterator>::value_type::View
        >;

      typedef typename std::iterator_traits<Iterator>::value_type::View View;

    public:

      // Add support for returning non-references from iterator.
      // We need a little bit of magic to make operator->() work for this iterator
      // because we return a temporary object from dereference(), and the standard
      // implementation of operator->() in the facade tries to take the address of
      // that temporary, which the compiler will vehemently object to... ;-)
      //
      // So I borrowed the following neat little trick from Boost's iterator library:
      // The proxy object stores a copy of the temporary View object, and operator()->
      // returns the proxy object to the caller. As mandated by the standard, the compiler
      // will then attempt to repeat the operator->() on the returned object and get the
      // address of the copy stored in the (temporary) proxy object. That proxy object
      // is guaranteed to live until the next sequence point, and that is precisely as
      // long as we have to guarantee the validity of the pointer to our View object.
      // Problem solved - and another example of how difficult it is to get this low-level
      // stuff implemented on the same level as Boost...
      struct proxy
      {

        explicit proxy(const View& v)
          : _tmp(v)
        {}

        View* operator->()
        {
          return &_tmp;
        }

        View _tmp;
      };

      // The proxy object will stand in as a pointer
      typedef proxy pointer;

      DOFIndexViewIterator()
        : _iterator()
        , _tail_length(0)
      {}

      explicit DOFIndexViewIterator(Iterator it, std::size_t tail_length = 0)
        : _iterator(it)
        , _tail_length(tail_length)
      {}

      void cut_back()
      {
        ++_tail_length;
      }

      void restore_back()
      {
        --_tail_length;
      }

      const typename std::iterator_traits<Iterator>::reference raw_index() const
      {
        return *_iterator;
      }

      bool equals(const DOFIndexViewIterator& other) const
      {
        return _iterator == other._iterator;
      }

      void increment()
      {
        ++_iterator;
      }

      void decrement()
      {
        --_iterator;
      }

      void advance(int n)
      {
        _iterator += n;
      }

      std::ptrdiff_t distanceTo(DOFIndexViewIterator& other) const
      {
        return other._iterator - _iterator;
      }

      const View dereference() const
      {
        return _iterator->view(_iterator->treeIndex().size() - _tail_length);
      }

      pointer operator->() const
      {
        return pointer(dereference());
      }

    private:

      Iterator _iterator;
      std::size_t _tail_length;

    };


    template<typename Iterator>
    struct extract_lfs_leaf_size_visitor
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename LeafLFS, typename TreePath>
      void leaf(const LeafLFS& leaf_lfs, TreePath tp)
      {
        (*it) = leaf_lfs.size();
        ++it;
      }

      extract_lfs_leaf_size_visitor(Iterator leaf_size_container_iterator)
        : it(leaf_size_container_iterator)
      {}

      Iterator it;

    };

    template<typename LFS, typename Iterator>
    Iterator extract_lfs_leaf_sizes(const LFS& lfs, Iterator it)
    {
      extract_lfs_leaf_size_visitor<Iterator> visitor(it);
      TypeTree::applyToTree(lfs,visitor);
      return visitor.it;
    }


    template<typename DOFIterator,
             typename ContainerIterator,
             typename LeafSizeIterator,
             std::size_t tree_depth,
             bool fast>
    struct map_dof_indices_to_container_indices
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename Ordering, typename TreePath>
      void leaf(const Ordering& ordering, TreePath tp)
      {
        std::size_t leaf_size = *(leaf_size_pos++);
        if (fast)
          dof_end += 1;
        else
          dof_end += leaf_size;
        ordering.map_lfs_indices(dof_pos,dof_end,container_pos);
        dof_pos = dof_end;
        container_pos += leaf_size;
      }

      template<typename Ordering, typename TreePath>
      void post(const Ordering& ordering, TreePath tp)
      {
        if (Ordering::consume_tree_index)
          {
            dof_pos.restore_back();
            dof_end.restore_back();
          }
        ordering.map_lfs_indices(dof_stack.top(),dof_end,container_stack.top());
        dof_stack.pop();
        container_stack.pop();
      }

      template<typename Ordering, typename TreePath>
      void pre(const Ordering& ordering, TreePath tp)
      {
        dof_stack.push(dof_pos);
        container_stack.push(container_pos);
        if (Ordering::consume_tree_index)
          {
            dof_pos.cut_back();
            dof_end.cut_back();
          }
      }

      map_dof_indices_to_container_indices(DOFIterator dof_begin,
                                           ContainerIterator container_begin,
                                           LeafSizeIterator leaf_size_begin,
                                           std::size_t dof_index_tail_length = 0)
        : dof_pos(dof_begin,dof_index_tail_length)
        , dof_end(dof_begin,dof_index_tail_length)
        , container_pos(container_begin)
        , leaf_size_pos(leaf_size_begin)
      {}


      DOFIndexViewIterator<DOFIterator> dof_pos;
      DOFIndexViewIterator<DOFIterator> dof_end;
      ContainerIterator container_pos;
      LeafSizeIterator leaf_size_pos;
      std::stack<DOFIndexViewIterator<DOFIterator>,ReservedVector<DOFIndexViewIterator<DOFIterator>,tree_depth> > dof_stack;
      std::stack<ContainerIterator,ReservedVector<ContainerIterator,tree_depth> > container_stack;

    };



    template<typename LFS, typename C, typename CacheTag, bool fast>
    class LFSIndexCacheBase
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
      typedef typename Ordering::Traits::ContainerIndex ContainerIndex;
      typedef ContainerIndex CI;
      typedef typename Ordering::Traits::DOFIndex DOFIndex;
      typedef DOFIndex DI;
      typedef std::size_t size_type;

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
        if(fast) {
          _container_indices[0].resize(2);
          _container_indices[0][0] = 0;
          _container_indices[0][1] = LFS::Traits::GridFunctionSpace::Ordering::Traits::DOFIndexAccessor::entityIndex(_lfs.dofIndex(0));
        }
        else {

          // clear out existing state
          _container_index_map.clear();
          for (typename CIVector::iterator it = _container_indices.begin(); it != _container_indices.end(); ++it)
            it->clear();

          _inverse_map.clear();
          _inverse_cache_built = false;

          // extract size for all leaf spaces (into a flat list)
          typedef ReservedVector<size_type,TypeTree::TreeInfo<LFS>::leafCount> LeafSizeVector;
          LeafSizeVector leaf_sizes;
          leaf_sizes.resize(TypeTree::TreeInfo<LFS>::leafCount);
          extract_lfs_leaf_sizes(_lfs,leaf_sizes.begin());

          // perform the actual mapping
          map_dof_indices_to_container_indices<
            typename LFS::Traits::DOFIndexContainer::const_iterator,
            typename CIVector::iterator,
            typename LeafSizeVector::const_iterator,
            TypeTree::TreeInfo<Ordering>::depth,
            fast
            > index_mapper(_lfs._dof_indices->begin(),_container_indices.begin(),leaf_sizes.begin(),_lfs.subSpaceDepth());
          TypeTree::applyToTree(_lfs.gridFunctionSpace().ordering(),index_mapper);

          if (_enable_constraints_caching)
            {
              _constraints.resize(0);
              std::vector<std::pair<size_type,typename C::const_iterator> > non_dirichlet_constrained_dofs;
              size_type constraint_entry_count = 0;
              size_type end = fast ? 1 : _lfs.size();
              for (size_type i = 0; i < end; ++i)
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
        }
      }

      const DI& dofIndex(size_type i) const
      {
        return _lfs.dofIndex(i);
      }

      const CI& containerIndex(size_type i) const
      {
        return _container_indices[i];
      }

      const CI& containerIndex(const DI& i) const
      {
        // look up DOFIndex i
        std::pair<typename CIMap::iterator,bool> r = _container_index_map.insert(std::make_pair(std::ref(i),CI()));

        // i did not exist in the cache, map it into the newly inserted container index
        if (r.second)
            _lfs.gridFunctionSpace().ordering().mapIndex(i.view(),r.first->second);

        // return cached container index
        return r.first->second;
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

      bool constraintsCachingEnabled() const
      {
        return _enable_constraints_caching;
      }

    private:

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

      const LFS& _lfs;
      const bool _enable_constraints_caching;
      CIVector _container_indices;
      std::vector<unsigned char> _dof_flags;
      std::vector<std::pair<ConstraintsIterator,ConstraintsIterator> > _constraints_iterators;
      mutable CIMap _container_index_map;
      ConstraintsVector _constraints;
      mutable std::array<size_type,TypeTree::StaticDegree<LFS>::value> _offsets;
      mutable std::array<size_type,TypeTree::StaticDegree<LFS>::value> _extended_offsets;
      mutable bool _inverse_cache_built;
      mutable InverseMap _inverse_map;

      const C& _gfs_constraints;

    };


    template<typename LFS, typename CacheTag, bool fast>
    class LFSIndexCacheBase<LFS,EmptyTransformation,CacheTag,fast>
    {

    public:

      typedef LFS LocalFunctionSpace;
      typedef typename LFS::Traits::GridFunctionSpace GFS;
      typedef typename GFS::Ordering Ordering;
      typedef typename Ordering::Traits::ContainerIndex ContainerIndex;
      typedef ContainerIndex CI;
      typedef typename Ordering::Traits::DOFIndex DOFIndex;
      typedef DOFIndex DI;
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
        , _container_indices(lfs.maxSize())
      {
      }

      template<typename C>
      LFSIndexCacheBase(const LFS& lfs, const C& c, bool enable_constraints_caching)
        : _lfs(lfs)
        , _container_indices(lfs.maxSize())
      {
      }


      void update()
      {
        if(fast) {
          _container_indices[0].resize(2);
          _container_indices[0][0] = 0;
          _container_indices[0][1] = LFS::Traits::GridFunctionSpace::Ordering::Traits::DOFIndexAccessor::entityIndex(_lfs.dofIndex(0));
        }
        else {

          // clear out existing state
          _container_index_map.clear();
          for (typename CIVector::iterator it = _container_indices.begin(); it != _container_indices.end(); ++it)
            it->clear();

          // extract size for all leaf spaces (into a flat list)
          typedef ReservedVector<size_type,TypeTree::TreeInfo<LFS>::leafCount> LeafSizeVector;
          LeafSizeVector leaf_sizes;
          leaf_sizes.resize(TypeTree::TreeInfo<LFS>::leafCount);
          extract_lfs_leaf_sizes(_lfs,leaf_sizes.begin());

          // perform the actual mapping
          map_dof_indices_to_container_indices<
            typename LFS::Traits::DOFIndexContainer::const_iterator,
            typename CIVector::iterator,
            typename LeafSizeVector::const_iterator,
            TypeTree::TreeInfo<Ordering>::depth,
            fast
            > index_mapper(_lfs._dof_indices->begin(),_container_indices.begin(),leaf_sizes.begin(),_lfs.subSpaceDepth());
          TypeTree::applyToTree(_lfs.gridFunctionSpace().ordering(),index_mapper);
        }
      }

      const DI& dofIndex(size_type i) const
      {
        return _lfs.dofIndex(i);
      }

      const CI& containerIndex(size_type i) const
      {
        return _container_indices[i];
      }

      const CI& containerIndex(const DI& i) const
      {
        // look up DOFIndex i
        std::pair<typename CIMap::iterator,bool> r = _container_index_map.insert(std::make_pair(std::ref(i),CI()));

        // i did not exist in the cache, map it into the newly inserted container index
        if (r.second)
            _lfs.gridFunctionSpace().ordering().mapIndex(i.view(),r.first->second);

        // return cached container index
        return r.first->second;
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

    private:

      const LFS& _lfs;
      CIVector _container_indices;
      mutable CIMap _container_index_map;
      const ConstraintsVector _constraints;

    };



    template<typename LFS, typename C, bool fast>
    class LFSIndexCacheBase<LFS,C,SimpleLFSCacheTag,fast>
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
        if(fast) {
          DUNE_THROW(Dune::Exception, "This function shouldn't be called in fast mode");
        }
        else {

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


    template<typename LFS, bool fast>
    class LFSIndexCacheBase<LFS,EmptyTransformation,SimpleLFSCacheTag,fast>
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


    template<typename LFS, typename C = EmptyTransformation, bool fast = false>
    class LFSIndexCache
      : public LFSIndexCacheBase<LFS,C,typename LFS::Traits::GridFunctionSpace::Ordering::CacheTag,fast>
    {

    public:

      template<typename CC>
      LFSIndexCache(const LFS& lfs, const CC& c, bool enable_constraints_caching = !std::is_same<C,EmptyTransformation>::value)
        : LFSIndexCacheBase<LFS,C,typename LFS::Traits::GridFunctionSpace::Ordering::CacheTag,fast>(lfs,c,enable_constraints_caching)
      {
      }

      explicit LFSIndexCache(const LFS& lfs)
        : LFSIndexCacheBase<LFS,C,typename LFS::Traits::GridFunctionSpace::Ordering::CacheTag,fast>(lfs)
      {
      }

    };


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_LFSINDEXCACHE_HH
