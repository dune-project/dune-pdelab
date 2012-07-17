// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_LFSCONTAINERINDEXCACHE_HH
#define DUNE_PDELAB_LFSCONTAINERINDEXCACHE_HH

#include <vector>
#include <stack>
#include <unordered_map>

#include <dune/common/reservedvector.hh>

#include <dune/pdelab/common/typetree.hh>
#include <dune/pdelab/constraints/constraintstransformation.hh>


#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/functional/hash.hpp>

namespace Dune {
  namespace PDELab {

    template<typename Iterator>
    class DOFIndexViewIterator
      : public boost::iterator_adaptor<DOFIndexViewIterator<Iterator>,
                                       Iterator,
                                       const typename std::iterator_traits<Iterator>::value_type::View,
                                       boost::use_default,
                                       const typename std::iterator_traits<Iterator>::value_type::View
                                       >
    {

    public:

      DOFIndexViewIterator()
        : DOFIndexViewIterator::iterator_adaptor_()
        , _tail_length(0)
      {}

      explicit DOFIndexViewIterator(Iterator it)
        : DOFIndexViewIterator::iterator_adaptor_(it)
        , _tail_length(0)
      {}

      void cut_back()
      {
        ++_tail_length;
      }

      void restore_back()
      {
        --_tail_length;
      }

    private:

      friend class boost::iterator_core_access;
      typedef typename std::iterator_traits<Iterator>::value_type::View View;

      std::size_t _tail_length;

      const View dereference() const
      {
        return this->base()->view(this->base()->treeIndex().size() - _tail_length);
      }

    };


    template<typename Container>
    struct extract_lfs_leaf_sizes
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename LeafLFS, typename TreePath>
      void leaf(const LeafLFS& leaf_lfs, TreePath tp)
      {
        leaf_sizes.push_back(leaf_lfs.size());
      }

      extract_lfs_leaf_sizes(Container& leaf_sizes_)
        : leaf_sizes(leaf_sizes_)
      {}

      Container& leaf_sizes;

    };


    template<typename DOFIterator,
             typename ContainerIterator,
             typename LeafSizeIterator,
             std::size_t tree_depth>
    struct map_dof_indices_to_container_indices
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename Ordering, typename TreePath>
      void leaf(const Ordering& ordering, TreePath tp)
      {
        std::size_t leaf_size = *(leaf_size_pos++);
        dof_end += leaf_size;
        ordering.map_indices(dof_pos,dof_end,container_pos);
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
        ordering.map_indices(dof_stack.top(),dof_end,container_stack.top());
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
                                           LeafSizeIterator leaf_size_begin)
        : dof_pos(dof_begin)
        , dof_end(dof_begin)
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

    template<typename LFS, typename C, typename CacheTag>
    class LFSContainerIndexCacheBase
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
      typedef std::unordered_map<DI,CI,boost::hash<DI> > CIMap;

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

      LFSContainerIndexCacheBase(const LFS& lfs, const C& constraints)
        : _lfs(lfs)
        , _container_indices(lfs.maxSize())
        , _dof_flags(lfs.maxSize())
        , _constraints_iterators(lfs.maxSize())
        , _gfs_constraints(constraints)
      {
      }

      void update()
      {
        // clear out existing state
        _container_index_map.clear();
        for (typename CIVector::iterator it = _container_indices.begin(); it != _container_indices.end(); ++it)
          it->clear();

        // extract size for all leaf spaces (into a flat list)
        typedef ReservedVector<size_type,TypeTree::TreeInfo<LFS>::leafCount> LeafSizeVector;
        LeafSizeVector leaf_sizes;
        extract_lfs_leaf_sizes<LeafSizeVector> leaf_size_extractor(leaf_sizes);
        TypeTree::applyToTree(_lfs,leaf_size_extractor);

        // perform the actual mapping
        map_dof_indices_to_container_indices<
          typename LFS::Traits::DOFIndexContainer::const_iterator,
          typename CIVector::iterator,
          typename LeafSizeVector::const_iterator,
          TypeTree::TreeInfo<Ordering>::depth
          > index_mapper(_lfs._dof_indices->begin(),_container_indices.begin(),leaf_sizes.begin());
        TypeTree::applyToTree(_lfs.gridFunctionSpace().ordering(),index_mapper);

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

      const DI& dof_index(size_type i) const
      {
        return _lfs.dofIndex(i);
      }

      const CI& container_index(size_type i) const
      {
        return _container_indices[i];
      }

      const CI& container_index(const DI& i) const
      {
        // look up DOFIndex i
        std::pair<typename CIMap::iterator,bool> r = _container_index_map.insert(std::make_pair(std::ref(i),CI()));

        // i did not exist in the cache, map it into the newly inserted container index
        if (r.second)
            _lfs.gridFunctionSpace().ordering().map_index(i.view(),r.first->second);

        // return cached container index
        return r.first->second;
      }

      bool constrained(size_type i) const
      {
        return _dof_flags[i] & DOF_CONSTRAINED;
      }

      bool dirichlet_constraint(size_type i) const
      {
        return _dof_flags[i] & DOF_DIRICHLET;
      }

      ConstraintsIterator constraints_begin(size_type i) const
      {
        assert(constrained(i));
        return _constraints_iterators[i].first;
      }

      ConstraintsIterator constraints_end(size_type i) const
      {
        assert(constrained(i));
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


    template<typename LFS, typename CacheTag>
    class LFSContainerIndexCacheBase<LFS,EmptyTransformation,CacheTag>
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
      typedef std::unordered_map<DI,CI,boost::hash<DI> > CIMap;

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

      explicit LFSContainerIndexCacheBase(const LFS& lfs)
        : _lfs(lfs)
        , _container_indices(lfs.maxSize())
      {
      }

      template<typename C>
      LFSContainerIndexCacheBase(const LFS& lfs, const C& c)
        : _lfs(lfs)
        , _container_indices(lfs.maxSize())
      {
      }


      void update()
      {
        // clear out existing state
        _container_index_map.clear();
        for (typename CIVector::iterator it = _container_indices.begin(); it != _container_indices.end(); ++it)
          it->clear();

        // extract size for all leaf spaces (into a flat list)
        typedef ReservedVector<size_type,TypeTree::TreeInfo<LFS>::leafCount> LeafSizeVector;
        LeafSizeVector leaf_sizes;
        extract_lfs_leaf_sizes<LeafSizeVector> leaf_size_extractor(leaf_sizes);
        TypeTree::applyToTree(_lfs,leaf_size_extractor);

        // perform the actual mapping
        map_dof_indices_to_container_indices<
          typename LFS::Traits::DOFIndexContainer::const_iterator,
          typename CIVector::iterator,
          typename LeafSizeVector::const_iterator,
          TypeTree::TreeInfo<Ordering>::depth
          > index_mapper(_lfs._dof_indices->begin(),_container_indices.begin(),leaf_sizes.begin());
        TypeTree::applyToTree(_lfs.gridFunctionSpace().ordering(),index_mapper);
      }

      const CI& container_index(size_type i) const
      {
        return _container_indices[i];
      }

      const CI& container_index(const DI& i) const
      {
        // look up DOFIndex i
        std::pair<typename CIMap::iterator,bool> r = _container_index_map.insert(std::make_pair(std::ref(i),CI()));

        // i did not exist in the cache, map it into the newly inserted container index
        if (r.second)
            _lfs.gridFunctionSpace().ordering().map_index(i.view(),r.first->second);

        // return cached container index
        return r.first->second;
      }

      bool constrained(size_type i) const
      {
        return false;
      }

      bool dirichlet_constraint(size_type i) const
      {
        return false;
      }

      ConstraintsIterator constraints_begin(size_type i) const
      {
        return _constraints.begin();
      }

      ConstraintsIterator constraints_end(size_type i) const
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
      CIVector _container_indices;
      mutable CIMap _container_index_map;
      const ConstraintsVector _constraints;

    };



    template<typename LFS, typename C>
    class LFSContainerIndexCacheBase<LFS,C,SimpleLFSCacheTag>
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
      typedef std::unordered_map<DI,CI,boost::hash<DI> > CIMap;

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

      LFSContainerIndexCacheBase(const LFS& lfs, const C& constraints)
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

      const DI& dof_index(size_type i) const
      {
        return _lfs.dofIndex(i);
      }

      CI container_index(size_type i) const
      {
        return CI(_lfs.dofIndex(i)[0]);
      }

      const CI& container_index(const DI& i) const
      {
        return CI(i[0]);
      }

      bool constrained(size_type i) const
      {
        return _dof_flags[i] & DOF_CONSTRAINED;
      }

      bool dirichlet_constraint(size_type i) const
      {
        return _dof_flags[i] & DOF_DIRICHLET;
      }

      ConstraintsIterator constraints_begin(size_type i) const
      {
        assert(constrained(i));
        return _constraints_iterators[i].first;
      }

      ConstraintsIterator constraints_end(size_type i) const
      {
        assert(constrained(i));
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
    class LFSContainerIndexCacheBase<LFS,EmptyTransformation,SimpleLFSCacheTag>
    {

    public:

      typedef LFS LocalFunctionSpace;
      typedef typename LFS::Traits::GridFunctionSpace GFS;
      typedef typename GFS::Ordering Ordering;
      typedef typename Ordering::Traits::ContainerIndex CI;
      typedef typename Ordering::Traits::DOFIndex DI;
      typedef std::size_t size_type;

      typedef std::vector<CI> CIVector;
      typedef std::unordered_map<DI,CI,boost::hash<DI> > CIMap;

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

      explicit LFSContainerIndexCacheBase(const LFS& lfs)
        : _lfs(lfs)
      {
      }

      template<typename C>
      LFSContainerIndexCacheBase(const LFS& lfs, const C& c)
        : _lfs(lfs)
      {
      }


      void update()
      {
        // there's nothing to do here...
      }

      CI container_index(size_type i) const
      {
        return CI(_lfs.dofIndex(i)[0]);
      }

      CI container_index(const DI& i) const
      {
        return CI(i[0]);
      }

      bool constrained(size_type i) const
      {
        return false;
      }

      bool dirichlet_constraint(size_type i) const
      {
        return false;
      }

      ConstraintsIterator constraints_begin(size_type i) const
      {
        return _constraints.begin();
      }

      ConstraintsIterator constraints_end(size_type i) const
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
    class LFSContainerIndexCache
      : public LFSContainerIndexCacheBase<LFS,C,typename LFS::Traits::GridFunctionSpace::Ordering::CacheTag>
    {

    public:

      template<typename CC>
      LFSContainerIndexCache(const LFS& lfs, const CC& c)
        : LFSContainerIndexCacheBase<LFS,C,typename LFS::Traits::GridFunctionSpace::Ordering::CacheTag>(lfs,c)
      {
      }

      LFSContainerIndexCache(const LFS& lfs)
        : LFSContainerIndexCacheBase<LFS,C,typename LFS::Traits::GridFunctionSpace::Ordering::CacheTag>(lfs)
      {
      }

    };


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LFSCONTAINERINDEXCACHE_HH
