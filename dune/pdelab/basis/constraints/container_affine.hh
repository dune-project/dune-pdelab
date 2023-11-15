#ifndef DUNE_PDELAB_BASIS_CONSTRAINTS_CONTAINER_AFFINE_HH
#define DUNE_PDELAB_BASIS_CONSTRAINTS_CONTAINER_AFFINE_HH

#include <dune/pdelab/concepts/multiindex.hh>
#include <dune/pdelab/concepts/basis.hh>

#include <dune/pdelab/common/multiindex.hh>
#include <dune/pdelab/common/entity_cast.hh>
#include <dune/pdelab/common/container_entry.hh>

#include <dune/grid/concepts/gridview.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/typetree/leafnode.hh>
#include <dune/typetree/treepath.hh>

#include <dune/common/float_cmp.hh>

#ifdef HAVE_TBB
#include <tbb/parallel_for.h>
#include <tbb/concurrent_unordered_map.h>
#endif

#if __cpp_lib_execution
#include <execution>
#include <algorithm>
#elif __has_include(<oneapi/dpl/execution>)
#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#endif
#include <tuple>
#include <span>
#include <unordered_map>
#include <vector>
#include <utility>

namespace Dune::PDELab::inline Experimental {


namespace Impl {

  // store constraints of the form Cx+b
  template<class Value, Concept::MultiIndex ContainerIndex>
  struct AffineConstraintsContainerBase {
    // thread-safe if HAVE_TBB is defined
    void addAffineConstraint(ContainerIndex ci, Value translation, std::span<const std::pair<ContainerIndex, double>> linear_coefficients) {
      // invariant: linear coefficents are always sorted lexicograpically by the multi-index
      // sanity checks
      for (auto [col, weight] : linear_coefficients) {
        if (FloatCmp::eq<double>(weight,0.))
          DUNE_THROW(RangeError, "Weights cannot be zero");
        if (col == ci)
          DUNE_THROW(RangeError, "Constraint cannot constraint to itself");
      }
      auto it = _map.find(ci);
      if (it == _map.end()) {
        auto vec = std::vector(begin(linear_coefficients), end(linear_coefficients));
        std::sort(begin(vec),end(vec), [](const auto& l, const auto& r){
          // (move to dynamic multi-index to make comparison easier)
          return Dune::PDELab::MultiIndex(l.first) < Dune::PDELab::MultiIndex(r.first);
        });
        _map.emplace(std::piecewise_construct,
            std::forward_as_tuple(ci),
            std::forward_as_tuple(translation, std::move(vec)));
      } else {
        auto& ltrans = it->second.first;
        if (FloatCmp::ne<Value>(ltrans, translation)) // repeated constraint with different translation
          DUNE_THROW(RangeError, "Constraint [" << ci << "] added twice with different translation");

        auto& lcoef = it->second.second;
        // ci is already constraints, we hace to add constraints one-by-one
        for (auto [cci, coeff] : linear_coefficients) {
          // find repeated constraints
          auto iit = std::lower_bound(begin(lcoef), end(lcoef), Dune::PDELab::MultiIndex(cci) ,[](const auto& l, const auto& r){
            // (move to dynamic multi-index to make comparison easier)
            return Dune::PDELab::MultiIndex(l.first) < r;
          });
          if ((iit == end(lcoef)) or (iit->first == cci)) // new constraint
            lcoef.insert(iit, {cci, coeff});
          else if (FloatCmp::ne<double>(iit->second, coeff)) // repeated constraint with different weight
            DUNE_THROW(RangeError, "Constraint [" << ci << "] -> [" << cci <<  "] added twice with different weight");
        }
      }
    }

    void addLinearConstraint(ContainerIndex ci, std::span<const std::pair<ContainerIndex, double>> linear_coefficients) {
      addAffineConstraint(ci, 0., linear_coefficients);
    }

    void addTranslationConstraint(ContainerIndex ci, Value translation) {
      addAffineConstraint(ci, translation, {});
    }


    [[nodiscard]] bool isConstrained(ContainerIndex ci) const {
      return (_map.count(ci) == 1);
    }

    void applyAffineConstraints(auto& x) const {
      const auto y = x;
      forEachSparseRow([&](std::size_t i) {
        auto [row, value, offset] = _global_offset[i];
        auto next_offset = std::get<2>(_global_offset[i+1]);
        for (; offset != next_offset; ++offset) {
          auto [col, weight] =  _global_constraints[offset];
          value += weight * containerEntry(y, col);
        }
        containerEntry(x, row) = value;
      });
    }

    void applyLinearConstraints(auto& x) const {
      const auto y = x;
      forEachSparseRow([&](std::size_t i) {
        auto [row, _, offset] = _global_offset[i];
        auto next_offset = std::get<2>(_global_offset[i+1]);
        containerEntry(x, row) = Value{0.};
        for (; offset != next_offset; ++offset) {
          auto [col, weight] =  _global_constraints[offset];
          containerEntry(x, row) += weight * containerEntry(y, col);
        }
      });
    }

    void applyLinearTransposedConstraints(auto& x) const {
      const auto y = x;
      forEachSparseRow([&](std::size_t i) {
        auto [row, _, offset] = _global_offset[i];
        containerEntry(x, row) = Value{0.};
      });
      forEachSparseRow([&](std::size_t i) {
        auto [row, _, offset] = _global_offset[i];
        auto next_offset = std::get<2>(_global_offset[i+1]);
        for (; offset != next_offset; ++offset) {
          auto [col, weight] =  _global_constraints[offset];
          containerEntry(x, col) += weight * containerEntry(y, row);
        }
      });
    }

    std::size_t size() const {
      return _map.size();
    }

    void clear() {
      _map.clear();
      _global_offset.clear();
      _global_constraints.clear();
    }

    void globalCompress() {
#if __cpp_lib_execution
      using std::execution::par;
      using std::sort;
#elif __has_include(<oneapi/dpl/execution>)
      using oneapi::dpl::sort;
      using oneapi::dpl::execution::par;
#endif

      _global_offset.clear();
      _global_constraints.clear();

      // sort values in map w.r.t container indices
      std::vector<std::pair<ContainerIndex, std::pair<Value,std::vector<std::pair<ContainerIndex, double>>>>> sorted_map;
      for(const auto& [row, constraints] : _map)
        sorted_map.push_back({row, constraints});
      // sort multi-indices with a lexicographic compare
      sort(par, begin(sorted_map), end(sorted_map),
        [](const auto& l, const auto& r){
          // (move to dynamic multi-index to make comparison easier)
          return Dune::PDELab::MultiIndex(l.first) < Dune::PDELab::MultiIndex(r.first);
        });

      // move sorted map to a compressed by row representation
      auto offset = 0;
      for(const auto& [row, constraints] : sorted_map) {
        const auto& [b_i, C_row] = constraints;
        _global_offset.push_back({row, b_i, offset});
        _global_constraints.insert(end(_global_constraints), begin(C_row), end(C_row));
        offset += C_row.size();
      }
      // last row only contains the final offset
      _global_offset.push_back({ContainerIndex{}, Value{}, offset});
    }

    void forEachSparseRow(auto&& at_row) const {
      assert(_global_offset.size() == _map.size() + 1 && "[Invalid State]: Constrains are not compressed");
#ifdef HAVE_TBB
      tbb::parallel_for(tbb::blocked_range<std::size_t>{0, _global_offset.size()-1}, [&](auto range){
        for (std::size_t i = std::begin(range); i != std::end(range); ++i)
          at_row(i);
      });
#else
      for(std::size_t i = 0; i != _global_offset.size()-1; ++i)
        at_row(i);
#endif
    }


    template<class T>
    class Hash;

    template<class... T>
    class Hash<TypeTree::HybridTreePath<T...>> {

    public:

      std::size_t operator()(const TypeTree::HybridTreePath<T...>& mi) const {
        std::size_t hash = 9999;
        using Tuple = typename TypeTree::HybridTreePath<T...>::Data;
        Hybrid::forEach(mi.enumerate(), [&](auto i) {
          using Index = std::tuple_element_t<i, Tuple>;
          if constexpr (not Dune::IsIntegralConstant<Index>::value) {
            hash_combine(hash, mi[i]);
          }
        });
        return hash;
      }
    };

    template<class... T>
    struct Hash<const TypeTree::HybridTreePath<T...>> : public Hash<TypeTree::HybridTreePath<T...>> {};

    // slow global sparse matrix of constraints
    using Row = std::pair<Value,std::vector<std::pair<ContainerIndex, double>>>;
#ifdef HAVE_TBB
    tbb::concurrent_unordered_map<ContainerIndex, Row, Hash<ContainerIndex>> _map;
#else
    std::unordered_map<ContainerIndex, Row, Hash<ContainerIndex>> _map;
#endif

    // fast global sparse matrix of constraints
    std::vector<std::tuple<ContainerIndex, Value, std::size_t>> _global_offset; // sparse constraint offsets
    std::vector<std::pair<ContainerIndex, double>> _global_constraints; // sparse constraints
  };

} // namespace Impl



  template<class Value, Concept::MultiIndex ContainerIndex, Dune::Concept::GridView EntitySet>
  class AffineConstraintsContainer
    : public TypeTree::LeafNode
    , private Impl::AffineConstraintsContainerBase<Value, ContainerIndex>
  {
    using Base = Impl::AffineConstraintsContainerBase<Value, ContainerIndex>;

    class LocalViewData {
    public:
      using size_type = std::size_t;

      LocalViewData(std::shared_ptr<const AffineConstraintsContainer> container)
        : _container{std::move(container)}
      {}

      void bind(const Dune::Concept::Entity auto& entity) noexcept {
        assert(_loffset.empty());
        assert(_lconstrained.empty());
        assert(_ltranslation.empty());
        assert(_llinear.empty());
        if (not _container->_mapper.gridView().indexSet().contains(entity)) return;
        const auto& casted_entity = entityCast(_container->_mapper.gridView(), entity);
        auto entity_index = _container->_mapper.index(casted_entity);
        auto [lebegin, leoffset] = _container->_local_entity_offset[entity_index];
        // local offsets are compressed on the global container
        auto _loffset_compressed = std::span(_container->_local_offset).subspan(lebegin, leoffset);
        if (_loffset_compressed.empty()) return;
        // we need to reconstruct the offsets locally
        auto max_dof = _loffset_compressed.back()[0];
        assert(max_dof < std::numeric_limits<size_type>::max()-2);
        _loffset.assign(max_dof+2, {0,0});
        _lconstrained.assign(max_dof+2, false);
        auto it = begin(_loffset_compressed);
        auto [_, clo_begin] = (*it);

        size_type cto = 0;
        size_type clo = 0;
        for (std::size_t dof = 0; dof+1 != _loffset.size(); ++dof) {
          _loffset[dof] = {cto, clo};
          if (dof == (*it)[0]) {
            _lconstrained[dof] = true;
            ++cto; // here we take advantage that _local_offset has the same size as _local_translation!
            clo = (*it)[1] - clo_begin;
            ++it;
          }
        }
        assert(it == end(_loffset_compressed));
        // notice that the compressed local offset range does not contain the last offset of linear constraints
        // so here we take advantage that (i) offsets where writen consecutively so the next entity will
        // have our last offset, and (ii) the local offset writen one past the last entity
        clo = _container->_local_offset[lebegin+leoffset][1] - clo_begin;
        _loffset.back() = {cto, clo};

        _ltranslation = std::span(_container->_local_translation).subspan(lebegin, leoffset);
        _llinear = std::span(_container->_local_linear).subspan(clo_begin, _loffset.back()[1]);
      }

      void unbind() {
        _loffset.clear();
        _lconstrained.clear();
        _ltranslation = {};
        _llinear = {};
      }

    protected:

      std::span<const Value> _ltranslation;
      std::span<const std::pair<ContainerIndex,double>> _llinear;
      std::vector<std::array<size_type,2>> _loffset; // dof->(_ltranslation, _llinear)
      std::vector<bool> _lconstrained;
      std::shared_ptr<const AffineConstraintsContainer> _container;
    };

  public:
    AffineConstraintsContainer(const EntitySet& entity_set)
      : _mapper{entity_set, mcmgElementLayout()}
    {}

    using Base::clear;
    using Base::addAffineConstraint;
    using Base::addLinearConstraint;
    using Base::addTranslationConstraint;
    using Base::applyAffineConstraints;
    using Base::applyLinearConstraints;
    using Base::applyLinearTransposedConstraints;



    class LeafLocalView
      : public TypeTree::LeafNode
      , public LocalViewData
    {
      using LocalViewData::_lconstrained;
      using LocalViewData::_loffset;
      using LocalViewData::_ltranslation;
      using LocalViewData::_llinear;
    public:
      using size_type = std::size_t;

      LeafLocalView(std::shared_ptr<const AffineConstraintsContainer> container)
        : LocalViewData{std::move(container)}
      {}

      [[nodiscard]] bool isConstrained(size_type dof) const noexcept {
        return dof >= _lconstrained.size() ? false : _lconstrained[dof];
       }

      [[nodiscard]] auto linearCoefficients(size_type dof) const noexcept {
        assert(isConstrained(dof));
        return _llinear.subspan(_loffset[dof][1], _loffset[dof+1][1] - _loffset[dof][1]);
      }

      [[nodiscard]] Value translationValue(size_type dof) const noexcept {
        assert(isConstrained(dof));
        return _ltranslation[_loffset[dof][0]];
      }
    };

    template<Concept::Tree Tree>
    static auto makeLocalViewNode(const Tree& tree, std::shared_ptr<const AffineConstraintsContainer> container) {
      if constexpr (Concept::LeafTreeNode<Tree>) {
        return std::make_shared<LeafLocalView>(std::move(container));
      } else { // skeleton finite element case
        static_assert(Concept::VectorTreeNode<Tree> and Concept::LeafTreeNode<typename Tree::ChildType>);
        //! @todo Implement sub entity local view on affine constraints
        static_assert(AlwaysFalse<Tree>{}, "Constrained skeleton finite elements is not yet implemeted");
        // typename SubEntityLocalView::NodeStorage storage(tree.tree());
        // for (std::size_t i = 0; i != storage.size(); ++i)
        //   storage[i] = std::make_shared<LeafLocalView>(container);
        // return std::make_shared<SubEntityLocalView>(std::move(storage));
      }
    }

    void globalCompress(Concept::Basis auto space) {
      Base::globalCompress();
      // initialize variables for local compression
      _local_entity_offset.assign(_mapper.size(), {});
    }

    // build the data structures for local constraints (not thread safe)
    template<Concept::LocalBasisTree LocalBasisTree>
    void localCompress(const LocalBasisTree& lbasis_tree) {
      if (lbasis_tree.size() == 0) return;
      std::size_t entity_index = _mapper.index(lbasis_tree.element());
      std::size_t _local_offset_begin = _local_offset.size();
      std::size_t node_offset = 0;
      // lbasis_tree may have children nodes if it refers to skeleton finite elements
      forEachLeafNode(lbasis_tree, [&](const auto& lbasis_leaf) {
        for (std::size_t dof = 0; dof != lbasis_leaf.size(); ++dof) {
          auto gdof = lbasis_leaf.index(dof);
          auto it = _map.find(gdof);
          if (it != _map.end()) {
            auto& [val, vec] = it->second;
            _local_offset.push_back({node_offset+dof,_local_linear.size()});
            _local_translation.emplace_back(val);
            _local_linear.insert(end(_local_linear), begin(vec), end(vec));
            // check if linear constraints are within the local function spaces!
            for (auto [cci, _] : vec) {
              bool cci_is_local = false;
              for (std::size_t jdof = 0; jdof != lbasis_leaf.size(); ++jdof)
                cci_is_local |= (cci == lbasis_leaf.index(jdof));
              // if (lbasis_tree.memoryRegion()) // TODO!
              if (not cci_is_local)
                DUNE_THROW(NotImplemented, "Non local constraints do not work with concurrent entity sets yet!");
            }
          }
        }
        node_offset += lbasis_leaf.size();
      });

      _local_entity_offset[entity_index] = {_local_offset_begin, _local_offset.size() - _local_offset_begin};

      // write one past the last offset to allow to recover offset of linear coefficients vector
      if (entity_index + 1 == _local_entity_offset.size())
        _local_offset.push_back({std::numeric_limits<std::size_t>::max(),_local_linear.size()});

    }

    // print all local offests in the terminal
    void debug() {
      std::cout << "\nLEO " << _local_entity_offset.size() << std::endl;
      for (std::size_t i = 0; i != _local_entity_offset.size(); ++i)
        std::cout << i << ": " << _local_entity_offset[i][0] << " " << _local_entity_offset[i][1] << std::endl;

      std::cout << "\nLO " << _local_offset.size() << std::endl;
      for (std::size_t i = 0; i != _local_offset.size(); ++i)
        std::cout << i << ": " << _local_offset[i][0] << " " << _local_offset[i][1] << std::endl;

      std::cout << "\nLT " << _local_translation.size() << std::endl;
      for (std::size_t i = 0; i != _local_translation.size(); ++i)
        std::cout << i << ": " << _local_translation[i] << std::endl;

      std::cout << "\nLL " << _local_linear.size() << std::endl;
      for (std::size_t i = 0; i != _local_linear.size(); ++i)
        std::cout << i << ": (" << _local_linear[i].first << " " << _local_linear[i].second << ")" << std::endl;
    }

  private:

    using Base::_global_offset;
    using Base::_map;

    MultipleCodimMultipleGeomTypeMapper<EntitySet> _mapper;
    // fast local matrix of constraints
    std::vector<std::array<std::size_t,2>> _local_entity_offset;
    std::vector<std::array<std::size_t,2>> _local_offset;
    std::vector<Value> _local_translation;
    std::vector<std::pair<ContainerIndex, double>> _local_linear;
  };

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_BASIS_CONSTRAINTS_CONTAINER_AFFINE_HH
