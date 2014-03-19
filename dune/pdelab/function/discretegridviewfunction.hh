// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_FUNCTION_DISCRETEGRIDVIEWFUNCTION_HH
#define DUNE_PDELAB_FUNCTION_DISCRETEGRIDVIEWFUNCTION_HH

#include <cstdlib>
#include <vector>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/static_assert.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/common/jacobiantocurl.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>

namespace Dune {
namespace PDELab {

template<typename DGVF>
class DiscreteLocalGridViewFunction;

template<typename GFS, typename F>
class DiscreteGridViewFunction
  : public Functions::GridViewFunction<typename GFS::Traits::GridView,
                                       typename GFS::Traits::FiniteElement::Traits::LocalBasisType::Traits::RangeType
                                       >
{

  typedef typename BackendVectorSelector<GFS,F>::Type V;
  typedef Functions::GridViewFunction<
    typename GFS::Traits::GridView,
    typename GFS::Traits::FiniteElement::Traits::LocalBasisType::Traits::RangeType
    > Base;

public:

  typedef typename Base::Element Element;
  typedef GFS GridFunctionSpace;
  typedef typename BackendVectorSelector<GFS,F>::Type Vector;

  typedef DiscreteLocalGridViewFunction<
    DiscreteGridViewFunction
    > LocalFunction;

  typedef typename Base::Domain Domain;
  typedef typename Base::Range Range;

  virtual typename Base::LocalFunctionBasePointer localFunction() const  DUNE_FINAL
  {
    return make_shared<LocalFunction>(_pgfs, _v);
  }

  virtual typename Base::DerivativeBasePointer derivative() const DUNE_FINAL
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  virtual void evaluate(const Domain& domain, Range& r) const DUNE_FINAL
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  DiscreteGridViewFunction(const GFS& gfs, const V& v)
    : Base(gfs.gridView())
    , _pgfs(stackobject_to_shared_ptr(gfs))
    , _v(stackobject_to_shared_ptr(v))
  {}

  DiscreteGridViewFunction(std::shared_ptr<const GFS> pgfs, std::shared_ptr<const V> v)
    : Base(pgfs->gridView())
    , _pgfs(pgfs)
    , _v(v)
  {}

  const GridFunctionSpace& gridFunctionSpace() const
  {
    return *_pgfs;
  }

  const V& dofs() const
  {
    return *_v;
  }

private:

  friend LocalFunction;

  const shared_ptr<const GFS> _pgfs;
  const shared_ptr<const V> _v;

};

template<typename DGVF>
class DiscreteLocalGridViewFunction
  : public DGVF::Base::LocalFunction
{

  typedef DGVF DiscreteGridViewFunction;
  typedef typename DGVF::Vector Vector;
  typedef typename DGVF::GridFunctionSpace GFS;
  typedef typename DGVF::Base::LocalFunction EBase;

public:

  typedef typename EBase::LocalContext Element;

  typedef typename EBase::Domain Domain;
  typedef typename EBase::Range Range;

  DiscreteLocalGridViewFunction(const shared_ptr<const GFS> gfs, const shared_ptr<const Vector> v)
    : _pgfs(gfs)
    , _v(v)
    , _lfs(*_pgfs)
    , _lfs_cache(_lfs)
    , _x_view(*_v)
    , _xl(_pgfs->maxLocalSize())
    , _yb(_pgfs->maxLocalSize())
    , _element(nullptr)
  {}

  virtual typename EBase::DerivativeBasePointer derivative() const DUNE_FINAL
  {
    DUNE_THROW(NotImplemented,"bla");
  }

  virtual void bind(const Element& element) DUNE_FINAL
  {
    _element = &element;
    _lfs.bind(element);
    _lfs_cache.update();
    _x_view.bind(_lfs_cache);
    _x_view.read(_xl);
    _x_view.unbind();
  }

  virtual void unbind() DUNE_FINAL
  {
    _element = nullptr;
  }

  virtual void evaluate(const Domain& coord, Range& r) const DUNE_FINAL
  {
    _lfs.finiteElement().localBasis().evaluateFunction(coord,_yb);
    r = 0;
    for (unsigned int i = 0; i < _yb.size(); ++i)
      {
        r.axpy(_xl[i],_yb[i]);
      }
  }

  virtual const Element& localContext() const  DUNE_FINAL
  {
#ifndef NDEBUG
    if (!_element)
      DUNE_THROW(InvalidStateException,"bla");
#endif
    return *_element;
  }

private:

  typedef LocalFunctionSpace<GFS> LFS;
  typedef LFSIndexCache<LFS> LFSCache;
  typedef typename Vector::template ConstLocalView<LFSCache> XView;

  const shared_ptr<const GFS> _pgfs;
  const shared_ptr<const Vector> _v;
  LFS _lfs;
  LFSCache _lfs_cache;
  XView _x_view;
  mutable std::vector<typename Vector::ElementType> _xl;
  mutable std::vector<Range> _yb;
  const Element* _element;

};

} // end of namespace Dune::Functions
} // end of namespace Dune

#endif // DUNE_PDELAB_FUNCTION_DISCRETEGRIDVIEWFUNCTION_HH
