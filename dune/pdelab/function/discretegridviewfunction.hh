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

template<typename LocalFunction>
class DiscreteGridViewFunctionBase;

template<typename GFS, typename F>
class DiscreteLocalGridViewFunctionBase;

template<typename GFS, typename F>
class DiscreteLocalGridViewFunction;

template<typename GFS, typename F>
struct DiscreteGridViewFunctionTraits
{
  typedef GFS GridFunctionSpace;
  typedef typename GFS::Traits::GridView GridView;
  typedef Functions::GridViewEntitySet<GridView, 0> EntitySet;

  typedef typename EntitySet::GlobalCoordinate Domain;
  typedef typename GFS::Traits::FiniteElement::Traits::LocalBasisType::Traits::RangeType Range;
  typedef Functions::GridViewFunction<GridView, Range> FunctionInterface;

  typedef F VectorFieldType;
  typedef typename BackendVectorSelector<GridFunctionSpace,VectorFieldType>::Type Vector;
};

template<typename GFS, typename F>
class DiscreteGridViewFunction
  : public DiscreteGridViewFunctionBase< DiscreteLocalGridViewFunction<GFS,F> >
{
  typedef DiscreteGridViewFunctionBase< DiscreteLocalGridViewFunction<GFS,F> > Base;
public:

  typedef typename Base::GridFunctionSpace GridFunctionSpace;
  typedef typename Base::Vector Vector;

  DiscreteGridViewFunction(const GridFunctionSpace& gfs, const Vector& v)
    : Base(stackobject_to_shared_ptr(gfs),stackobject_to_shared_ptr(v))
  {}

  DiscreteGridViewFunction(std::shared_ptr<const GridFunctionSpace> pgfs, std::shared_ptr<const Vector> v)
    : Base(pgfs,v)
  {}
};

template<typename LocalFnkt>
class DiscreteGridViewFunctionBase
  : public LocalFnkt::Traits::FunctionInterface
{

  typedef typename LocalFnkt::Traits Traits;
  typedef typename Traits::FunctionInterface Base;

public:

  typedef typename Base::Element Element;
  typedef typename Base::Domain Domain;
  typedef typename Base::Range Range;

  typedef LocalFnkt LocalFunction;

  typedef typename Traits::GridFunctionSpace GridFunctionSpace;
  typedef typename Traits::Vector Vector;

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

  DiscreteGridViewFunctionBase(std::shared_ptr<const GridFunctionSpace> pgfs, std::shared_ptr<const Vector> v)
    : Base(pgfs->gridView())
    , _pgfs(pgfs)
    , _v(v)
  {}

  const GridFunctionSpace& gridFunctionSpace() const
  {
    return *_pgfs;
  }

  const Vector& dofs() const
  {
    return *_v;
  }

private:

  const shared_ptr<const GridFunctionSpace> _pgfs;
  const shared_ptr<const Vector> _v;

};

template<typename GFS, typename V>
class DiscreteLocalGridViewFunctionBase
  : public DiscreteGridViewFunctionTraits<GFS,V>::FunctionInterface::LocalFunction
{

public:
  typedef DiscreteGridViewFunctionTraits<GFS,V> Traits;

private:
  typedef typename Traits::FunctionInterface::LocalFunction Base;

public:

  typedef typename Base::LocalContext Element;
  typedef typename Base::Domain Domain;
  typedef typename Base::Range Range;

  typedef typename Traits::Vector Vector;
  typedef typename Traits::GridFunctionSpace GridFunctionSpace;

  DiscreteLocalGridViewFunctionBase(const shared_ptr<const GridFunctionSpace> gfs, const shared_ptr<const Vector> v)
    : _pgfs(gfs)
    , _v(v)
    , _lfs(*_pgfs)
    , _lfs_cache(_lfs)
    , _x_view(*_v)
    , _xl(_pgfs->maxLocalSize())
    , _yb(_pgfs->maxLocalSize())
    , _element(nullptr)
  {}

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

  virtual const Element& localContext() const  DUNE_FINAL
  {
#ifndef NDEBUG
    if (!_element)
      DUNE_THROW(InvalidStateException,"bla");
#endif
    return *_element;
  }

protected:

  typedef LocalFunctionSpace<GridFunctionSpace> LFS;
  typedef LFSIndexCache<LFS> LFSCache;
  typedef typename Vector::template ConstLocalView<LFSCache> XView;

  const shared_ptr<const GridFunctionSpace> _pgfs;
  const shared_ptr<const Vector> _v;
  LFS _lfs;
  LFSCache _lfs_cache;
  XView _x_view;
  mutable std::vector<typename Vector::ElementType> _xl;
  mutable std::vector<Range> _yb;
  const Element* _element;

};

template<typename GFS, typename V>
class DiscreteLocalGridViewFunction
  : public DiscreteLocalGridViewFunctionBase<GFS,V>
{

  typedef DiscreteLocalGridViewFunctionBase<GFS,V> Base;

  using Base::_lfs;
  using Base::_yb;
  using Base::_xl;

public:

  typedef typename Base::Domain Domain;
  typedef typename Base::Range Range;

  typedef typename Base::Vector Vector;
  typedef typename Base::GridFunctionSpace GridFunctionSpace;

  DiscreteLocalGridViewFunction(const shared_ptr<const GridFunctionSpace> gfs, const shared_ptr<const Vector> v)
    : Base(gfs,v)
  {}

  virtual typename Base::DerivativeBasePointer derivative() const DUNE_FINAL
  {
    DUNE_THROW(NotImplemented,"bla");
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

};

} // end of namespace Dune::Functions
} // end of namespace Dune

#endif // DUNE_PDELAB_FUNCTION_DISCRETEGRIDVIEWFUNCTION_HH
