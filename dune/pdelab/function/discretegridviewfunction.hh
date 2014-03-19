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

template<typename GFS, typename V>
class DiscreteGridViewFunction
  : public Functions::GridViewFunction<typename GFS::Traits::GridView,
                            typename GFS::Traits::FiniteElement::Traits::LocalBasisType::Traits::RangeType
                            >
{

  typedef Functions::GridViewFunction<
    typename GFS::Traits::GridView,
    typename GFS::Traits::FiniteElement::Traits::LocalBasisType::Traits::RangeType
    > Base;

public:

  typedef typename Base::Element Element;

  class LocalFunction
    : public Base::ElementFunction
  {

    typedef typename Base::LocaFunction EBase;

  public:

    typedef typename EBase::LocalContext Element;

    typedef typename EBase::Domain Domain;
    typedef typename EBase::Range Range;

    LocalFunction(shared_ptr<const DiscreteGridViewFunction> dgvf)
      : _dgvf(dgvf)
      , _lfs(dgvf->gridFunctionSpace())
      , _lfs_cache(_lfs)
      , _x_view(dgvf->dofs())
      , _xl(dgvf->gridFunctionSpace().maxLocalSize())
      , _yb(dgvf->gridFunctionSpace().maxLocalSize())
      , _element(nullptr)
    {}

    virtual typename EBase::Derivative* derivative() const DUNE_FINAL
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

    typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
    typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
    typedef typename V::template ConstLocalView<LFSCache> XView;

    shared_ptr<DiscreteGridViewFunction const> _dgvf;
    LFS _lfs;
    LFSCache _lfs_cache;
    XView _x_view;
    mutable std::vector<typename V::ElementType> _xl;
    mutable std::vector<Range> _yb;
    const Element* _element;

  };

  typedef typename Base::Domain Domain;
  typedef typename Base::Range Range;

  virtual typename Base::LocalFunctionBasePointer localFunction() const  DUNE_FINAL
  {
    return make_shared<LocalFunction>(std::static_pointer_cast<const DiscreteGridViewFunction>(this->shared_from_this()));
  }

  virtual typename Base::Derivative* derivative() const DUNE_FINAL
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  virtual void evaluate(const Domain& domain, Range& r) const DUNE_FINAL
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  DiscreteGridViewFunction(const GFS& gfs, const V& v)
    : Base(gfs.gridView())
    , _gfs(gfs)
    , _v(v)
  {}

  const GFS& gridFunctionSpace() const
  {
    return _gfs;
  }

  const V& dofs() const
  {
    return _v;
  }

private:

  const GFS& _gfs;
  const V& _v;

};

} // end of namespace Dune::Functions
} // end of namespace Dune

#endif // DUNE_PDELAB_FUNCTION_DISCRETEGRIDVIEWFUNCTION_HH
