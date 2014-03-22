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

template<class DT, class RT, int N>
struct EvaluateDerivativeTraits
{
  typedef typename Functions::DerivativeTraits<
    DT,
    typename EvaluateDerivativeTraits<DT,RT,N-1>::RangeType
    >::DerivativeRange DerivativeRange;
};

template<class DT, class RT>
struct EvaluateDerivativeTraits<DT,RT,1>
{
  typedef typename Functions::DerivativeTraits<
    DT, RT>::DerivativeRange DerivativeRange;
};

template<class DT, class RT>
struct EvaluateDerivativeTraits<DT,RT,0>
{
  typedef RT DerivativeRange;
};

template<typename LocalFunction>
class DiscreteGridViewFunctionBase;

template<typename Traits>
class DiscreteLocalGridViewFunctionBase;

template<typename GFS, typename F>
class DiscreteLocalGridViewFunction;

template<typename GFS, typename V>
class DiscreteLocalGridViewFunctionJacobian;

template<typename GFS, typename F>
struct DiscreteGridViewFunctionTraits
{
  typedef GFS GridFunctionSpace;
  typedef typename GFS::Traits::GridView GridView;
  typedef Functions::GridViewEntitySet<GridView, 0> EntitySet;

  typedef typename EntitySet::GlobalCoordinate Domain;
  //! range type of the function
  typedef typename GFS::Traits::FiniteElement::Traits::LocalBasisType::Traits::RangeType Range;
  typedef Functions::GridViewFunction<GridView, Range> FunctionInterface;

  typedef F VectorFieldType;
  typedef typename BackendVectorSelector<GridFunctionSpace,VectorFieldType>::Type Vector;

  typedef Range LFERange;
  enum { RangeExists = 1 };
};

template<typename GFS, typename F, int N>
struct DiscreteGridViewFunctionDerivativeTraits :
    public DiscreteGridViewFunctionTraits<GFS,F>
{
  //! range type of the N'th derivative
  typedef typename EvaluateDerivativeTraits<
    typename DiscreteGridViewFunctionTraits<GFS,F>::Domain,
    typename DiscreteGridViewFunctionTraits<GFS,F>::Range,
    N
    >::DerivativeRange DerivativeRange;
  typedef Functions::GridViewFunction<
    typename DiscreteGridViewFunctionTraits<GFS,F>::GridView,
    DerivativeRange> FunctionInterface;
  enum { RangeExists =
         ((GFS::Traits::FiniteElement::Traits::LocalBasisType::Traits::diffOrder > N)
           &&
           std::is_same<DerivativeRange,Functions::InvalidRange>::value
           ) ? 1 : 0 };
};

template<typename GFS, typename F>
struct DiscreteGridViewFunctionJacobianTraits :
    public DiscreteGridViewFunctionDerivativeTraits<GFS,F,1>
{
  typedef typename GFS::Traits::FiniteElement::Traits::LocalBasisType::Traits::JacobianType Jacobian;
  typedef Jacobian LFERange;
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
    return make_shared<LocalFunction>(pgfs_, v_);
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
    , pgfs_(pgfs)
    , v_(v)
  {}

  const GridFunctionSpace& gridFunctionSpace() const
  {
    return *pgfs_;
  }

  const Vector& dofs() const
  {
    return *v_;
  }

private:

  const shared_ptr<const GridFunctionSpace> pgfs_;
  const shared_ptr<const Vector> v_;

};

template<typename T>
class DiscreteLocalGridViewFunctionBase
  : public T::FunctionInterface::LocalFunction
{

  typedef typename T::FunctionInterface::LocalFunction Base;

public:
  typedef T Traits;

  typedef typename Base::LocalContext Element;
  typedef typename Base::Domain Domain;
  typedef typename Base::Range Range;

  typedef typename Traits::Vector Vector;
  typedef typename Traits::GridFunctionSpace GridFunctionSpace;

  DiscreteLocalGridViewFunctionBase(const shared_ptr<const GridFunctionSpace> gfs, const shared_ptr<const Vector> v)
    : pgfs_(gfs)
    , v_(v)
    , lfs_(*pgfs_)
    , lfs_cache_(lfs_)
    , x_view_(*v_)
    , xl_(pgfs_->maxLocalSize())
    , yb_(pgfs_->maxLocalSize())
    , element_(nullptr)
  {}

  virtual void bind(const Element& element) DUNE_FINAL
  {
    element_ = &element;
    lfs_.bind(element);
    lfs_cache_.update();
    x_view_.bind(lfs_cache_);
    x_view_.read(xl_);
    x_view_.unbind();
  }

  virtual void unbind() DUNE_FINAL
  {
    element_ = nullptr;
  }

  virtual const Element& localContext() const  DUNE_FINAL
  {
#ifndef NDEBUG
    if (!element_)
      DUNE_THROW(InvalidStateException,"bla");
#endif
    return *element_;
  }

protected:

  typedef LocalFunctionSpace<GridFunctionSpace> LFS;
  typedef LFSIndexCache<LFS> LFSCache;
  typedef typename Vector::template ConstLocalView<LFSCache> XView;

  const shared_ptr<const GridFunctionSpace> pgfs_;
  const shared_ptr<const Vector> v_;
  LFS lfs_;
  LFSCache lfs_cache_;
  XView x_view_;
  mutable std::vector<typename Vector::ElementType> xl_;
  typedef typename Traits::LFERange LFERange;
  mutable std::vector<LFERange> yb_;
  const Element* element_;

};

template<typename GFS, typename V>
class DiscreteLocalGridViewFunction
  : public DiscreteLocalGridViewFunctionBase< DiscreteGridViewFunctionTraits<GFS,V> >
{

  typedef DiscreteLocalGridViewFunctionBase< DiscreteGridViewFunctionTraits<GFS,V> > Base;

  using Base::lfs_;
  using Base::yb_;
  using Base::xl_;
  using Base::pgfs_;
  using Base::v_;
  using Base::element_;

public:

  typedef typename Base::Domain Domain;
  typedef typename Base::Range Range;

  typedef typename Base::Vector Vector;
  typedef typename Base::GridFunctionSpace GridFunctionSpace;

  typedef DiscreteLocalGridViewFunctionJacobian<GFS,V> Derivative;

  DiscreteLocalGridViewFunction(const shared_ptr<const GridFunctionSpace> gfs, const shared_ptr<const Vector> v)
    : Base(gfs,v)
  {}

  virtual typename Base::DerivativeBasePointer derivative() const DUNE_FINAL
  {
    shared_ptr<Derivative> diff = make_shared<Derivative>(pgfs_,v_);
    // TODO: do we really want this?
    if (element_) diff->bind(*element_);
    return diff;
  }

  virtual void evaluate(const Domain& coord, Range& r) const DUNE_FINAL
  {
    lfs_.finiteElement().localBasis().evaluateFunction(coord,yb_);
    r = 0;
    for (unsigned int i = 0; i < yb_.size(); ++i)
      {
        r.axpy(xl_[i],yb_[i]);
      }
  }

};

template<typename GFS, typename V>
class DiscreteLocalGridViewFunctionJacobian
  : public DiscreteLocalGridViewFunctionBase< DiscreteGridViewFunctionJacobianTraits<GFS,V> >
{

  typedef DiscreteLocalGridViewFunctionBase< DiscreteGridViewFunctionJacobianTraits<GFS,V> > Base;

  using Base::lfs_;
  using Base::yb_;
  using Base::xl_;
  using Base::element_;

public:

  typedef typename Base::Domain Domain;
  typedef typename Base::Range Range;

  typedef typename Base::Vector Vector;
  typedef typename Base::GridFunctionSpace GridFunctionSpace;

  DiscreteLocalGridViewFunctionJacobian(const shared_ptr<const GridFunctionSpace> gfs, const shared_ptr<const Vector> v)
    : Base(gfs,v)
  {}

  virtual typename Base::DerivativeBasePointer derivative() const DUNE_FINAL
  {
    DUNE_THROW(NotImplemented,"sorry, no further derivatives");
  }

  virtual void evaluate(const Domain& coord, Range& r) const DUNE_FINAL
  {
    // get Jacobian of geometry
    const typename Base::Element::Geometry::JacobianInverseTransposed
      JgeoIT = element_->geometry().jacobianInverseTransposed(coord);

    // get local Jacobians/gradients of the shape functions
    lfs_.finiteElement().localBasis().evaluateJacobian(coord,yb_);

    typename Base::Range gradphi;
    r = 0;
    for(unsigned int i = 0; i < yb_.size(); ++i) {
      // compute global gradient of shape function i
      gradphi = 0;
      // TODO: in general this must be a matrix matrix product
      JgeoIT.umv(yb_[i][0], gradphi);

      // sum up global gradients, weighting them with the appropriate coeff
      r.axpy(xl_[i], gradphi);
    }
  }

};

} // end of namespace Dune::PDELab
} // end of namespace Dune

#endif // DUNE_PDELAB_FUNCTION_DISCRETEGRIDVIEWFUNCTION_HH
