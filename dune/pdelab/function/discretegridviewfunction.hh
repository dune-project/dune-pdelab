// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_FUNCTION_DISCRETEGRIDVIEWFUNCTION_HH
#define DUNE_PDELAB_FUNCTION_DISCRETEGRIDVIEWFUNCTION_HH

#include <cstdlib>
#include <vector>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/common/jacobiantocurl.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>

#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

namespace Dune {
namespace PDELab {

  namespace Imp {
/**
   \todo use GridViewTraits
   \todo fix relation to LocalDerivativeTraits
 */

    template<class Signature, int NthDerivative, template<class> class DerivativeTraits>
    struct EvaluateDerivativeTraits;

    template<class DT, class RT, int NthDerivative, template<class> class DerivativeTraits>
    struct EvaluateDerivativeTraits<RT(DT), NthDerivative, DerivativeTraits>
    {
      using InputRange =
        typename EvaluateDerivativeTraits<DT(RT),NthDerivative-1,DerivativeTraits>::Range;
      using Range = typename DerivativeTraits<InputRange(DT)>::Range;
    };

    template<class DT, class RT, template<class> class DerivativeTraits>
    struct EvaluateDerivativeTraits<RT(DT), 1, DerivativeTraits>
    {
      using Range = typename DerivativeTraits<RT(DT)>::Range;
    };

    template<class DT, class RT, template<class> class DerivativeTraits>
    struct EvaluateDerivativeTraits<RT(DT), 0, DerivativeTraits>
    {
      using Range = RT;
    };

    template<typename GFS, typename F, int NthDerivative,
             template<class> class DerivativeTraits=Functions::DefaultDerivativeTraits,
             std::size_t bufferSize = 16>
    struct DiscreteGridViewFunctionTraits
    {
      //! the GridFunctionSpace we are operating on
      using GridFunctionSpace = GFS;
      //! the underlying GridView
      using GridView = typename GFS::Traits::GridView;
      //! the type of the corresponding codim-0 EntitySet
      using EntitySet = Functions::GridViewEntitySet<GridView, 0>;
      //! the type of the codim-0 entity
      using Element = typename EntitySet::Element;

      //! domain type (aka. coordinate type of the world dimensions)
      using Domain = typename EntitySet::GlobalCoordinate;
      //! range type of the underlying local function
      using LocalFiniteElementRange = typename GFS::Traits::FiniteElement::Traits::LocalBasisType::Traits::RangeType;
      //! range type of the initial function
      using BasicRange = LocalFiniteElementRange; // LocalFiniteElementRange
      //! data type of the vector container
      using VectorFieldType = F;
      //! type of the vector container
      using Vector = typename BackendVectorSelector<GridFunctionSpace,VectorFieldType>::Type;

      //! range type of the N'th derivative
      using Range = typename EvaluateDerivativeTraits<BasicRange(Domain), NthDerivative, DerivativeTraits>::Range;

      //! Signature of the function
      using Signature = Range(Domain);

      //! the function interface we are providing
      using FunctionInterface = Functions::GridViewFunction<GridView, Range>;

      //! Signature of the derivative
      using LocalSignature = Range(typename EntitySet::LocalCoordinate);

      //! DerivativeTraits for the local functions
      template<class R>
        using LocalDerivativeTraits = typename Dune::Functions::LocalDerivativeTraits<EntitySet, DerivativeTraits>::template Traits<R>;

      //! LocalFunctionTraits associated with this type
      using LocalFunctionTraits = typename Dune::Functions::Imp::LocalFunctionTraits<LocalSignature, Element, LocalDerivativeTraits, bufferSize>;

      //! Interface type of the local function
      using LocalFunctionInterface = Dune::Functions::LocalFunction<LocalSignature, Element, LocalDerivativeTraits, bufferSize>;

      //! export how often the initial function can be differentiated
      enum { maxDiffOrder = 1 }; // GFS::Traits::FiniteElement::Traits::LocalBasisType::Traits::diffOrder };
      //! export which derivative we are currently evaluating
      enum { diffOrder = NthDerivative };
      //! export whether the range for the N'th derivative exists
      enum { RangeExists =
             ( (std::size_t(diffOrder) > std::size_t(maxDiffOrder))
               ||
               std::is_same<Range,Functions::InvalidRange>::value
               ) ? 0 : 1 };
    };

  } // end namespace Imp

  /* Forward Declarations */

template<typename LocalFunction>
class DiscreteGridViewFunctionBase;

template<typename Traits, typename LFERange>
class DiscreteLocalGridViewFunctionBase;

template<typename GFS, typename F>
class DiscreteLocalGridViewFunction;

template<typename GFS, typename V, int N,
         bool DerivativeExists = Imp::DiscreteGridViewFunctionTraits<GFS,V,N+1>::RangeExists>
class DiscreteLocalGridViewFunctionDerivative;

template<typename GFS, typename F>
class DiscreteGridViewFunction
  : public DiscreteGridViewFunctionBase< DiscreteLocalGridViewFunction<GFS,F> >
{
  using Base = DiscreteGridViewFunctionBase< DiscreteLocalGridViewFunction<GFS,F> >;
public:

  using GridFunctionSpace = typename Base::GridFunctionSpace;
  using Vector = typename Base::Vector;

  DiscreteGridViewFunction(const GridFunctionSpace& gfs, const Vector& v)
    : Base(stackobject_to_shared_ptr(gfs),stackobject_to_shared_ptr(v))
  {}

  DiscreteGridViewFunction(std::shared_ptr<const GridFunctionSpace> pgfs, std::shared_ptr<const Vector> v)
    : Base(pgfs,v)
  {}
};

template<typename LocalFnkt>
class DiscreteGridViewFunctionBase
{

public:

  using Traits = typename LocalFnkt::GlobalTraits;

  using Domain = typename Traits::Domain;
  using Range = typename Traits::Range;
  using EntitySet = typename Traits::EntitySet;
  using Element = typename Traits::Element;

  using Vector = typename Traits::Vector;
  using GridFunctionSpace = typename Traits::GridFunctionSpace;

  using LocalFunction = LocalFnkt;

  using LocalDerivative = typename LocalFunction::LocalDerivative;
  using Derivative = DiscreteGridViewFunctionBase<LocalDerivative>;

  friend LocalFunction localFunction(const DiscreteGridViewFunctionBase & f)
  {
    return LocalFunction(f.pgfs_, f.v_);
  }

  friend Derivative derivative(const DiscreteGridViewFunctionBase & f)
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  Range operator()(const Domain& domain) const
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  DiscreteGridViewFunctionBase(std::shared_ptr<const GridFunctionSpace> pgfs, std::shared_ptr<const Vector> v)
    : pgfs_(pgfs)
    , v_(v)
    , entitySet_(pgfs->gridView())
  {}

  const GridFunctionSpace& gridFunctionSpace() const
  {
    return *pgfs_;
  }

  const Vector& dofs() const
  {
    return *v_;
  }

  /**
   * \brief Get associated EntitySet
   */
  const EntitySet& entitySet() const
  {
    return entitySet_;
  }

private:

  const shared_ptr<const GridFunctionSpace> pgfs_;
  const shared_ptr<const Vector> v_;
  EntitySet entitySet_;

};

template<typename T, typename LFERange>
class DiscreteLocalGridViewFunctionBase
{
public:
  using GlobalTraits = T;
  using Traits = typename GlobalTraits::LocalFunctionTraits;

  using GlobalFunction = int;
  using Domain = typename Traits::Domain;
  using Range = typename Traits::Range;
  using LocalContext = typename Traits::LocalContext;
  using Element = LocalContext;

  using Vector = typename GlobalTraits::Vector;
  using GridFunctionSpace = typename GlobalTraits::GridFunctionSpace;

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

  void bind(const Element& element)
  {
    element_ = &element;
    lfs_.bind(element);
    lfs_cache_.update();
    x_view_.bind(lfs_cache_);
    x_view_.read(xl_);
    x_view_.unbind();
  }

  void unbind()
  {
    element_ = nullptr;
  }

  const Element& localContext() const
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
  mutable std::vector<LFERange> yb_;
  const Element* element_;

};

template<typename GFS, typename V>
class DiscreteLocalGridViewFunction
  : public DiscreteLocalGridViewFunctionBase<
  Imp::DiscreteGridViewFunctionTraits<GFS,V,0>,
  typename Imp::DiscreteGridViewFunctionTraits<GFS,V,0>::Range
  >
{

  typedef DiscreteLocalGridViewFunctionBase<
    Imp::DiscreteGridViewFunctionTraits<GFS,V,0>,
    typename Imp::DiscreteGridViewFunctionTraits<GFS,V,0>::Range
    > Base;

  using Base::lfs_;
  using Base::yb_;
  using Base::xl_;
  using Base::pgfs_;
  using Base::v_;
  using Base::element_;

public:

  using GlobalFunction = DiscreteGridViewFunction<GFS,V>;
  using Element = typename Base::Element;
  using Domain = typename Base::Domain;
  using Range = typename Base::Range;

  using Vector = typename Base::Vector;
  using GridFunctionSpace = typename Base::GridFunctionSpace;

  #warning to work around buggy dune-localfunctions implementations we force the first derivative to be available
  using LocalDerivative = DiscreteLocalGridViewFunctionDerivative<GFS,V,1,true>;

  DiscreteLocalGridViewFunction(const shared_ptr<const GridFunctionSpace> gfs, const shared_ptr<const Vector> v)
    : Base(gfs,v)
  {}

  Range operator()(const Domain& coord) const
  {
    lfs_.finiteElement().localBasis().evaluateFunction(coord,yb_);
    Range r(0);
    for (unsigned int i = 0; i < yb_.size(); ++i)
    {
      r.axpy(xl_[i],yb_[i]);
    }
    return r;
  }

  friend LocalDerivative derivative(const DiscreteLocalGridViewFunction& f)
  {
    return LocalDerivative(f.pgfs_,f.v_);
  }

};

template<typename GFS, typename V, bool DerivativeExists>
class DiscreteLocalGridViewFunctionDerivative<GFS,V,1,DerivativeExists>
  : public DiscreteLocalGridViewFunctionBase<
  Imp::DiscreteGridViewFunctionTraits<GFS,V,1>,
  typename Imp::DiscreteGridViewFunctionTraits<GFS,V,1>::Range
  >
{

  typedef DiscreteLocalGridViewFunctionBase<
    Imp::DiscreteGridViewFunctionTraits<GFS,V,1>,
    typename Imp::DiscreteGridViewFunctionTraits<GFS,V,1>::Range
    > Base;

  using Base::lfs_;
  using Base::yb_;
  using Base::xl_;
  using Base::pgfs_;
  using Base::v_;
  using Base::element_;

public:

  using GlobalFunction = DiscreteGridViewFunctionBase<DiscreteLocalGridViewFunctionDerivative>;
  using Element = typename Base::Element;
  using Domain = typename Base::Domain;
  using Range = typename Base::Range;

  using Vector = typename Base::Vector;
  using GridFunctionSpace = typename Base::GridFunctionSpace;

  DiscreteLocalGridViewFunctionDerivative(const shared_ptr<const GridFunctionSpace> gfs, const shared_ptr<const Vector> v)
    : Base(gfs,v)
  {}

  using LocalDerivative = DiscreteLocalGridViewFunctionDerivative<GFS,V,2>;

  Range operator()(const Domain& coord) const
  {
    // get Jacobian of geometry
    const typename Base::Element::Geometry::JacobianInverseTransposed
      JgeoIT = element_->geometry().jacobianInverseTransposed(coord);

    // get local Jacobians/gradients of the shape functions
    lfs_.finiteElement().localBasis().evaluateJacobian(coord,yb_);

    Range gradphi;
    Range r(0);
    for(std::size_t i = 0; i < yb_.size(); ++i) {
      assert(gradphi.size() == yb_[i].size());
      for(std::size_t j = 0; j < gradphi.size(); ++j) {
        // compute global gradient of shape function i
        // graphi += {J^{-1}}^T * yb_i0
        JgeoIT.mv(yb_[i][j], gradphi[j]);

        // sum up global gradients, weighting them with the appropriate coeff
        // r \in R^{1,dim}
        // r_0 += xl_i * grad \phi
        r[j].axpy(xl_[i], gradphi[j]);
      }
    }

    return r;
  }

  friend LocalDerivative derivative(const DiscreteLocalGridViewFunctionDerivative & f)
  {
    return LocalDerivative(f.pgfs_,f.v_);
  }

};

template<typename GFS, typename V, int N, bool>
class DiscreteLocalGridViewFunctionDerivative
  : public DiscreteLocalGridViewFunctionBase<
  Imp::DiscreteGridViewFunctionTraits<GFS,V,N>,
  typename Imp::DiscreteGridViewFunctionTraits<GFS,V,N>::LocalFiniteElementRange
  >
{

  typedef DiscreteLocalGridViewFunctionBase<
    Imp::DiscreteGridViewFunctionTraits<GFS,V,N>,
    typename Imp::DiscreteGridViewFunctionTraits<GFS,V,N>::LocalFiniteElementRange
    > Base;

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

  using LocalDerivative = DiscreteLocalGridViewFunctionDerivative<GFS,V,N+1>;

  DiscreteLocalGridViewFunctionDerivative(const shared_ptr<const GridFunctionSpace> gfs, const shared_ptr<const Vector> v)
    : Base(gfs,v)
  {}

  Range operator()(const Domain& coord) const
  {
    // TODO: we currently require affine geometries.
    if (! element_->geometry().affine())
      DUNE_THROW(NotImplemented, "Due to missing higher derivatives in the Geometry interface, "
        "the computation of higher derivatives (>=2) works only for affine transformations.");
    // get Jacobian of geometry
    const typename Base::Element::Geometry::JacobianInverseTransposed
      JgeoIT = element_->geometry().jacobianInverseTransposed(coord);

    // TODO: we currently only implement the hessian...
    //       a proper implementation will require TMP magic.
    static const unsigned int dim = Base::Traits::GridView::dimensionworld;
    // static_assert(
    //   isHessian<Range>::value,
    //   "We currently only higher order derivative we support is the Hessian of scalar functions");

    // get local hessian of the shape functions
    Range r(0);
    array<std::size_t, dim> directions;
    for(std::size_t i = 0; i < dim; ++i) {
      for(std::size_t j = i; j < dim; ++j) {
        directions[0] = 0;
        directions[1] = 0;
        directions[i]++;
        directions[j]++;
        lfs_.finiteElement().localBasis().evaluate(directions,coord,yb_);
        assert( yb_.size() == 1); // TODO: we only implement the hessian of scalar functions
        for(std::size_t n = 0; n < yb_.size(); ++n) {
          // sum up derivatives, weighting them with the appropriate coeff
          r[i][j] += xl_[i] * yb_[j];
        }
        // use symmetry of the hessian
        if (i != j) r[i][j] = r[j][i];
      }
    }
    // transform back to global coordinates
    for(std::size_t i = 0; i < dim; ++i)
      for(std::size_t j = i; j < dim; ++j)
        r[i][j] *= JgeoIT[i][j] * JgeoIT[i][j];
  }

  friend LocalDerivative derivative(const DiscreteLocalGridViewFunctionDerivative & f)
  {
    return LocalDerivative(f.pgfs_,f.v_);
  }

};

  template<typename GFS, typename V, int N>
class DiscreteLocalGridViewFunctionNoDerivative
  : public DiscreteLocalGridViewFunctionBase<
  Imp::DiscreteGridViewFunctionTraits<GFS,V,N>,
  typename Imp::DiscreteGridViewFunctionTraits<GFS,V,N>::LocalFiniteElementRange
  >
{

  typedef DiscreteLocalGridViewFunctionBase<
    Imp::DiscreteGridViewFunctionTraits<GFS,V,N>,
    typename Imp::DiscreteGridViewFunctionTraits<GFS,V,N>::LocalFiniteElementRange
    > Base;

public:

  typedef typename Base::Domain Domain;
  typedef typename Base::Range Range;

  typedef typename Base::Vector Vector;
  typedef typename Base::GridFunctionSpace GridFunctionSpace;

  using LocalDerivative = DiscreteLocalGridViewFunctionNoDerivative;

  DiscreteLocalGridViewFunctionNoDerivative(const shared_ptr<const GridFunctionSpace> gfs, const shared_ptr<const Vector> v)
    : Base(gfs,v)
  {
    DUNE_THROW(InvalidStateException, N << "th derivative not available");
  }

  friend LocalDerivative derivative(const DiscreteLocalGridViewFunctionNoDerivative &)
  {
    DUNE_THROW(InvalidStateException, N << "th derivative not available, thus you can't call any methods.");
  }

  Range operator()(const Domain& coord) const
  {
    DUNE_THROW(InvalidStateException, N << "th derivative not available, thus you can't call any methods.");
  }

};

template<typename GFS, typename V, int N>
class DiscreteLocalGridViewFunctionDerivative<GFS,V,N,false> :
    public DiscreteLocalGridViewFunctionNoDerivative<GFS,V,N>
{
public:
  using Base = DiscreteLocalGridViewFunctionNoDerivative<GFS,V,N>;
  using Base::Base;
};

template<typename GFS, typename V>
class DiscreteLocalGridViewFunctionDerivative<GFS,V,1,false> :
    public DiscreteLocalGridViewFunctionNoDerivative<GFS,V,1>
{
public:
  using Base = DiscreteLocalGridViewFunctionNoDerivative<GFS,V,1>;
  using Base::Base;
};

} // end of namespace Dune::PDELab
} // end of namespace Dune

#endif // DUNE_PDELAB_FUNCTION_DISCRETEGRIDVIEWFUNCTION_HH
