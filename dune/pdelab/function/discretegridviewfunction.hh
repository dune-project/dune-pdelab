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

#include <dune/functions/gridfunctions/gridviewfunction.hh>

namespace Dune {
namespace PDELab {

namespace Imp {
  template<typename R, typename R2>
  struct DiscreteGridViewFunctionRange;
  template<typename T, int N, typename R2>
  struct DiscreteGridViewFunctionRange<FieldVector<T,N>, R2>
  {
    typedef FieldVector<R2,N> type;
    typedef FieldVector<R2,N> vector;
  };
  template<typename T, typename R2>
  struct DiscreteGridViewFunctionRange<FieldVector<T,1>, R2>
  {
    typedef R2 type;
    typedef FieldVector<R2,1> vector;
  };
  template<typename T, int N, typename R2>
  struct DiscreteGridViewFunctionRange<FieldVector<T,N>, FieldVector<R2,N>>
  {
    typedef FieldVector<R2,N> type;
    typedef FieldVector<R2,N> vector;
  };
};

template<typename GFS, typename V, unsigned int diffOrder = 0>
class DiscreteGridViewFunction
{
public:
  using GridView = typename GFS::Traits::GridView;
  using EntitySet = Functions::GridViewEntitySet<GridView, 0>;

  using Domain = typename EntitySet::GlobalCoordinate;
  using LocalBasisRange = typename GFS::Traits::FiniteElementMap::Traits::FiniteElement::Traits::LocalBasisType::Traits::RangeType;
  using RangeField = typename V::ElementType;
  using Range = typename Imp::DiscreteGridViewFunctionRange<LocalBasisRange, RangeField>::type;
  using RangeVector = typename Imp::DiscreteGridViewFunctionRange<LocalBasisRange, RangeField>::vector;

  using LocalDomain = typename EntitySet::LocalCoordinate;
  using Element = typename EntitySet::Element;

  using Traits =
    Functions::Imp::GridFunctionTraits<Range(Domain), EntitySet,
                                       Functions::DefaultDerivativeTraits, 16>;

  using Basis = GFS;
  using GridFunctionSpace = GFS;
  using Vector = V;

  class LocalFunction
  {
    using LFS = LocalFunctionSpace<GridFunctionSpace>;
    using LFSCache = LFSIndexCache<LFS>;
    using XView = typename Vector::template ConstLocalView<LFSCache>;

  public:

    using GlobalFunction = DiscreteGridViewFunction;
    using Domain = LocalDomain;
    using Range = GlobalFunction::Range;
    using RangeVector = GlobalFunction::RangeVector;
    using Element = GlobalFunction::Element;
    using size_type = std::size_t;

    using GridFunctionSpace = GFS;
    using Vector = V;

    LocalFunction(const shared_ptr<const GridFunctionSpace> gfs, const shared_ptr<const Vector> v)
      : pgfs_(gfs)
      , v_(v)
      , lfs_(*pgfs_)
      , lfs_cache_(lfs_)
      , x_view_(*v_)
      , xl_(pgfs_->maxLocalSize())
      , yb_(pgfs_->maxLocalSize())
      , element_(nullptr)
    {}

    /**
     * \brief Bind LocalFunction to grid element.
     *
     * You must call this method before evaluate()
     * and after changes to the coefficient vector.
     */
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
        DUNE_THROW(InvalidStateException,"can't call localContext on unbound DiscreteGridViewFunction::LocalFunction");
#endif
      return *element_;
    }

    /**
     * \brief Evaluate LocalFunction at bound element.
     *
     * The result of this method is undefined if you did
     * not call bind() beforehand or changed the coefficient
     * vector after the last call to bind(). In the latter case
     * you have to call bind() again in order to make operator()
     * usable.
     */
    std::enable_if<diffOrder == 0, Range>
    operator()(const Domain& coord) const
    {
      RangeVector r(0);
      auto& basis = lfs_.finiteElement().localBasis();
      basis.evaluateFunction(coord,yb_);
      for (size_type i = 0; i < yb_.size(); ++i)
      {
        r.axpy(xl_[i],yb_[i]);
      }
      return r;
    }

    std::enable_if<diffOrder == 1, Jacobian>
    operator()(const Domain& coord) const
    {
      Jacobian r(0);
      // get Jacobian of geometry
      const typename Base::Element::Geometry::JacobianInverseTransposed
        JgeoIT = element_->geometry().jacobianInverseTransposed(coord);

      // get local Jacobians/gradients of the shape functions
      lfs_.finiteElement().localBasis().evaluateJacobian(coord,yb_);

      Range gradphi;
      r = 0;
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

    std::enable_if<diffOrder == 2, Hessian>
    operator()(const Domain& coord) const
    {
      // TODO: we currently require affine geometries.
      if (! element_->geometry().affine())
        DUNE_THROW(NotImplemented, "Due to missing features in the Geometry interface, "
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
      Hessian r = 0;
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

      return r;
    }

    friend typename Traits::LocalFunctionTraits::DerivativeInterface derivative(const LocalFunction& t)
    {
      DUNE_THROW(NotImplemented,"not implemented");
      // shared_ptr<Derivative> diff = make_shared<Derivative>(pgfs_,v_);
      // // TODO: do we really want this?
      // if (element_) diff->bind(*element_);
      // return diff;
    }

  protected:

    const shared_ptr<const GridFunctionSpace> pgfs_;
    const shared_ptr<const Vector> v_;
    LFS lfs_;
    LFSCache lfs_cache_;
    XView x_view_;
    mutable std::vector<RangeVector> xl_;
    mutable std::vector<RangeVector> yb_;
    const Element* element_;
  };

  DiscreteGridViewFunction(const GridFunctionSpace& gfs, const Vector& v)
    : pgfs_(stackobject_to_shared_ptr(gfs)),v_(stackobject_to_shared_ptr(v))
  {}

  DiscreteGridViewFunction(std::shared_ptr<const GridFunctionSpace> pgfs, std::shared_ptr<const Vector> v)
    : pgfs_(pgfs),v_(v)
  {}

  // this is part of the interface in dune-functions
  const Basis& basis() const
  {
    return *pgfs_;
  }
  const GridFunctionSpace& gridFunctionSpace() const
  {
    return *pgfs_;
  }

  const V& dofs() const
  {
    return *v_;
  }

  // TODO: Implement this using hierarchic search
  Range operator() (const Domain& x) const
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  friend typename Traits::DerivativeInterface derivative(const DiscreteGridViewFunction& t)
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  friend LocalFunction localFunction(const DiscreteGridViewFunction& t)
  {
    return LocalFunction(t.pgfs_, t.v_);
  }

  /**
   * \brief Get associated EntitySet
   */
  EntitySet entitySet() const
  {
    return pgfs_->gridView();
  }

private:

  const shared_ptr<const GridFunctionSpace> pgfs_;
  const shared_ptr<const Vector> v_;

};

} // end of namespace Dune::PDELab
} // end of namespace Dune

#endif // DUNE_PDELAB_FUNCTION_DISCRETEGRIDVIEWFUNCTION_HH
