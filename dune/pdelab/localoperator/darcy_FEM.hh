// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_DARCY_FEM_HH
#define DUNE_PDELAB_LOCALOPERATOR_DARCY_FEM_HH

#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

//! \brief Provide Darcy velocity as a vector-valued grid function
/**
 * The function values should be single-component vectors.  The Gradient
 * will be a dimDomain-component function.
 *
 * \tparam T Type of GridFunctionSpace.  The LocalBasis must provide the
 *           evaluateJacobian() method.
 * \tparam X Type of coefficients vector
 */
template<typename P, typename T, typename X>
class DarcyVelocityFromHeadFEM
  : public Dune::PDELab::GridFunctionInterface<
  Dune::PDELab::GridFunctionTraits<
    typename T::Traits::GridViewType,
    typename T::Traits::FiniteElementType::Traits::LocalBasisType
    ::Traits::RangeFieldType,
    T::Traits::FiniteElementType::Traits::LocalBasisType::Traits
    ::dimDomain,
    Dune::FieldVector<
      typename T::Traits::FiniteElementType::Traits
      ::LocalBasisType::Traits::RangeFieldType,
      T::Traits::FiniteElementType::Traits::LocalBasisType::Traits
      ::dimDomain> >,
  DarcyVelocityFromHeadFEM<P,T,X> >
{
  typedef T GFS;
  typedef typename GFS::Traits::FiniteElementType::Traits::
  LocalBasisType::Traits LBTraits;

  typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
  typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
  typedef typename X::template LocalView<LFSCache> LView;

public:
  typedef Dune::PDELab::GridFunctionTraits<
  typename GFS::Traits::GridViewType,
  typename LBTraits::RangeFieldType,
  LBTraits::dimDomain,
  Dune::FieldVector<
    typename LBTraits::RangeFieldType,
    LBTraits::dimDomain> > Traits;

private:
  typedef Dune::PDELab::GridFunctionInterface<
  Traits,
  DarcyVelocityFromHeadFEM<P,T,X> > BaseT;

public:
  /** \brief Construct a DarcyVelocityFromHeadFEM
   *
   * \param gfs The GridFunctionsSpace
   * \param x_  The coefficients vector
   */
  DarcyVelocityFromHeadFEM (const P& p, const GFS& gfs, X& x_)
    : pgfs(stackobject_to_shared_ptr(gfs))
    , pxg(stackobject_to_shared_ptr(x_))
    , pp(stackobject_to_shared_ptr(p))
    , lfs(pgfs)
    , lfs_cache(lfs)
    , lview(*pxg)
  {}

  // Evaluate
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    // get and bind local functions space
    lfs.bind(e);
    lfs_cache.update();

    // get local coefficients
    std::vector<typename Traits::RangeFieldType> xl(lfs.size());
    lview.bind(lfs_cache);
    lview.read(xl);
    lview.unbind();

    // get Jacobian of geometry
    const typename Traits::ElementType::Geometry::JacobianInverseTransposed
      JgeoIT(e.geometry().jacobianInverseTransposed(x));

    // get local Jacobians/gradients of the shape functions
    std::vector<typename LBTraits::JacobianType> J(lfs.size());
    lfs.finiteElement().localBasis().evaluateJacobian(x,J);

    typename Traits::RangeType gradphi;
    typename Traits::RangeType minusgrad(0);
    for(unsigned int i = 0; i < lfs.size(); ++i) {
      // compute global gradient of shape function i
      gradphi = 0;
      JgeoIT.umv(J[i][0], gradphi);

      // sum up global gradients, weighting them with the appropriate coeff
      minusgrad.axpy(-xl[i], gradphi);
    }

    // multiply with permeability tensor
    typedef typename Traits::DomainFieldType DF;
    const int dim = LBTraits::dimDomain;
    const Dune::FieldVector<DF,dim>
      inside_cell_center_local = Dune::ReferenceElements<DF,dim>::
      general(e.type()).position(0,0);
    typename P::Traits::PermTensorType A(pp->A(e,inside_cell_center_local));
    A.mv(minusgrad,y);
  }

  //! get a reference to the GridView
  inline const typename Traits::GridViewType& getGridView () const
  {
    return pgfs->gridView();
  }

private:
  std::shared_ptr<const GFS> pgfs;
  std::shared_ptr<X> pxg;
  std::shared_ptr<const P> pp;
  mutable LFS lfs;
  mutable LFSCache lfs_cache;
  mutable LView lview;
};

#endif
