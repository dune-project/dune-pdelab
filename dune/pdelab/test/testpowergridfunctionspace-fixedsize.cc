// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/finiteelementmap/qkdg.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/ordering/chunkedblockordering.hh>

// Wrapped Qk local basis which can have a zero size on some elements
template <class T, int k, int d> struct WrappedBasis {
  using Basis = Dune::QkStuff::QkLocalBasis<T, T, k, d>;

  Basis basis;
  bool setZero;

  WrappedBasis() : setZero(false) {}

  unsigned int size() const { return setZero ? 0 : basis.size(); }

  inline void
  evaluateFunction(const typename Basis::Traits::DomainType &in,
                   std::vector<typename Basis::Traits::RangeType> &out) const {
    if (setZero)
      out.clear();
    else
      basis.evaluateFunction(in, out);
  }

  void evaluateJacobian(
      const typename Basis::Traits::DomainType &in,
      std::vector<typename Basis::Traits::JacobianType> &out) const {
    if (setZero)
      out.clear();
    else
      basis.evaluateFunction(in, out);
  }

  unsigned int order() const { return setZero ? 0 : basis.order(); }
};

// wrapped QkDG local coefficients which can have a zero size on some elements
template <int k, int d> struct WrappedCoefficients {
  using Coefficients = Dune::QkStuff::QkDGLocalCoefficients<k, d>;

  Coefficients coeffs;
  bool setZero;

  WrappedCoefficients() : setZero(false) {}

  std::size_t size() const { return setZero ? 0 : coeffs.size(); }

  const Dune::LocalKey &localKey(std::size_t i) const {
    if (setZero) {
      DUNE_THROW(Dune::Exception, "set zero");
    } else {
      return coeffs.localKey(i);
    }
  }
};

// wrapped Qk local interpolation which can have a zero size on some elements
template <class T, int k, int d> struct WrappedInterpolation {
  using Interpolation = Dune::QkStuff::QkLocalInterpolation<
      k, d, Dune::QkStuff::QkLocalBasis<T, T, k, d>>;

  Interpolation interpolation;
  bool setZero;

  WrappedInterpolation() : setZero(false) {}

  template <typename F, typename C>
  void interpolate(const F &f, std::vector<C> &out) const {
    if (setZero)
      out.clear();
    else
      interpolation.interpolate(f, out);
  }
};

// wrapped Qk local finite element which can have a zero size on some elements
template <class T, int k, int d> class WrappedFiniteElement {
  using Basis = WrappedBasis<T, k, d>;
  using Coefficients = WrappedCoefficients<k, d>;
  using Interpolation = WrappedInterpolation<T, k, d>;

public:
  typedef Dune::LocalFiniteElementTraits<Basis, Coefficients, Interpolation>
      Traits;

  WrappedFiniteElement()
    : gt(Dune::GeometryTypes::cube(d))
  {}

  const typename Traits::LocalBasisType &localBasis() const { return basis; }

  const typename Traits::LocalCoefficientsType &localCoefficients() const {
    return coefficients;
  }

  const typename Traits::LocalInterpolationType &localInterpolation() const {
    return interpolation;
  }

  Dune::GeometryType type() const { return gt; }

  WrappedFiniteElement *clone() const {
    return new WrappedFiniteElement(*this);
  }

  void setZero(bool value) {
    basis.setZero = value;
    coefficients.setZero = value;
    interpolation.setZero = value;
  }

private:
  Basis basis;
  Coefficients coefficients;
  Interpolation interpolation;
  Dune::GeometryType gt;
};

// wrapped QkDG local finite element map which has no local basis functions on
// elements with an even
// index
template <class GV, class T, int k, int d>
class WrappedFiniteElementMap
    : public Dune::PDELab::LocalFiniteElementMapInterface<
          Dune::PDELab::LocalFiniteElementMapTraits<
              WrappedFiniteElement<T, k, d>>,
          WrappedFiniteElementMap<GV, T, k, d>> {

public:
  using Traits =
      Dune::PDELab::LocalFiniteElementMapTraits<WrappedFiniteElement<T, k, d>>;

  WrappedFiniteElementMap(const GV &gv_, unsigned int rank_)
      : gv(gv_), rank(rank_) {}

  bool fixedSize() const { return false; }

  bool hasDOFs(int codim) const { return codim == 0; }

  std::size_t size(Dune::GeometryType gt) const {
    DUNE_THROW(Dune::PDELab::VariableElementSize,
               "the fem has variable element size");
  }

  std::size_t maxLocalSize() const {
    return Dune::QkStuff::QkSize<k, d>::value;
  }

  template <class EntityType>
  const typename Traits::FiniteElementType &find(const EntityType &e) const {
    fe.setZero(gv.indexSet().index(e) % rank != 0);
    return fe;
  }

private:
  mutable WrappedFiniteElement<T, k, d> fe;
  GV gv;
  unsigned int rank;
};

int main(int argc, char **argv) {
  try {
    // Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    const int dim = 2;

    Dune::FieldVector<double, dim> L(1.0);
    Dune::array<int, dim> N(Dune::fill_array<int, dim>(2));
    Dune::YaspGrid<dim> grid(L, N);

    const int degree = 1;
    const int blockSize = Dune::QkStuff::QkSize<degree, dim>::value;

    using GV = Dune::YaspGrid<dim>::LeafGridView;
    using FEM = WrappedFiniteElementMap<GV, double, degree, dim>;
    using VBE =
        Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::fixed,
                                          blockSize>;
    using GFS =
        Dune::PDELab::GridFunctionSpace<GV, FEM, Dune::PDELab::NoConstraints,
                                        VBE>;
    using PVBE = Dune::PDELab::istl::VectorBackend<>;
    using POT =
        Dune::PDELab::ordering::Chunked<Dune::PDELab::EntityBlockedOrderingTag>;
    using PGFS = Dune::PDELab::PowerGridFunctionSpace<
        GFS, 2, Dune::PDELab::istl::VectorBackend<>, POT>;

    GV gv = grid.leafGridView();
    FEM fem0(gv, 2);
    FEM fem1(gv, 3);
    GFS gfs0(gv, fem0);
    GFS gfs1(gv, fem1);
    PGFS pgfs(gfs0, gfs1, PVBE(), POT(blockSize));
    pgfs.ordering();
    return 0;
  } catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
