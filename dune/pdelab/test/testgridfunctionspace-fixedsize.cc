// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/finiteelementmap/qkdg.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>

// Wrapped Qk local basis which can have a zero size on some elements
template <class T, int k, int d>
struct WrappedBasis {
  using Basis = Dune::QkStuff::QkLocalBasis<T, T, k, d>;

  Basis basis;
  bool setZero;

  WrappedBasis()
      : setZero(false)
  {
  }

  unsigned int size() const
  {
    return setZero ? 0 : basis.size();
  }

  inline void evaluateFunction(const typename Basis::Traits::DomainType& in,
      std::vector<typename Basis::Traits::RangeType>& out) const
  {
    if (setZero)
      out.clear();
    else
      basis.evaluateFunction(in, out);
  }

  void evaluateJacobian(const typename Basis::Traits::DomainType& in,
      std::vector<typename Basis::Traits::JacobianType>& out) const
  {
    if (setZero)
      out.clear();
    else
      basis.evaluateFunction(in, out);
  }

  unsigned int order() const
  {
    return setZero ? 0 : basis.order();
  }
};

// wrapped QkDG local coefficients which can have a zero size on some elements
template <int k, int d>
struct WrappedCoefficients {
  using Coefficients = Dune::QkStuff::QkDGLocalCoefficients<k, d>;

  Coefficients coeffs;
  bool setZero;

  WrappedCoefficients()
      : setZero(false)
  {
  }

  std::size_t size() const
  {
    return setZero ? 0 : coeffs.size();
  }

  const Dune::LocalKey& localKey(std::size_t i) const
  {
    if (setZero) {
      DUNE_THROW(Dune::Exception, "set zero");
    } else {
      return coeffs.localKey(i);
    }
  }
};

// wrapped Qk local interpolation which can have a zero size on some elements
template <class T, int k, int d>
struct WrappedInterpolation {
  using Interpolation
      = Dune::QkStuff::QkLocalInterpolation<k, d, Dune::QkStuff::QkLocalBasis<T, T, k, d> >;

  Interpolation interpolation;
  bool setZero;

  WrappedInterpolation()
      : setZero(false)
  {
  }

  template <typename F, typename C>
  void interpolate(const F& f, std::vector<C>& out) const
  {
    if (setZero)
      out.clear();
    else
      interpolation.interpolate(f, out);
  }
};

// wrapped Qk local finite element which can have a zero size on some elements
template <class T, int k, int d>
class WrappedFiniteElement
{
  using Basis = WrappedBasis<T, k, d>;
  using Coefficients = WrappedCoefficients<k, d>;
  using Interpolation = WrappedInterpolation<T, k, d>;

public:
  typedef Dune::LocalFiniteElementTraits<Basis, Coefficients, Interpolation> Traits;

  WrappedFiniteElement()
    : gt(Dune::GeometryTypes::cube(d))
  {}

  const typename Traits::LocalBasisType& localBasis() const
  {
    return basis;
  }

  const typename Traits::LocalCoefficientsType& localCoefficients() const
  {
    return coefficients;
  }

  const typename Traits::LocalInterpolationType& localInterpolation() const
  {
    return interpolation;
  }

  Dune::GeometryType type() const
  {
    return gt;
  }

  WrappedFiniteElement* clone() const
  {
    return new WrappedFiniteElement(*this);
  }

  void setZero(bool value)
  {
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

// wrapped QkDG local finite element map which has no local basis functions on elements with an even
// index
template <class GV, class T, int k, int d>
class WrappedFiniteElementMap
    : public Dune::PDELab::
          LocalFiniteElementMapInterface<Dune::PDELab::
                                             LocalFiniteElementMapTraits<WrappedFiniteElement<T, k,
                                                 d> >,
              WrappedFiniteElementMap<GV, T, k, d> >
{

public:
  using Traits = Dune::PDELab::LocalFiniteElementMapTraits<WrappedFiniteElement<T, k, d> >;

  WrappedFiniteElementMap(const GV& gv_, std::function<bool(std::size_t)> zeroPredicate_)
      : gv(gv_), zeroPredicate(zeroPredicate_)
  {
  }

  bool fixedSize() const
  {
    return false;
  }

  bool hasDOFs(int codim) const
  {
    return codim == 0;
  }

  std::size_t size(Dune::GeometryType gt) const
  {
    DUNE_THROW(Dune::PDELab::VariableElementSize, "the fem has variable element size");
  }

  std::size_t maxLocalSize() const
  {
    return Dune::QkStuff::QkSize<k, d>::value;
  }

  template <class EntityType>
  const typename Traits::FiniteElementType& find(const EntityType& e) const
  {
    fe.setZero(zeroPredicate(gv.indexSet().index(e)));
    return fe;
  }

private:
  mutable WrappedFiniteElement<T, k, d> fe;
  GV gv;
  std::function<bool(std::size_t)> zeroPredicate;
};

/**
 * check if the size of the ordering corresponds to the number of basis functions reported by the
 * finite element map
 */
template <class GFS, class FEM>
bool testSize(const GFS& gfs, const FEM& fem)
{
  std::cout << "ordering fixedSize: " << gfs.ordering().fixedSize() << "\n";
  std::cout << "ordering size: " << gfs.ordering().size() << "\n";
  std::cout << "ordering blockCount: " << gfs.ordering().blockCount() << "\n";

  std::size_t size = 0;
  for (const auto& e : Dune::elements(gfs.gridView())) {
    size += fem.find(e).localBasis().size();
  }
  std::cout << "size by summing: " << size << "\n";

  return size == gfs.ordering().size();
}

/**
 * create a grid function space on a structured grid and test wether the size of the ordering is
 * correct
 */
template <int dim>
bool testPredicate(std::size_t cells, bool should_be_fixed_size,  std::function<bool(std::size_t)> predicate)
{
  Dune::FieldVector<double, dim> L(1.0);
  std::array<int, dim> N(Dune::filledArray<dim, int>(cells));
  Dune::YaspGrid<dim> grid(L, N);

  using GV = typename Dune::YaspGrid<dim>::LeafGridView;
  using FEM = WrappedFiniteElementMap<GV, double, 1, dim>;
  using GFS = Dune::PDELab::GridFunctionSpace<GV, FEM>;

  GV gv = grid.leafGridView();
  FEM fem(gv, predicate);
  GFS gfs(gv, fem);
  gfs.update();
  return testSize(gfs, fem) and gfs.ordering().fixedSize() == should_be_fixed_size;
}

int main(int argc, char** argv)
{
  try {
    // Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    const std::size_t cells = 15;

    // In these two basic scenarios, PDELab should detect that the ordering actually is fixed
    // size and switch to fixed size mode after collecting the size information
    std::cout << "testing 2D all cells active" << std::endl;
    if (!testPredicate<2>(cells, true, [](auto i) { return true; })) {
      return -1;
    }
    std::cout << "testing 2D all cells disabled" << std::endl;
    if (!testPredicate<2>(cells, true, [](auto i) { return false; })) {
      return -1;
    }


    // in the following we test different scenarios of a grid function space with a finite element
    // that has 0 or c number of DOFs on different entities, where c is a constant. Wether an entity
    // has 0 DOFs is determined by the predicates provided below.

    std::cout << "testing 1D alternating sizes" << std::endl;
    if (!testPredicate<1>(cells, false, [](auto i) { return i % 2 == 0; })) {
      return -1;
    }
    std::cout << "testing 2D alternating sizes" << std::endl;
    if (!testPredicate<2>(cells, false, [](auto i) { return i % 2 == 0; })) {
      return -1;
    }
    std::cout << "testing 1D empty fem at the end" << std::endl;
    if (!testPredicate<1>(cells, false, [cells](auto i) { return i > cells / 2; })) {
      return -1;
    }
    std::cout << "testing 1D empty fem at the beginning" << std::endl;
    if (!testPredicate<1>(cells, false, [cells](auto i) { return i <= cells / 2; })) {
      return -1;
    }
    return 0;

  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
