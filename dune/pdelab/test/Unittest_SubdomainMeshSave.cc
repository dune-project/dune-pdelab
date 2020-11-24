#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab.hh>

#include <dune/grid/utility/globalindexset.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/common/gridinfo.hh> // visualize grid information -> could be remove at the end


#include <dune/grid/utility/parmetisgridpartitioner.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/SubDomainGmshWriter.hh>
#include <dune/pdelab/backend/istl/geneo/OfflineOnline/partitioner.hh>

/*
 * Defining a Darcy problem with alternating layers of permeability and a high contrast
 */
template<typename GV, typename RF>
class GenericEllipticProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);

    RF perm1 = 1e0;
    RF perm2 = 1e0; // FIXME we want high contrast
    RF layer_thickness = 1.0 / 40.0;

    RF coeff = (int)std::floor(xglobal[0] / layer_thickness) % 2 == 0 ? perm1 : perm2;

    typename Traits::PermTensorType I;
    I[0][0] = coeff;
    I[0][1] = 0;
    I[1][0] = 0;
    I[1][1] = coeff;
    return I;
  }

  //! tensor diffusion constant per cell? return false if you want more than one evaluation of A per cell.
  static constexpr bool permeabilityIsConstantPerCell()
  {
    return true;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType v(0.0);
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 1.0;
  }

  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);
    if (!(xglobal[0]<1E-6 || xglobal[0]>1.0-1E-6))
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    /*if (xglobal[0] > 1.0-1E-6)
      return 1.0;
    else*/
      return 0.0;
  }

  //! flux boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }
};

template <typename V_>
struct AddGatherScatter
{
    static typename V_::value_type gather(const V_ &a, int i)
    {
        return a[i]; // I am sending my value
    }
    static void scatter(V_ &a, typename V_::value_type v, int i)
    {
        a[i] += v; // add what I receive to my value
    }
};

void driver(Dune::MPIHelper& helper) {

  // ~~~~~~~~~~~~~~~~~~
//  Grid set up
  // ~~~~~~~~~~~~~~~~~~
  // define parameters
  const unsigned int dim = 2;
  const unsigned int degree = 1;
  const std::size_t nonzeros = std::pow(2*degree+1,dim);
  typedef double NumberType;

  typedef Dune::UGGrid<dim> GRID;
  Dune::GridFactory<GRID> factory;
  Dune::GmshReader<GRID>::read(factory, "grids/24x24.msh", true, true);
  std::unique_ptr<GRID> grid (factory.createGrid());

  typedef typename GRID::LeafGridView GV;
  auto gv = grid->leafGridView();

  // std::vector<unsigned> partition(grid->leafGridView().size(0));
  // std::cout << grid->leafGridView().size(0) << std::endl;
  // for (int i=0; i< partition.size(); i++) {
  //   if (i < partition.size()/2 ) {partition[i]=0;}
  //   else {partition[i]=1;}
  // }

  unsigned subdomains = 3; // grid->comm().size();
  unsigned uoverlap = 1u;
  std::string extensionmethod = "vertex";

  auto partition = parmetis_partitioning(gv,helper,subdomains);

  // ~~~~~~~~~~~~~~~~~~
  // Save domain decomposition
  // ~~~~~~~~~~~~~~~~~~
  // First recreate the overlapped subdomain
  auto overlappingsubdomains = grow_subdomains(gv,subdomains,partition,uoverlap,extensionmethod);

  // test overlap generation
  auto& indexset = gv.indexSet(); // to attach data to elements
  std::cout << indexset.size(0) << " elements" << std::endl;
  std::cout << indexset.size(dim) << " vertices" << std::endl;
  std::cout << subdomains << " subdomains" << std::endl;
  int vp = 0;
  std::vector<unsigned> adomain(indexset.size(0));
  std::vector<unsigned> k0(indexset.size(0));
  size_t k_0 = 0;
  unsigned count_elements=0;
  for (const auto& e : elements(gv))
  {
    adomain[indexset.index(e)] = 0;
    if (overlappingsubdomains[indexset.index(e)].count(vp)>0) adomain[indexset.index(e)] = 1;
    if (partition[indexset.index(e)]==vp) adomain[indexset.index(e)] = 2;
    k0[indexset.index(e)] = overlappingsubdomains[indexset.index(e)].size();
    k_0 = std::max(k_0,overlappingsubdomains[indexset.index(e)].size());
    count_elements += overlappingsubdomains[indexset.index(e)].size();
  }
  std::cout << k_0 << " k_0" << std::endl;
  std::cout << indexset.size(0) << " elements, " << count_elements << " in all subdomains, ratio=" << ((double)count_elements)/indexset.size(0) << std::endl;

  // Write subdomain as a msh file
  Dune::SubDomainGmshWriter<GV> sdwriter(gv, overlappingsubdomains, subdomains);
  sdwriter.write("", "subdomain.msh");
}

int main(int argc, char **argv)
{
  using Dune::PDELab::Backend::native;

  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc,argv);

  driver(helper);

  return 0;
}
