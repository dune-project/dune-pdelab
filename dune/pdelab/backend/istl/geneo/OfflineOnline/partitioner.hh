#ifndef Udune_ftworth_partitioner_HH
#define Udune_ftworth_partitioner_HH

#if HAVE_PARMETIS

#include <parmetis.h>

#include <algorithm>
#include <vector>
#include <dune/grid/utility/globalindexset.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>

// only enable for ParMETIS because the implementation uses functions that
// are not emulated by scotch
#ifdef PARMETIS_MAJOR_VERSION

/** \brief Returns a vector with a partition number for each element
 *
 * \author Benjamin Bykowski and adapted by Peter Bastian
 *
 * \param gv The grid view to be partitioned
 * \param mpihelper The MPIHelper object, needed to get the MPI communicator. This is needed by the function ParMETIS_V3_PartMeshKway and can unfortunately not be omitted
 * \param parts number of subdomains desired
 *
 * \return std::vector with one uint per All_Partition element.  For each element, the entry is the
 *    number of the partition the element is assigned to. This number is greater or equal zero and smaller as parts.
 *    No partitioning is done, only this vector is computed
 */
template<class GridView>
std::vector<unsigned> parmetis_partitioning (const GridView& gv, const Dune::MPIHelper& mpihelper, int parts) {

#if PARMETIS_MAJOR_VERSION > 3
  typedef idx_t idx_type;
  typedef ::real_t real_type;
#else
  typedef int idx_type;
  typedef float real_type;
#endif // PARMETIS_MAJOR_VERSION > 3

  const unsigned numElements = gv.size(0);

  std::vector<unsigned> part(numElements);

  // Setup parameters for ParMETIS
  idx_type wgtflag = 0;                                  // we don't use weights
  idx_type numflag = 0;                                  // we are using C-style arrays
  idx_type ncon = 1;                                     // number of balance constraints
  idx_type ncommonnodes = 2;                             // number of nodes elements must have in common to be considered adjacent to each other
  idx_type options[4] = {0, 0, 0, 0};                    // use default values for random seed, output and coupling
  idx_type edgecut;                                      // will store number of edges cut by partition
  idx_type nparts = parts;                               // number of partitions to create is a parameter
  std::vector<real_type> tpwgts(ncon*nparts, 1./nparts); // load per subdomain and weight (same load on every process)
  std::vector<real_type> ubvec(ncon, 1.05);              // weight tolerance (same weight tolerance for every weight there is)

  // The difference elmdist[i+1] - elmdist[i] is the number of nodes that are on process i
  std::vector<idx_type> elmdist(nparts+1);
  elmdist[0] = 0;
  std::fill(elmdist.begin()+1, elmdist.end(), gv.size(0)); // all elements are on process zero

  // Create and fill arrays "eptr", where eptr[i] is the number of vertices that belong to the i-th element, and
  // "eind" contains the vertex-numbers of the i-the element in eind[eptr[i]] to eind[eptr[i+1]-1]
  std::vector<idx_type> eptr, eind;
  int numVertices = 0;
  eptr.push_back(numVertices);

  for (const auto& element : elements(gv, Dune::Partitions::interior)) {
    const size_t curNumVertices = Dune::referenceElement<double,GridView::dimension>(element.type()).size(GridView::dimension);

    numVertices += curNumVertices;
    eptr.push_back(numVertices);

    for (size_t k = 0; k < curNumVertices; ++k)
      eind.push_back(gv.indexSet().subIndex(element, k, GridView::dimension));
  }

  // Partition mesh using ParMETIS
  if (0 == mpihelper.rank()) {
    MPI_Comm comm = Dune::MPIHelper::getLocalCommunicator();

#if PARMETIS_MAJOR_VERSION >= 4
    const int OK =
#endif
      ParMETIS_V3_PartMeshKway(elmdist.data(), eptr.data(), eind.data(), NULL, &wgtflag, &numflag,
                               &ncon, &ncommonnodes, &nparts, tpwgts.data(), ubvec.data(),
                               options, &edgecut, reinterpret_cast<idx_type*>(part.data()), &comm);

#if PARMETIS_MAJOR_VERSION >= 4
    if (OK != METIS_OK)
      DUNE_THROW(Dune::Exception, "ParMETIS returned an error code.");
#endif
  }

  return part;
}

#else // PARMETIS_MAJOR_VERSION
#warning "You seem to be using the ParMETIS emulation layer of scotch, which does not work with this file."
#endif

#else // HAVE_PARMETIS
#warning "PARMETIS was not found, please check your configuration"
#endif

/**
 * \brief Takes a partition and extends all subdomains by a number of layers given by overlap
 *
 * \param gv The grid view to be operated on
 * \param partitions number of partitions (or subdomains) the grid has been partitioned into (using the function parmetis_partition)
 * \param partition partition information for each element
 * \param overlap overlap that should be added (zero is fine, then nothing is done)
 * \param mode determines how partitions are grown. Can have the value "vertex" or "element", default is "vertex"
 *
 * If mode has the value "element" the extension is done via element faces. Partitions are extended to neighbors of elements in overlap rounds.
 * If mode has the value "vertex" the extension is done via vertices of the grid. First, partitioning is converted from elements to
 * vertices. Then partitions are extended to neighbors of vertices in overlap rounds. Finally, partitioning
 * is converted back to elements.
 * If mode has neither the value "element" or "vertex" the original partitioning is returned, converted to a set for each element
 *
 * \return std::vector<std::set<unsigned>> for each element stores a set of subomains containing this element
 */
template<class GV>
std::vector<std::set<unsigned>> grow_subdomains (const GV& gv, unsigned partitions, std::vector<unsigned> partition, unsigned overlap, std::string mode="vertex")
{
  const int dim = GV::dimension; // extract dimension (codim of vertices)
  auto& indexset = gv.indexSet(); // to attach data to elements

  if (mode=="vertex")
    {
      std::vector<std::set<unsigned>> subdomainsv(indexset.size(dim)); // set of subdomains for each vertex

      // initialize subdomain list for each vertex by the partition
      for (const auto& e : elements(gv))
        for (unsigned int i=0; i<e.subEntities(dim); ++i)
          {
            auto v = e.template subEntity<dim>(i);
            subdomainsv[indexset.index(v)].insert(partition[indexset.index(e)]);
          }

      // in each round extend overlap by one
      for (int rounds=0; rounds<overlap; rounds++)
      {
        std::vector<std::set<unsigned>> old(subdomainsv); // copy current state
        for (const auto& e : elements(gv))
          {
            // build union of all partitions in all vertices of the element
            std::set<unsigned> unification;
            for (unsigned int i=0; i<e.subEntities(dim); ++i)
              for (const auto& j : old[indexset.index(e.template subEntity<dim>(i))])
                unification.insert(j);
            // now add union to all vertices (a clique)
            for (const auto& j : unification)
              for (unsigned int i=0; i<e.subEntities(dim); ++i)
                subdomainsv[indexset.index(e.template subEntity<dim>(i))].insert(j);
          }
      }

      // now convert again to elements: element is in subdomain if *all* vertices are in subdomain
      std::vector<std::set<unsigned>> subdomainse(indexset.size(0)); // set of subdomains for each element
      for (const auto& e : elements(gv))
      {
        std::set<unsigned> intersection(subdomainsv[indexset.index(e.template subEntity<dim>(0))]);
        for (unsigned int i=1; i<e.subEntities(dim); ++i)
          {
            std::set<unsigned> update;
            for (const auto& j : subdomainsv[indexset.index(e.template subEntity<dim>(i))])
        if (intersection.count(j)>0) update.insert(j);
            intersection = update;
          }
        subdomainse[indexset.index(e)] = intersection;
      }
      // and we are done
      return subdomainse;
    }

  // now the element mode
  std::vector<std::set<unsigned>> subdomains(indexset.size(0)); // set of subdomains for each element

  // initialize subdomain list for each element by the partition
  for (const auto& e : elements(gv))
    subdomains[indexset.index(e)].insert(partition[indexset.index(e)]);

  if (mode=="element")
    {
      // in each round extend overlap by one
      for (int rounds=0; rounds<overlap; rounds++)
      {
        std::vector<std::set<unsigned>> old(subdomains); // copy current state
        for (const auto& e : elements(gv))
          for (const auto& is : intersections(gv,e))
            if (is.neighbor())
        for (const auto& i : old[indexset.index(is.outside())])
          subdomains[indexset.index(e)].insert(i);
      }
    }
  // and we are done
  return subdomains;
}



#endif // Udune_ftworth_HH
