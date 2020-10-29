#ifndef GRIDPARAMETERS_HH
#define GRIDPARAMETERS_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <bitset>
#include <array>

#include <dune/common/parametertree.hh>


namespace Utility
{

  /* The grid parameters are needed several times throughout the testcodes for
   * GenEO. We use this container for easier reading and passing of arguments.
   * It is only supposed to be used in 2 or 3 dimensions, so only the
   * corresponding template specializations can be used.
   */

  template<typename GVCoord, int DIM>
  class GridParameters
  {
  public:

    // only the 2 and 3 dimensional versions should be createable
    GridParameters() = delete;
  };


  template<typename GVCoord>
  class GridParameters<GVCoord, 2>
  {
  public:

    const int ovlp;
    const Dune::FieldVector<GVCoord, 2> upperRight;
    const std::bitset<2> isPeriodic;
    const std::array<int, 2> nCells;
    const std::array<int, 2> subdomLayout;

    GridParameters(const Dune::ParameterTree& ptree) :
      ovlp { ptree.get<int>("grid.overlap") },
      upperRight { ptree.get<double>("grid.LX"),
                   ptree.get<double>("grid.LY") },
      isPeriodic("00"),
      nCells { ptree.get<int>("grid.nCellsX"),
               ptree.get<int>("grid.nCellsY") },
      subdomLayout { ptree.get<int>("grid.nSubdomX"),
                     ptree.get<int>("grid.nSubdomY") }
    {}
  };


  // partial template specialization for 3 dimensions
  template<typename GVCoord>
  class GridParameters<GVCoord, 3>
  {
  public:

    const int ovlp;
    const Dune::FieldVector<GVCoord, 2> upperRight;
    const std::bitset<3> isPeriodic;
    const std::array<int, 3> nCells;
    const std::array<int, 3> subdomLayout;

    GridParameters(const Dune::ParameterTree& ptree) :
      ovlp { ptree.get<int>("grid.overlap") },
      upperRight { ptree.get<double>("grid.LX"),
                   ptree.get<double>("grid.LY"),
                   ptree.get<double>("grid.LZ") },
      isPeriodic("000"),
      nCells { ptree.get<int>("grid.nCellsX"),
               ptree.get<int>("grid.nCellsY"),
               ptree.get<int>("grid.nCellsZ") },
      subdomLayout { ptree.get<int>("grid.nSubdomX"),
                     ptree.get<int>("grid.nSubdomY"),
                     ptree.get<int>("grid.nSubdomZ") }
    {}
  };

} // namespace Utility

#endif
