#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <memory>
#include <string>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/common/vtkexport.hh>


namespace Utility
{

  template<class GV, class GFS, class V>
  void write_gridfunction_vtk(const std::string& filename, const GV& gv,
                              const GFS& gfs, const V& dofVec, const int ord)
  {
    using DGF = Dune::PDELab::DiscreteGridFunction<GFS, V>;
    using VTKF = Dune::PDELab::VTKGridFunctionAdapter<DGF>;
    using VTKWRITER = Dune::SubsamplingVTKWriter<GV>;

    DGF dofVecDGF(gfs, dofVec);

    // prepare the VTKWriter and write to file
    int subsampling { ord };
    VTKWRITER vtkwriter(gv, Dune::refinementIntervals(subsampling));

    std::string outputname { "gridfunction" };

    vtkwriter.addVertexData(
      std::shared_ptr<VTKF>(new VTKF(dofVecDGF, outputname)));
    vtkwriter.write(filename, Dune::VTK::ascii);

    return;
  }

} // namespace Utility
