// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_VTK_PVDWRITER_HH
#define DUNE_VTK_PVDWRITER_HH

#include <vector>
#include <fstream>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include "onestep.hh"

#include <stdlib.h>

namespace Dune
{

    template< class GridView, class VTK = VTKWriter<GridView> >
    class PVDWriter : public VTK
    {
        const GridView & gv;
        std::string basename;
        PDELab::FilenameHelper fn;
        std::string path;
        std::vector<double> timesteps;
        Dune::VTK::OutputType outputtype;
        unsigned int offset;

    public:

        PVDWriter(const GridView & gv_, std::string basename_,
                  Dune::VTK::DataMode datamode_ = Dune::VTK::conforming,
                  Dune::VTK::OutputType outputtype_ = Dune::VTK::appendedraw,
                  std::string path_="vtk", unsigned int offset_=0) :
            VTK(gv_,datamode_), gv(gv_),
            basename(basename_), fn(basename_,offset_),
            path(path_), outputtype(outputtype_),
            offset(offset_){}

        void write(double time)
        {
            /* remember current time step */
            timesteps.push_back(time);
            /* make sure the directory exists */
            // mkdir("vtk", 777);
            /* write VTK file */
            VTK::pwrite(fn.getName(),path,"",outputtype);
            /* write pvd file */
            std::string pvdname = basename + ".pvd";
            std::ofstream pvd(pvdname.c_str());
            //std::cout << "WRITE PVD FILE " << pvdname << std::endl;
            assert(pvd.is_open());
            pvd << std::fixed;
            pvd << "<?xml version=\"1.0\"?>\n"
                << "<VTKFile type=\"Collection\" version=\"0.1\">\n"
                << "<Collection>\n";
            PDELab::FilenameHelper fnloop(basename,offset);
            for (unsigned int i=0; i<timesteps.size(); i++)
                {
                    std::string fname = this->getParallelHeaderName(fnloop.getName(), path, gv.comm().size());
                    pvd << "  <DataSet timestep=\"" << timesteps[i]
                        << "\" file=\"" << fname << "\"/>\n";
                    fnloop.increment();
                }
            pvd << "</Collection>\n"
                << "</VTKFile>\n";
            pvd.close();

            /* increment counter */
            fn.increment();
        }
    };

} // end namespace Dune

#endif // DUNE_VTK_PVDWRITER_HH
