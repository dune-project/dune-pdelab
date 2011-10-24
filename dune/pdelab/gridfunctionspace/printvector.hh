// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_PRINTVECTOR_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_PRINTVECTOR_HH

#include <iomanip>
#include <iosfwd>
#include <ostream>
#include <string>

#include <unistd.h>

#include <dune/common/fvector.hh>
#include <dune/common/ios_state.hh>

#include <dune/pdelab/backend/backendselector.hh>
#include <dune/pdelab/gridfunctionspace/dofinfo.hh>

namespace Dune {
  namespace PDELab {

    //! Print a vector with associciated dof positions
    /**
\code
#include <dune/pdelab/gridfunctionspace/printvector.hh>
\endcode
     * This prints a vector in the following form:
\verbatim
title " [rank=" rank ", size=" size "]" NEWLINE
position TAB vector entry NEWLINE
.
.
.
\endverbatim
     * The coordinate components in "position" are printed in fix-width,
     * seperated by a single space, with the precision as determined by the \c
     * coordPrecision parameter.  The "vector entry" is printed in whatever
     * format the output stream was before entering this function.
     *
     * When running in parallel, each process will print its part of the
     * matrix, starting with rank 0.  To avoid mixing the output, the
     * processes will syncronize using barrier(), once for each process before
     * it starts writing and once more after all processes have finished.
     * This is sometimes not enough however, since the unix tty can still read
     * the output of different processes out-of-order, even though these
     * processes have flushed their buffers and syncronized.  To avoid this,
     * the processes can sleep for some time (default 1 second) after each
     * barrier() operation.
     *
     * \param outStream      The stream to write output to.  It should be set
     *                       up to correctly format vector entries.
     * \param gfs            The GridFunctionSpace the vector came from.  Used
     *                       to extract the "positions" of the dofs.
     * \param v              The vector to print.
     * \param title          A title to print before the vector.
     * \param coordPrecision Precision used to print the coordinates.
     * \param delay          How many seconds to wait after each barrier
     *                       operation to give the tty time to settle.
     *
     * \note Currently, this function does not work with trees of
     *       GridFunctionSpaces.  Even if it would work, there would be the
     *       possibility of multiple dofs at the same position.
     * \note Currently, this function does not distinguish multiple dofs at
     *       the same entity.
     */
    template<class GFS, class Vector>
    void printVector(std::ostream &outStream, const GFS& gfs, const Vector &v,
                     const std::string &title, int coordPrecision = 3,
                     unsigned delay = 1)
    {
      ios_base_all_saver stateSaver(outStream);

      typedef typename GFS::Traits::GridViewType GV;
      const GV &gv = gfs.gridview();
      typedef FieldVector<typename GV::ctype, GV::dimensionworld> DomainW;

      typedef typename Dune::PDELab::BackendVectorSelector<GFS, DomainW>::Type
        PosV;
      PosV posV(gfs);
      getDofEntityPositions(gfs, posV);

      typedef typename GFS::Traits::BackendType Backend;

      outStream << std::flush;

      for(int p = 0; p < gv.comm().size(); ++p) {
        gv.comm().barrier();
        sleep(delay);
        if(p != gv.comm().rank()) continue;

        outStream << title << " [rank=" << p << ", "
                  << "size=" << gfs.size() << "]\n"
                  << "pos\tvalue\n";

        for(typename GFS::Traits::SizeType i = 0; i < gfs.size(); ++i) {
          const DomainW &pos = Backend::access(posV, i);
          outStream << std::fixed << std::setprecision(coordPrecision)
                    << pos << "\t";
          stateSaver.restore();
          outStream << Backend::access(v, i) << "\n";
        }
        outStream << std::flush;
      }
      gv.comm().barrier();
      sleep(delay);
    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_PRINTVECTOR_HH
