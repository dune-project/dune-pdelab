// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDOPERATORSPACE_PRINTMATRIX_HH
#define DUNE_PDELAB_GRIDOPERATORSPACE_PRINTMATRIX_HH

#include <iomanip>
#include <iosfwd>
#include <ostream>
#include <string>

#include <unistd.h>

#include <dune/common/fvector.hh>
#include <dune/common/ios_state.hh>

#include <dune/pdelab/gridfunctionspace/dofinfo.hh>

namespace Dune {
  namespace PDELab {

    //! Print a matrix in sparse format with associciated dof positions
    /**
     * This prints a matrix in the following form:
     * \pre
title " (rank=" rank ")" NEWLINE
row position TAB column position TAB matrix entry NEWLINE
.
.
.
     * \endpre
     * The coordinate components in "row position" and "column position" are
     * printed in fix-width, seperated by a single space, with the precision
     * as determined by the \c coordPrecision parameter.  The "matrix entry"
     * is printed in whatever format the output stream was before entering
     * this function.
     *
     * When running in parallel, each process will print it's part of the
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
     *                       up to correctly format matrix entries.
     * \param gos            The GridOperatorSpace the matrix came from.  Used
     *                       to extract the GridFunctionSpaces and the
     *                       "positions" of the dofs.
     * \param m              The matrix to print.
     * \param title          A title to print before the matrix.
     * \param coordPrecision Precision used to print the coordinates.
     * \param delay          How many seconds to wait after each barrier
     *                       operation to give the tty time to settle.
     *
     * \note Currently, this function is restricted to matrices from the
     *       ISTLBCRSMatrixBackend with a blocksize of 1.
     * \note Currently, this function does not work with trees of
     *       GridFunctionSpaces.  Even if it would work, there would be the
     *       possibility of multiple dofs at the same position.
     * \note Currently, this function does not distinguish multiple dofs at
     *       the same entity.
     */
    template<class GOS, class Matrix>
    void printMatrix(std::ostream &outStream, const GOS& gos, const Matrix &m,
                     const std::string &title, int coordPrecision = 3,
                     unsigned delay = 1)
    {
      ios_base_all_saver stateSaver(outStream);

      typedef typename GOS::Traits::TrialGridFunctionSpace TrialGFS;
      const TrialGFS &trialGFS = gos.trialGridFunctionSpace();
      typedef typename GOS::Traits::TestGridFunctionSpace TestGFS;
      const TestGFS &testGFS = gos.testGridFunctionSpace();

      typedef typename TrialGFS::Traits::GridViewType GV;
      const GV &gv = trialGFS.gridview();
      typedef FieldVector<typename GV::ctype, GV::dimensionworld> DomainW;

      typedef typename TrialGFS::template VectorContainer<DomainW>::Type
        TrialPosV;
      TrialPosV trialPosV(trialGFS);
      getDofEntityPositions(trialGFS, trialPosV);
      typedef typename TestGFS::template VectorContainer<DomainW>::Type
        TestPosV;
      TestPosV testPosV(testGFS);
      getDofEntityPositions(testGFS, testPosV);

      typedef typename Matrix::BaseT ISTLMatrix;
      const ISTLMatrix &istlMatrix = m.base();
      typedef typename ISTLMatrix::ConstRowIterator RowIterator;
      typedef typename ISTLMatrix::ConstColIterator ColIterator;

      outStream << std::flush;

      for(int p = 0; p < gv.comm().size(); ++p) {
        gv.comm().barrier();
        sleep(delay);
        if(p != gv.comm().rank()) continue;

        outStream << title << "(rank=" << p << ")\n"
                  << "row pos\tcol pos\tvalue\n";

        const RowIterator &rowEnd = istlMatrix.end();
        for(RowIterator rowIt = istlMatrix.begin(); rowIt != rowEnd; ++rowIt)
        {
          const DomainW &rowPos = testPosV[rowIt.index()];
          const ColIterator &colEnd = rowIt->end();
          for(ColIterator colIt = rowIt->begin(); colIt != colEnd;
              ++colIt)
          {
            const DomainW &colPos = trialPosV[colIt.index()];
            outStream << std::fixed << std::setprecision(coordPrecision)
                      << rowPos << "\t" << colPos << "\t";
            stateSaver.restore();
            outStream << *colIt << "\n";
          }
        }
        outStream << std::flush;
      }
      gv.comm().barrier();
      sleep(delay);
    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDOPERATORSPACE_PRINTMATRIX_HH
