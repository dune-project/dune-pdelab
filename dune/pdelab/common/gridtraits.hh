// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_GRIDTRAITS_HH
#define DUNE_PDELAB_COMMON_GRIDTRAITS_HH

namespace Dune {
  namespace PDELab {

    /** \brief Specialize with 'true' if implementation requires communication when running on a single process
        \ingroup PDELabGridCapabilities
     */
    template<class Grid>
    struct requiresCommOnSequential
    {
      static bool v(const Grid& grid) {
        return false;
      };
    };

    /** \brief YaspGrid needs communication for periodic boundaries to work (even when running sequentially)
        \ingroup PDELabGridCapabilities
     */
    template<int dim, class Coordinates>
    struct requiresCommOnSequential< YaspGrid<dim, Coordinates> >
    {
      static bool v(const YaspGrid<dim, Coordinates>& grid) {
        for (int i = 0; i < dim; ++i) {
          if (grid.isPeriodic(i))
            return true;
        }
        return false;
      };
    };

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_COMMON_GRIDTRAITS_HH
