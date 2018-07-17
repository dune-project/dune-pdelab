// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_FLAVOR_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_FLAVOR_HH

namespace Dune {
  namespace PDELab {

    namespace Flavor {

      struct Generic
      {};

      struct Inside
      {};

      struct Outside
      {};

      struct Test
      {};

      struct Trial
      {};

      struct InsideTest
        : public Inside
        , public Test
      {};

      struct InsideTrial
        : public Inside
        , public Trial
      {};

      struct OutsideTest
        : public Outside
        , public Test
      {};

      struct OutsideTrial
        : public Outside
        , public Trial
      {};

    } // namespace Flavor
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_FLAVOR_HH
