// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <fenv.h>

#include <dune/localfunctions/test/test-localfe.hh>

#include <dune/pdelab/finiteelement/pk1d.hh>

int main (int argc, char *argv[])
{
#if __linux__ \
  && (!defined __INTEL_COMPILER || __INTEL_COMPILER >= 1010) \
  && (!defined __clang__ || __clang_major__ >= 4)
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

  bool success = true;

  for (int i=1; i<5; i++)
  {
    Dune::Pk1dLocalFiniteElement<double,double> pk1d(i);
    TEST_FE(pk1d);
  }

  return success ? 0 : 1;
}
