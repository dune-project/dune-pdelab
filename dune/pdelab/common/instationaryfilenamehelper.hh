// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_INSTATIONARYFILENAMEHELPER_HH
#define DUNE_PDELAB_COMMON_INSTATIONARYFILENAMEHELPER_HH

#include <iostream>

namespace Dune {
  namespace PDELab {

    /**
     *  @addtogroup OneStepMethod
     *  @{
     */
    class FilenameHelper
    {
    public:
      FilenameHelper(const char *basename_, int i_=0)
        : i(i_)
      {
        sprintf(basename,"%s",basename_);
      }

      FilenameHelper(const std::string & basename_, int i_=0)
        : i(i_)
      {
        sprintf(basename,"%s",basename_.c_str());
      }

      const char *getName (int i_)
      {
        sprintf(fname,"%s-%05d",basename,i_);
        return fname;
      }

      const char *getName ()
      {
        sprintf(fname,"%s-%05d",basename,i);
        return fname;
      }

      void increment ()
      {
        i++;
      }

    private:
      char fname[255];
      char basename[255];
      int i;
    }; // end class FilenameHelper

    /** @} */
  } // end namespace PDELab
} // end namespace Dune
#endif // DUNE_PDELAB_COMMON_INSTATIONARYFILENAMEHELPER_HH
