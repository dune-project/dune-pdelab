// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_COMMON_TAGS_HH
#define DUNE_PDELAB_BACKEND_COMMON_TAGS_HH

#include <dune/common/deprecated.hh>

/** \file
 * \brief Various tags for influencing backend behavior.
 * \ingroup Backend
 */

namespace Dune {

  namespace PDELab{

    namespace Backend {

      //! \addtogroup backend_tags Tags
      //! \brief Tags for controlling behavior within the Backend subsystem.
      //! \ingroup Backend
      //! \{

      //! Tag for requesting a vector or matrix container without a pre-attached underlying object.
      struct unattached_container
      {};

      //! Tag for requesting a vector or matrix container with a pre-attached underlying object.
      struct attached_container
      {};

      //! \} group backend_tags

    }

  } // namespace PDELab

} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_COMMON_TAGS_HH
