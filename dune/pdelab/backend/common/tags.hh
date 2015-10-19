// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_TAGS_HH
#define DUNE_PDELAB_BACKEND_TAGS_HH

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

    namespace tags {

      using unattached_container DUNE_DEPRECATED_MSG("unattached_containers was moved to the namespace Dune::PDELab::Backend for PDELab 2.4. Please update your code accordingly.") = Backend::unattached_container;
      using attached_container DUNE_DEPRECATED_MSG("unattached_containers was moved to the namespace Dune::PDELab::Backend for PDELab 2.4. Please update your code accordingly.") = Backend::attached_container;

    }

  } // namespace PDELab

} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_TAGS_HH
