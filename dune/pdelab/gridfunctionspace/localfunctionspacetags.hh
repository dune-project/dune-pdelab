// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_LOCALFUNCTIONSPACETAGS_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_LOCALFUNCTIONSPACETAGS_HH

/** \file
    \author Christian Engwer
    Define tags to mark local trial and test spaces.
 */

namespace Dune {
  namespace PDELab {
    
    /**
       \addtogroup PDELAB_StrictTrialAndTest Strict Trial and Test space handling
       \ingroup GridFunctionSpace
       \{
       
       When compiling a program with 
       \code
       -DDUNE_PDELAB_CHECK_LOCALINDEX
       \endcode

       additional checks are enabled. For these checks the
       GridFunctionSpace::LocalFunctionSpace types are either tagged
       with TrialSpaceTag or with TestSpaceTag. Usually they are
       tagged with AnySpaceTag, which means that no additional test
       are performed.

       To make code aware of the difference between trial and test space, compile with
       \code
       -DDUNE_PDELAB_CHECK_LOCALINDEX_STRICT_TRIAL_AND_TEST_SPACE
       \endcode
       note, that this will report errors for all local operators, where lfsu==lfsv is assumed.
     */

    /**
       The I-don't-care-which-space-it-is tag
     */
    struct AnySpaceTag {};

#if defined DOXYGEN \
    || defined DUNE_PDELAB_CHECK_LOCALINDEX \
    || defined DUNE_PDELAB_CHECK_LOCALINDEX_STRICT_TRIAL_AND_TEST_SPACE
    /**
       Tag to mark a LocalFunctionSpace as trial-space
     */
    struct TrialSpaceTag {};
#if defined DOXYGEN \
    || defined DUNE_PDELAB_CHECK_LOCALINDEX_STRICT_TRIAL_AND_TEST_SPACE
    /**
       Tag to mark a LocalFunctionSpace as test-space
     */
    struct TestSpaceTag {};
#else
    typedef TrialSpaceTag TestSpaceTag;
#endif    
#else
    typedef AnySpaceTag TrialSpaceTag;
    typedef AnySpaceTag TestSpaceTag;
#endif

    /**
       \}
     */

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_LOCALFUNCTIONSPACETAGS_HH
