// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_UTILITY_HH
#define DUNE_PDELAB_COMMON_TYPETREE_UTILITY_HH

#include <dune/common/shared_ptr.hh>

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      /** \addtogroup TypeTree
       *  \ingroup PDELab
       *  \{
       */

#ifndef DOXYGEN

      template<typename T>
      shared_ptr<T> convert_arg(const T& t)
      {
        return make_shared<T>(t);
      }

      template<typename T>
      shared_ptr<T> convert_arg(T& t)
      {
        return stackobject_to_shared_ptr(t);
      }

#if HAVE_RVALUE_REFERENCES

      // only bind to real rvalues
      template<typename T>
      typename enable_if<!std::is_lvalue_reference<T>::value,shared_ptr<T> >::type convert_arg(T&& t)
      {
      return make_shared<T>(std::forward<T>(t));
    }

#endif

#endif // DOXYGEN

      //! \} group TypeTree

    } // namespace TypeTree
  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TYPETREE_UTILITY_HH
