// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_UTILITY_HH
#define DUNE_PDELAB_COMMON_UTILITY_HH

#include <algorithm>
#include <cctype>
#include <memory>
#include <string_view>
#include <vector>

#include <dune/common/shared_ptr.hh>

namespace Dune {
  namespace PDELab {


    /** \addtogroup common Common Utilities
     *  \ingroup PDELab
     *  \{
     */

    //! Ensures that t is wrapped in a shared_ptr<T>
    /**
     * You have to consider three situations:
     *
     * a) t is of type T&
     *  t is a stack object and must not be deleted.
     *  You create a shared_ptr<T>(&t) with a null_deleter.
     *  b) t is of type T*
     *  t is a raw pointer and the user os assumed to own this pointer.
     *  You create a shared_ptr<T>(t) with a null_deleter.
     *  c) t is of type shared_ptr<T>
     *  t is already a shared_ptr<T>.
     *  You don't have to do anything.
     */
    template<typename T>
    std::shared_ptr<T> ensure_shared_ptr(T & t)
    {
      return std::shared_ptr<T>(&t, null_deleter<T>());
    }

#ifndef DOXYGEN

    template<typename T>
    std::shared_ptr<T> ensure_shared_ptr(T * t)
    {
      return std::shared_ptr<T>(t, null_deleter<T>());
    }

    template<typename T>
    std::shared_ptr<T> & ensure_shared_ptr(std::shared_ptr<T> & t)
    {
      return t;
    }

#endif // DOXYGEN

    //! Trims a string_view of leading and trailing whitespace.
    std::string_view trim(std::string_view s)
    {
      auto isspace = [](unsigned char c) { return std::isspace(c); };
      auto front = std::find_if_not(begin(s),end(s),isspace);
      auto back = std::find_if_not(rbegin(s),rend(s),isspace).base();
      return front < back ? std::string_view(s.data() + (front - begin(s)),back-front) : std::string_view();
    }

    //! Parses a string with a list of delimited entries into a vector with the individual, trimmed entries.
    std::vector<std::string_view> parseConfigList(std::string_view list, char delim = ',')
    {
      std::vector<std::string_view> result;
      std::size_t front = 0;
      std::size_t back = 0;
      while ((back = list.find(delim,front)) != std::string::npos)
      {
        auto item = trim(list.substr(front,back-front));
        if (not item.empty())
          result.push_back(item);
        front = back + 1;
      }
      auto item = trim(list.substr(front));
      if (not item.empty())
        result.push_back(item);
      return result;
    }



    //! \}

  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_UTILITY_HH
