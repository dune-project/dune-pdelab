// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=8 sw=2 et sts=2:
#ifndef DUNE_PDELAB_COMMON_IOSTATESAVER_HH
#define DUNE_PDELAB_COMMON_IOSTATESAVER_HH

#include <ios>

namespace Dune {
  namespace PDELab {

    //! Exception safe temporary modifications to stream format attributes
    /**
     * This class was mostly stolen from boost.ios_state.  It is intended to
     * be drop-in compatible with the class from boost.ios_state.
     */
    class ios_base_all_saver {
    public:
      //! export base class of the stream holding format information
      typedef ::std::ios_base state_type;

      //! constructor
      /**
       * \param s_ The stream to save the format attributes for.  A reference
       *           to this objet is saved, so the object must be valid until
       *           the destructor of this saver class completes.
       */
      explicit inline ios_base_all_saver(state_type &s_)
        : s(s_), flags(s.flags()), precision(s.precision())
        , width(s.width())
      {}

      //! destructor
      /**
       * Calls restore() to restore the format attributes.
       */
      inline ~ios_base_all_saver()
      { this->restore(); }

      //! restore format attributes
      /**
       * Restore the format attributes of the stream to the state they had at
       * construction of this saver object.  This method maya be called as
       * often as desired.
       */
      inline void restore()
      {
        s.width(width);
        s.precision(precision);
        s.flags(flags);
      }

    private:
      state_type &                s;
      state_type::fmtflags const  flags;
      ::std::streamsize const     precision;
      ::std::streamsize const     width;

      ios_base_all_saver& operator=(const ios_base_all_saver&);
    };

  } // namespace PDELab
} // namespace Dune
#endif // DUNE_PDELAB_COMMON_IOSTATESAVER_HH
