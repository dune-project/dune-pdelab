// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_LOGTAG_HH
#define DUNE_PDELAB_COMMON_LOGTAG_HH

#include <memory>
#include <ostream>
#include <sstream>
#include <string>

namespace Dune {
  namespace PDELab {

    //! function that writes a log tag to some stream
    /**
     * Usage:
     * \code
#include <dune/pdelab/common/logtag.hh>
int main() {
  using Dune::PDELab::logtag;
  std::cout << logtag << "Hello world!" << std::endl;
}
     * \endcode
     *
     * By default hostPidWallUserLogtagFormatFunc() is used to generate the
     * tags.
     */
    extern std::ostream &logtag(std::ostream &s);

    //! \brief logtag format function that includes host name, pid, wall time
    //!        and CPU time
    extern std::ostream &hostPidWallUserLogtagFormatFunc(std::ostream &s);

    //! \brief logtag format function that includes hostname, rank (if
    //!        available), wall time and CPU time
    /**
     * For serial programs, the rank is omitted completely.  For parallel
     * programs before logtagSetupMPI() has been called, the rank will show as
     * '?'.
     */
    extern std::ostream &hostRankWallUserLogtagFormatFunc(std::ostream &s);

    //! logtag format function that does not write anything
    extern std::ostream &nullFormatFunc(std::ostream &s);

    //! collect MPI information for the logtag formatters
    /**
     * This function needs to be called explicitly after MPI has been set up.
     * Logtag formatters will work before this function has been called but
     * MPI information will be omitted.
     *
     * Besides determining MPI information, this function also syncronises the
     * widthes of the generated output, such that e.g. host names, ranks and
     * pids will print with the same width on each rank.  This can be switched
     * off by passing \c syncWidthes=false.
     */
    extern void logtagSetupMPI(bool syncWidthes = true);

    //! virtual base class for logger formatters
    /**
     * This is the virtual base class used for type erasure purposes.
     */
    struct LogtagFormatterBase {
      //! function that writes the tag to a stream
      /**
       * Implementations of this function can alter the stream format flags as
       * they like.  The flags are restored to their original values by the
       * caller.
       */
      virtual void writeTag(std::ostream &s) const = 0;
    };

    //! A log tag formatter that wraps a unary formatting function or functor
    template<class FormatFunc>
    class GeneralLogtagFormatter :
      public LogtagFormatterBase
    {
      FormatFunc formatFunc;

    public:
      //! constructor
      /**
       * \param formatFunc_ Formatting functor or function-pointer.  Must be
       *                    copy-constructible.
       */
      GeneralLogtagFormatter(const FormatFunc &formatFunc_) :
        formatFunc(formatFunc_)
      { }
      //! write the tag to the stream
      /**
       * This calls formatFunc(s) to write the tag.
       */
      virtual void writeTag(std::ostream &s) const { formatFunc(s); }
    };
    //! Convenience function to create a GeneralLogtagFormatter
    template<class FormatFunc>
    std::shared_ptr<LogtagFormatterBase>
    makeGeneralLogtagFormatter(const FormatFunc &formatFunc)
    { return std::make_shared<GeneralLogtagFormatter<FormatFunc> >(formatFunc); }
    //! Convenience function to create a GeneralLogtagFormatter
    extern std::shared_ptr<LogtagFormatterBase>
    makeGeneralLogtagFormatter(std::ostream &(&formatFunc)(std::ostream&));

    //! get the log tag formatter currently used by logtag()
    extern const std::shared_ptr<LogtagFormatterBase> &getLogtagFormatter();
    //! set a new log tag formatter to be used by logtag()
    /**
     * Calling this with a 0-pointer or no argument restores reinitializes to
     * the formatter in use at startup.
     */
    extern void
    setLogtagFormatter(const std::shared_ptr<LogtagFormatterBase> &formatter
                       = std::shared_ptr<LogtagFormatterBase>());
    //! set a new log tag format function to be used by logtag()
    /**
     * This automatically wraps the function into a GeneralLogtagFormatter
     * before passing it to setLogtagFormatter().
     */
    template<class FormatFunc>
    void setLogtagFormatFunc(const FormatFunc &formatFunc)
    { setLogtagFormatter(makeGeneralLogtagFormatter(formatFunc)); }

    //! temporarily use a different log tag format function
    /**
     * This class sets a different log tag format function and restores the
     * old one when the current scope exits, even if that happens via an
     * exception.
     */
    class WithLogtag {
      std::shared_ptr<LogtagFormatterBase> savedFormatter;

    public:
      template<class FormatFunc>
      WithLogtag(const FormatFunc &formatFunc) :
        savedFormatter(getLogtagFormatter())
      { setLogtagFormatFunc(formatFunc); }

      ~WithLogtag();
    };

    //! Insert standard boilerplate into log messages
    /**
     * This class can be used to create your own logtag locally, e.g. inside a
     * function.  When inserted into a std::ostream, it will insert the
     * standard logtag and some static boilerplate string.  The boilerplate
     * string can be constructed by inserting into the localtag itself.
     * Sample usage:
     * \code
void myfunc() {
  LocalTag mytag;
  mytag << "myfunc(): ";
  std::cout << mytag << "Starting..." << std::endl;

  for(unsigned stage = 0; stage < 42; ++stage) {
    LocalTag stagetag(mytag);
    stagetag << "stage " << stage << ": ";
    std::cout << stagetag << "Starting... << std::endl;

    do_somethin(stage);

    std::cout << stagetag << "Finished." << std::endl;
  }

  std::cout << mytag << "Finished." << std::endl;
}
     * \endcode
     */
    class LocalTag {
      std::string str_;

    public:
      //! extract the static boilerplate message
      inline const std::string &str() const { return str_; }

      //! append something to the static boilerplate message
      template<class V>
      LocalTag &operator<<(const V &v) {
        std::stringstream s;
        s << v;
        str_ += s.str();
        return *this;
      }
    };

    //! insert a localtag into a std::ostream
    inline std::ostream &operator<<(std::ostream &s, const LocalTag &tag)
    { return s << logtag << tag.str(); }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_COMMON_LOGTAG_HH
