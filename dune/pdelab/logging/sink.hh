#ifndef DUNE_PDELAB_LOGGING_SINK_HH
#define DUNE_PDELAB_LOGGING_SINK_HH

#include <string>
#include <string_view>

#include <dune/pdelab/logging/logmessage.hh>

namespace Dune::PDELab {

  /**
   * \addtogroup logging
   * \{
   */

  //! A sink is responsible for actually persisting a LogMessage.
  /**
   * Sink is the abstract base class for all sink implementations in the framework.
   */
  class Sink
  {

    friend class Logging;

  public:

    //! Returns the maximum LogLevel for which this Sink is enabled.
    LogLevel level() const
    {
      return _level;
    }

    //! Sets the maximum LogLevel for which this Sink is enabled.
    void setLevel(LogLevel level)
    {
      _level = level;
    }

    //! Returns the name of this Sink.
    const std::string_view name() const
    {
      return _name;
    }

    //! Returns the length of the longest logger name currently known to the system.
    /**
     * This information is automatically updated by the logging system and can be used for creating
     * a logger name log field of uniform width.
     */
    std::size_t widestLogger() const
    {
      return _widest_logger;
    }

  private:

    //! Sets the name of this Sink.
    void setName(std::string_view name)
    {
      _name = name;
    }

    //! Sets the length of the longest logger name currently known to the system.
    void setWidestLogger(std::size_t widest_logger)
    {
      _widest_logger = widest_logger;
    }

  public:

    //! Process the given LogMessage.
    /**
     * \note This method will only be called by the logging system if the LogLevel of the LogMessage
     * is not larger than the one of this sink, so implementations do not have to check for that.
     */
    virtual void process(const LogMessage& msg) = 0;

    //! Virtual destructor as required for classes with virtual methods.
    virtual ~Sink();

  protected:

    //! Initialize Sink base class with the given information.

    Sink(std::string_view name, LogLevel level, std::size_t widest_logger)
      : _level(level)
      , _widest_logger(widest_logger)
      , _name(name)
    {}

  private:

    LogLevel _level = LogLevel::all;
    std::size_t _widest_logger = 0;
    std::string _name;

  };


  //! This sink does nothing.
  class NullSink final
    : public Sink
  {

  public:

    //! Constructs a new NullSink with the given name.
    NullSink(std::string_view name)
      : Sink(name,LogLevel::off,0)
    {}

    //! Empty processing function.
    void process(const LogMessage&) override;

  };

  /**
   * \}
   */

} // end namespace Dune::PDELab

#endif // DUNE_PDELAB_COMMON_LOGGING_SINK_HH
