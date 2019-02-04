#ifndef DUNE_PDELAB_LOGGING_SINKMESSAGEITEMS_HH
#define DUNE_PDELAB_LOGGING_SINKMESSAGEITEMS_HH

#include <cassert>
#include <string>
#include <string_view>
#include <vector>

#include <dune/pdelab/logging/fmt.hh>
#include <dune/pdelab/logging/logger.hh>

namespace Dune::PDELab {

  namespace SinkMessageItems {

    inline auto localTime(const LogMessage& msg)
    {
      return LazyFormatArgument([&]()
        {
          return msg.localTime();
        });
    }

    inline auto relativeDays(const LogMessage& msg)
    {
      return LazyFormatArgument([&]()
        {
          return msg.relativeDays();
        });
    }

    template<typename Buffer>
    inline auto backend(const LogMessage& msg, Buffer& buffer, std::size_t width)
    {
      return LazyFormatArgument([&,width]() -> std::string_view
        {
          auto logger = msg.logger().name();
          auto size = logger.size();
          if (size < width)
          {
            buffer.clear();
            buffer.append(logger.data(),logger.data() + size);
            buffer.resize(width);
            std::fill_n(buffer.data() + size,width - size,' ');
            return {buffer.data(),buffer.size()};
          }
          else
          {
            return logger;
          }
        });
    }

  } // end namespace SinkMessageItems

} // end namespace Dune::PDELab


#endif // DUNE_PDELAB_LOGGING_SINKMESSAGEITEMS_HH
