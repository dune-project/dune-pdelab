#include "config.h"

#include <cassert>
#include <string>
#include <string_view>
#include <vector>

#include <dune/pdelab/logging/consolesink.hh>
#include <dune/pdelab/logging/logger.hh>

namespace Dune::PDELab {

  namespace SinkMessageItems {

    template<typename Factory>
    struct Item
    {

      Item(Factory f)
        : factory(f)
      {}

      Factory factory;
    };

    auto localTime(const LogMessage& msg)
    {
      return Item([&]()
        {
          return msg.localTime();
        });
    }

    using namespace std::literals;
    constexpr auto whitespace = "                                        "sv;

    template<typename Buffer>
    auto logger(const LogMessage& msg, Buffer& buffer, std::size_t width)
    {
      return Item([&,width]() -> std::string_view
        {
          auto logger = name(msg.level());
          buffer.clear();
          buffer.append(logger.data(),logger.data() + logger.size());
          assert(width - logger.size() < 40);
          buffer.append(whitespace.data(),whitespace.data() + width - logger.size());
          return {buffer.data(),buffer.size()};
        });
    }

  }

}

namespace fmt {

  template <typename Factory>
  struct formatter<Dune::PDELab::SinkMessageItems::Item<Factory>>
    : public formatter<std::decay_t<decltype(std::declval<Factory>()())>>
  {

    using Item = Dune::PDELab::SinkMessageItems::Item<Factory>;
    using Base = formatter<std::decay_t<decltype(std::declval<Factory>()())>>;

    // we just inherit the parsing from the original formatter

    template <typename FormatContext>
    auto format(const Item &item, FormatContext &ctx) {
      return Base::format(item.factory(),ctx);
    }

  };

}

namespace Dune::PDELab {

  void ConsoleSink::process(const LogMessage& msg)
  {
    fmt::print(_stream,_format,msg.payload(),SinkMessageItems::logger(msg,_logger_buffer,_widest_logger));
  }

} // end namespace Dune::PDELab
