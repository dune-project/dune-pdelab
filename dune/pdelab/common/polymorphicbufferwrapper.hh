// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMMON_POLYMORPHICBUFFERWRAPPER_HH
#define DUNE_PDELAB_COMMON_POLYMORPHICBUFFERWRAPPER_HH

#include <cstddef>
#include <cassert>

namespace Dune {
  namespace PDELab {

    //! Wrapper for message buffers of grid DataHandles that allows for sending different types of data.
    /**
     * The standard message buffers passed to the callbacks of the grid DataHandles
     * are templated on a specific data type and do not support writing data of other
     * types.
     *
     * This wrapper takes a MessageBuffer for char and allows you to write any kind of POD data
     * into it. It works by simply interpreting the data as a byte stream and serializing / deserializing
     * it. Be aware that this implementation may create problems on heterogeneous architectures with
     * different byte orderings.
     *
     * Moreover, the buffer can optionally take care of providing information about the MPI rank of
     * the sender through the method senderRank() if the rank of the current process and
     * transmit_rank = true are passed to the constructor of the buffer.
     *
     * \warning The underlying MessageBuffer *must* use char as its data type, otherwise you will probably
     *          get compile or runtime errors!
     */
    template<typename Buffer>
    class PolymorphicBufferWrapper
    {

    public:

      enum class Mode
      {
        send,
        receive,
      };

      template<typename T>
      void write(const T& data)
      {
        assert(_mode == Mode::send);
        const char* raw_data = reinterpret_cast<const char*>(&data);
        for (std::size_t i = 0; i < sizeof(T); ++i)
          _buffer.write(*(raw_data++));
      }

      template<typename T>
      void read(T& data)
      {
        assert(_mode == Mode::receive);
        char* raw_data = reinterpret_cast<char*>(&data);
        for (std::size_t i = 0; i < sizeof(T); ++i)
          _buffer.read(*(raw_data++));
      }

      PolymorphicBufferWrapper(Buffer& buffer, Mode mode, int rank = -1, bool transmit_rank = false)
        : _buffer(buffer)
        , _mode(mode)
        , _rank(rank)
        , _transmit_rank(transmit_rank)
        , _sender_rank(-1)
      {
        if (_transmit_rank)
          {
            assert(_rank >= 0);
            if (_mode == Mode::receive)
              read(_sender_rank);
            else
              write(_rank);
          }
      }

      Mode mode() const
      {
        return _mode;
      }

      bool transmitRank() const
      {
        return _transmit_rank;
      }

      int senderRank() const
      {
        assert(_sender_rank >= 0);
        return _sender_rank;
      }

    private:

      Buffer& _buffer;
      const Mode _mode;
      int _rank;
      const bool _transmit_rank;
      int _sender_rank;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_COMMON_POLYMORPHICBUFFERWRAPPER_HH
