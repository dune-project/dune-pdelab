// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_COMMON_POLYMORPHICBUFFERWRAPPER_HH
#define DUNE_PDELAB_COMMON_POLYMORPHICBUFFERWRAPPER_HH

#include <cstddef>

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
     * \warning The underlying MessageBuffer *must* use char as its data type, otherwise you will probably
     *          get compile or runtime errors!
     */
    template<typename Buffer>
    class PolymorphicBufferWrapper
    {

    public:

      template<typename T>
      void write(const T& data)
      {
        const char* raw_data = reinterpret_cast<const char*>(&data);
        for (std::size_t i = 0; i < sizeof(T); ++i)
          _buffer.write(*(raw_data++));
      }

      template<typename T>
      void read(T& data)
      {
        char* raw_data = reinterpret_cast<char*>(&data);
        for (std::size_t i = 0; i < sizeof(T); ++i)
          _buffer.read(*(raw_data++));
      }

      PolymorphicBufferWrapper(Buffer& buffer)
        : _buffer(buffer)
      {}

    private:

      Buffer& _buffer;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_COMMON_POLYMORPHICBUFFERWRAPPER_HH
