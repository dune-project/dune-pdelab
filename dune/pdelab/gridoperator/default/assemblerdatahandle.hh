#ifndef ASSEMBLERDATAHANDLE_HH
#define ASSEMBLERDATAHANDLE_HH

#include <dune/grid/common/datahandleif.hh>
#include <dune/pdelab/common/polymorphicbufferwrapper.hh>

namespace Dune{
  namespace PDELab{

    template <typename GFSU, typename GFSV, typename CU, typename CV, typename LocalAssemblerEngine>
    class AssemblerDataHandle
    {

      using EntitySet = typename GFSU::Traits::EntitySet;
      using Element = typename EntitySet::Element;
      using Intersection = typename EntitySet::Intersection;

      // local function spaces
      using LFSU = LocalFunctionSpace<GFSU, TrialSpaceTag>;
      using LFSV =  LocalFunctionSpace<GFSV, TestSpaceTag>;
      using LFSUCache = LFSIndexCache<LFSU,CU>;
      using LFSVCache = LFSIndexCache<LFSV,CV>;


      // Gather and scatter method only differ by function call in the
      // assembler engine. We use this function to avoid code duplication.
      template <typename MessageBuffer, typename EntityType>
      void gatherScatter(MessageBuffer& buff, const EntityType& e, bool gather) const
      {
        using Buf = PolymorphicBufferWrapper<MessageBuffer>;
        Buf bufWrapper(buff);

        Element element = e.inside();
        ElementGeometry<Element> eg(element);

        // Bind local test function space to element
        lfsv_.bind(element);
        lfsvCache_.update();

        // Notify assembler engine about bind
        assemblerEngine_.onBindLFSV(eg,lfsvCache_);

        // Bind local trial function space to element
        lfsu_.bind(element);
        lfsuCache_.update();

        // Notify assembler engine about bind
        assemblerEngine_.onBindLFSUV(eg,lfsuCache_,lfsvCache_);

        // Load coefficients of local functions
        assemblerEngine_.loadCoefficientsLFSUInside(lfsuCache_);

        // call gather or scatter method of local operator
        if (gather){
          assemblerEngine_.assembleUVProcessBoundaryGather(e,lfsuCache_,lfsvCache_,bufWrapper);
        }
        else{
          assemblerEngine_.assembleUVProcessBoundaryScatter(e,lfsuCache_,lfsvCache_,bufWrapper);
        }

        // Notify assembler engine about unbinds
        assemblerEngine_.onUnbindLFSUV(eg,lfsuCache_,lfsvCache_);

        // Notify assembler engine about unbinds
        assemblerEngine_.onUnbindLFSV(eg,lfsvCache_);
      }


    public:
      using DataType = char;

      AssemblerDataHandle(const GFSU& gfsu, const GFSV& gfsv,
                          LocalAssemblerEngine& assemblerEngine)
        : cdim_(1)
        , assemblerEngine_(assemblerEngine)
        , lfsu_(gfsu)
        , lfsv_(gfsv)
        , lfsuCache_(lfsu_)
        , lfsvCache_(lfsv_)
        , doAlpha_(doAlpha)
        , doLambda_(doLambda)
      {}

      //! returns true if data for this codim should be communicated
      bool contains (int dim, int codim) const
      {
        return (codim==cdim_);
      }

      //! returns true if size per entity of given dim and codim is a constant
      bool fixedsize (int dim, int codim) const
      {
        return assemblerEngine_.communicationFixedSize();
      }

      /*! how many objects of type DataType have to be sent for a given entity
        Note: Only the sender side needs to know this size.
      */
      template<class EntityType>
      size_t size (EntityType& e) const
      {
        return assemblerEngine_.communicationSize(e,lfsuCache_,lfsvCache_);
      }


      //! pack data from user to message buffer
      template<class MessageBuffer, class EntityType>
      void gather (MessageBuffer& buff, const EntityType& e) const
      {
        gatherScatter(buff,e,true);
      }

      // //! unpack data from message buffer to user
      template<class MessageBuffer, class EntityType>
      void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
      {
        gatherScatter(buff,e,false);
      }

    private:
      int cdim_;

      LocalAssemblerEngine& assemblerEngine_;

      // local function spaces in local cell
      mutable LFSU lfsu_;
      mutable LFSV lfsv_;

      // corresponding caches
      mutable LFSUCache lfsuCache_;
      mutable LFSVCache lfsvCache_;
    };

  } // namespace PDELab
} // namespace Dune

#endif
