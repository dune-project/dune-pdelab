// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_FUNCTION_SELECTCOMPONENT_HH
#define DUNE_PDELAB_FUNCTION_SELECTCOMPONENT_HH

#include <algorithm>
#include <cstddef>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/pdelab/common/function.hh>

namespace Dune {
  namespace PDELab {

    //! Select certain component(s) of a gridfunction
    /**
     * \tparam GF   Type for the GridFunction
     * \tparam dimR Number of components of the resulting GridFunction
     *
     * \sa SelectComponentAdaptor, BoundaryGridFunctionSelectComponentAdaptor.
     */
    template<typename GF, std::size_t dimR = 1>
    class SelectComponentGridFunctionAdapter
      : public GridFunctionBase<
          GridFunctionTraits<
            typename GF::Traits::GridViewType,
            typename GF::Traits::RangeFieldType, dimR,
            FieldVector<typename GF::Traits::RangeFieldType, dimR>
            >,
          SelectComponentGridFunctionAdapter<GF, dimR>
        >
    {
    public:
      typedef GridFunctionTraits<
        typename GF::Traits::GridViewType,
        typename GF::Traits::RangeFieldType, dimR,
        FieldVector<typename GF::Traits::RangeFieldType, dimR>
        > Traits;

    private:
      typedef GridFunctionBase<
        Traits, SelectComponentGridFunctionAdapter
        > Base;

      GF& gf;
      std::size_t remap[dimR];

      void checkRemap() const {
        for(std::size_t c = 0; c < dimR; ++c)
          if(remap[c] >= GF::Traits::dimRange)
            DUNE_THROW(InvalidStateException, "remap[" << c << "] = "
                       << remap[c] << " >= GF::Traits::dimRange = "
                       << GF::Traits::dimRange);
      }
    public:
      //! construct with a consecutive range of indices
      SelectComponentGridFunctionAdapter(GF& gf_, std::size_t first) :
        gf(gf_)
      {
        for(std::size_t c = 0; c < dimR; ++c)
          remap[c] = first+c;
        checkRemap();
      }

      //! construct with a full index map
      SelectComponentGridFunctionAdapter
      ( GF& gf_, const std::vector<std::size_t> remap_) :
        gf(gf_)
      {
        if(remap_.size() != dimR)
          DUNE_THROW(Exception, "Got an index map of size "
                     << remap_.size() << " but size " << dimR << " was "
                     "expected.");
        std::copy(remap_.begin(), remap_.end(), remap);
        checkRemap();
      }

      void evaluate(const typename Traits::ElementType &e,
                    const typename Traits::DomainType &x,
                    typename Traits::RangeType &y) const
      {
        typename GF::Traits::RangeType y_;
        gf.evaluate(e,x,y_);
        for(std::size_t c = 0; c < dimR; ++c)
          y[c] = y_[remap[c]];
      }

      const typename Traits::GridViewType& getGridView() const
      { return gf.getGridView(); }

      template<typename Time>
      void setTime(Time time) { gf.setTime(time); }
    };

  } // namspace PDELab
} // namspace Dune

#endif // DUNE_PDELAB_FUNCTION_SELECTCOMPONENT_HH
