#ifndef DUNE_PDELAB_BACKEND_COMMON_GFS_TRAITS_HH
#define DUNE_PDELAB_BACKEND_COMMON_GFS_TRAITS_HH

#include <dune/functions/functionspacebases/sizeinfo.hh>

namespace Dune::PDELab {

template<class GFS>
class GFSSizeInfo
{
public:
    using size_type = std::size_t;
    using SizePrefix = Dune::ReservedVector<std::size_t,
                                            GFS::Ordering::Traits::ContainerIndex::capacity()>;

    /**
     * \brief Construct from basis
     */
    GFSSizeInfo(const GFS& gfs) : ordering_(gfs.ordering())
    {}

    /**
     * \brief Return number possible values for next position in multi index
     */
    size_type operator()(const SizePrefix& prefix) const
    {
        return size(prefix);
    }

    /**
     * \brief Return number possible values for next position in multi index
     *
     * This shall vanish. It's just here such that this can be used
     * as size provider n place of the basis.
     */
    size_type size(const SizePrefix& prefix) const
    {
        return childOrderingSize(ordering_, prefix, 0);
    }

    operator size_type () const
    {
        return ordering_.blockCount();
    }

private:
    const typename GFS::Ordering & ordering_;

    template<typename Ordering>
    static size_type childOrderingSize(const Ordering & ordering, const SizePrefix& prefix, int pos)
    {
        if (pos < prefix.size())
        {
            if(ordering.childOrderingCount())
                return childOrderingSize(ordering.childOrdering(prefix[pos]), prefix, pos+1);
            else // this should only happen if we are the leaf...
                return 0; // indicates static size
                // return ordering.size()/ordering.blockCount();
        }
        return ordering.blockCount();
    }

};

template<typename GFS, typename Enable = void>
struct GFSTraits;

template<typename GFS>
struct GFSTraits<GFS,
                 std::enable_if_t<isGridFunctionSpace<GFS>()>>
{
    // Multiindex with which we can access entries in a DOF container
    using ContainerIndex = typename GFS::Ordering::Traits::ContainerIndex;

    auto static sizeInfo(const std::shared_ptr<const GFS> & gfs)
    {
        return GFSSizeInfo<GFS>(*gfs);
    }
};

template<typename FSB>
struct GFSTraits<FSB,
                 std::enable_if_t<isBasisInfo<FSB>()>>
{
    // Multiindex with which we can access entries in a DOF container
    using Basis = typename FSB::Basis;
    using LocalView = typename Basis::LocalView;
    using ContainerIndex = typename LocalView::MultiIndex;

    auto static sizeInfo(const std::shared_ptr<const FSB> & gfs)
    {
        return Dune::Functions::sizeInfo(gfs->basis());
    }
};

}

#endif // DUNE_PDELAB_BACKEND_COMMON_GFS_TRAITS_HH
