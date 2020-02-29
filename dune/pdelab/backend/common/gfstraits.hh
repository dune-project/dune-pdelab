namespace Dune::PDELab {

template<typename GFS, typename Enable = void>
struct GFSTraits;

template<typename GFS>
struct GFSTraits<GFS,
                 std::enable_if_t<isGridFunctionSpace<GFS>()>>
{
    // Multiindex with which we can access entries in a DOF container
    using ContainerIndex = typename GFS::Ordering::Traits::ContainerIndex;

    auto static blockCount(const std::shared_ptr<const GFS> & gfs)
    {
        return gfs->ordering().blockCount();
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

    auto static blockCount(const std::shared_ptr<const FSB> & gfs)
    {
        return 1;
    }
};

}
