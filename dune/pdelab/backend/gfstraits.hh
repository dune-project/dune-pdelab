template<typename GFS>
struct GFSTraits
{
    // Multiindex with which we can access entries in a DOF container
    using ContainerIndex = typename GFS::Ordering::Traits::ContainerIndex;

    auto static blockCount(const std::shared_ptr<const GFS> & gfs)
    {
        return gfs->ordering().blockCount();
    }
};
