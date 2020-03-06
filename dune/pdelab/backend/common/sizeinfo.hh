#ifndef DUNE_PDELAB_BACKEND_ISTL_SIZEINFO_HH
#define DUNE_PDELAB_BACKEND_ISTL_SIZEINFO_HH

#include <dune/common/reservedvector.hh>
#include <dune/functions/functionspacebases/sizeinfo.hh>
#include <dune/pdelab/common/concepts.hh>
#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/istl/tags.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>

namespace Dune {
namespace PDELab {
namespace Backend {

/** implementation of the Dune::Functions::SizeInfo interface for a GridFunctionSpace */
template<class GFS>
class GFSSizeInfo
{

    using block_vector = ISTL::tags::block_vector;
    using field_vector = ISTL::tags::field_vector;
    template<typename C>
    using container_t = ISTL::tags::container_t<C>;

    // ********************************************************************************
    // collect size information
    // ********************************************************************************
    static void setSize(block_vector, std::vector<std::size_t>& sizes, std::size_t offset, std::size_t size)
    {
        sizes[offset] = size;
    }

    static void setSize(field_vector, std::vector<std::size_t>&, std::size_t, std::size_t)
    {}

    template<typename V, typename Ordering>
    static void getSizes(field_vector, const Ordering& ordering, std::vector<std::size_t>& sizes, std::size_t offset)
    {}

    template<typename V, typename Ordering>
    static void getSizes(block_vector, const Ordering& ordering, std::vector<std::size_t>& sizes, std::size_t offset)
    {
        for (std::size_t i = 0; i < ordering.childOrderingCount(); ++i)
        {
            if (ordering.containerBlocked())
            {
                sizes.push_back(0);
                // dispatch on the block_type
                using block_type = typename V::block_type;
                using childtag = container_t<block_type>;
                setSize(childtag(),sizes,offset+1,ordering.childOrdering(i).blockCount());
                getSizes<block_type>(childtag(), ordering.childOrdering(i), sizes, offset+1);
            }
            else
                getSizes<V>(block_vector(), ordering.childOrdering(i), sizes, offset);
        }
    }

    template<typename V, typename Ordering>
    static void dispatchGetSizes(const Ordering& ordering, std::vector<std::size_t>& sizes,
        HierarchicContainerAllocationTag)
    {
        using tag = container_t<V>;
        getSizes<V>(tag(), ordering, sizes, 0);
    }

    template<typename V, typename Ordering>
    static void dispatchGetSizes(const Ordering& ordering, std::vector<std::size_t>& sizes,
        FlatContainerAllocationTag)
    {}

    static std::vector<std::size_t> getSizes(const GFS& gfs)
    {
        using Container = typename Backend::Vector<GFS,double>::Container;

        std::vector<std::size_t> sizes({gfs.ordering().blockCount()});
        dispatchGetSizes<Container>(gfs.ordering(), sizes, typename GFS::Ordering::ContainerAllocationTag());
        return sizes;
    }


public:
    using size_type = std::size_t;
    using SizePrefix = Dune::ReservedVector<std::size_t,
                                            GFS::Ordering::Traits::ContainerIndex::capacity()>;

    /**
     * \brief Construct from basis
     */
    GFSSizeInfo(const GFS& gfs) :
        //ordering_(gfs.ordering())
        sizes_(getSizes(gfs))
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
        if (prefix.size() >= sizes_.size())
            return 0;
        return sizes_[prefix.size()];
    }

    operator size_type () const
    {
        return sizes_[0];
    }

private:
    std::vector<std::size_t> sizes_;
};

template<typename GFS,
         std::enable_if_t<isGridFunctionSpace<GFS>(),int> = 0>
auto sizeInfo(const GFS & gfs) -> GFSSizeInfo<GFS> 
{
    return GFSSizeInfo<GFS>(gfs);
}

template<typename B,
         std::enable_if_t<models<Concept::BasisInfo, B>(),int> = 0>
auto sizeInfo(const B & basis) -> decltype(Dune::Functions::sizeInfo(basis.basis()))
{
    return Dune::Functions::sizeInfo(basis.basis());
}

}}} // Dune::PDELab::Backend

#endif // DUNE_PDELAB_BACKEND_ISTL_SIZEINFO_HH
