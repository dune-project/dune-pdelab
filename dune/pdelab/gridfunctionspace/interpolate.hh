// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_INTERPOLATE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_INTERPOLATE_HH

#include <vector>
#include <utility>

#include <dune/common/exceptions.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/localfunctions/common/virtualinterface.hh>

#include <dune/typetree/typetree.hh>
#include <dune/typetree/pairtraversal.hh>

#include <dune/pdelab/function/localfunction.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{
    namespace {

      template<typename LF, typename XG>
      struct InterpolateLeafFromScalarVisitor
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename LFS, typename TreePath>
        void leaf(const LFS& lfs, TreePath treePath) const
        {
          std::vector<typename XG::ElementType> xl(lfs.size());

          // call interpolate for the basis
          FiniteElementInterfaceSwitch<typename LFS::Traits::FiniteElement>::interpolation(lfs.finiteElement()).
            interpolate(lf, xl);

          // write coefficients into local vector
          xg.write_sub_container(lfs,xl);
        }

        InterpolateLeafFromScalarVisitor(const LF& lf_, XG& xg_)
          : lf(lf_)
          , xg(xg_)
        {}

        const LF& lf;
        XG& xg;
      };


      template<typename LF, typename XG>
      struct InterpolateLeafFromVectorVisitor
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {
        using Domain = typename Functions::SignatureTraits<LF>::Domain;
        using Range = typename Functions::SignatureTraits<LF>::Range;
        using RangeField = typename FieldTraits<Range>::field_type;

        template<typename LFS, typename TreePath>
        void leaf(const LFS& lfs, TreePath treePath) const
        {
          std::vector<typename XG::ElementType> xl(lfs.size());

          using LFSRange = typename LFS::Traits::FiniteElement::Traits::LocalBasisType::Traits::RangeType;
          static_assert(std::is_convertible<LFSRange, typename FieldTraits< LFSRange >::field_type>::value,
            "only interpolation into scalar leaf function spaces is implemented");

          // call interpolate for the basis
          auto f = [&](const Domain& x) -> RangeField { return lf(x)[index]; };

          FiniteElementInterfaceSwitch<typename LFS::Traits::FiniteElement>::interpolation(lfs.finiteElement()).
            interpolate(f, xl);

          // write coefficients into local vector
          xg.write_sub_container(lfs,xl);

          // increment index
          index++;
        }

        InterpolateLeafFromVectorVisitor(const LF& lf_, XG& xg_)
          : lf(lf_)
          , index(0)
          , xg(xg_)
        {}

        const LF& lf;
        mutable std::size_t index;
        XG& xg;
      };


      template<typename E, typename XG>
      struct InterpolateVisitor
        : public TypeTree::TreePairVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename F, typename LFS, typename TreePath>
        typename std::enable_if<F::isLeaf && LFS::isLeaf>::type
        leaf(const F& f, const LFS& lfs, TreePath treePath) const
        {
          std::vector<typename XG::ElementType> xl(lfs.size());
           // call interpolate for the basis
          FiniteElementInterfaceSwitch<typename LFS::Traits::FiniteElement>::interpolation(lfs.finiteElement()).
            interpolate(f, xl);

          // write coefficients into local vector
          xg.write_sub_container(lfs,xl);
        }

        // interpolate PowerLFS from scalar-valued function
        template<typename F, typename LFS, typename TreePath,
                 typename Range = typename Functions::SignatureTraits<F>::Range>
        typename std::enable_if<F::isLeaf &&
                           std::is_convertible<Range, typename FieldTraits< Range >::field_type>::value &&
                           (!LFS::isLeaf)>::type
        leaf(const F& f, const LFS& lfs, TreePath treePath) const
        {
          // call interpolate for the basis
          TypeTree::applyToTree(lfs,InterpolateLeafFromScalarVisitor<F,XG>(f, xg));
        }

        // interpolate PowerLFS from vector-valued function
        template<typename F, typename LFS, typename TreePath,
                 typename Range = typename Functions::SignatureTraits<F>::Range>
        typename std::enable_if<F::isLeaf &&
                          (!std::is_convertible<Range, typename FieldTraits< Range >::field_type>::value) &&
                          (!LFS::isLeaf)>::type
        leaf(const F& f, const LFS& lfs, TreePath treePath) const
        {
          static_assert(TypeTree::TreeInfo<LFS>::leafCount == Range::dimension,
                        "Number of leaves and dimension of range type " \
                        "must match for automatic interpolation of "    \
                        "vector-valued function");

          TypeTree::applyToTree(lfs,InterpolateLeafFromVectorVisitor<F,XG>(f,xg));
        }

        InterpolateVisitor(const E& e_, XG& xg_)
          : e(e_)
          , xg(xg_)
        {}

      private:
        const E& e;
        XG& xg;
      };

    } // anonymous namespace

    //! interpolation from a given grid function
    /**
     * \code
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
     * \endcode
     * \param f   Function to interpolate from.
     * \param gfs GridFunctionSpace to use for interpoaltion.
     * \param xg  Global vector of dofs to interpolate into.
     *
     * \note \c xg needs to be initialized to the correct size, but there is
     *       no need to initialize its contents.
     */
    template<typename F, typename GFS, typename XG>
    void interpolate (const F& f, const GFS& gfs, XG& xg)
    {
      // this is the leaf version now

      // get some types
      using EntitySet = typename GFS::Traits::EntitySet;
      using Element = typename EntitySet::Element;

      auto entity_set = gfs.entitySet();

      // make local function
      auto lf = makeLocalFunctionTree(f, gfs.gridView());

      // make local function space
      typedef LocalFunctionSpace<GFS> LFS;
      LFS lfs(gfs);
      typedef LFSIndexCache<LFS> LFSCache;
      LFSCache lfs_cache(lfs);
      typedef typename XG::template LocalView<LFSCache> XView;

      XView x_view(xg);

      // loop once over the grid
      for (const auto& element : elements(entity_set))
        {
          // bind local function space to element
          lf.bind(element);
          lfs.bind(element);
          lfs_cache.update();
          x_view.bind(lfs_cache);

          // call interpolate
          TypeTree::applyToTreePair(lf,lfs,InterpolateVisitor<Element,XView>(element,x_view));

          x_view.unbind();
          lf.unbind();
        }

      x_view.detach();
    }

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_INTERPOLATE_HH
