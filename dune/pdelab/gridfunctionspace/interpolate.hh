// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_INTERPOLATE_HH
#define DUNE_PDELAB_INTERPOLATE_HH

#include<vector>

#include<dune/common/exceptions.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/typetree/typetree.hh>
#include <dune/typetree/pairtraversal.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    // Backend for standard local interpolation
    struct InterpolateBackendStandard
    {
      template<typename FE, typename ElemFunction, typename XL>
      void interpolate(const FE &fe, const ElemFunction &elemFunction,
                       XL &xl) const
      {
        FiniteElementInterfaceSwitch<FE>::interpolation(fe).
          interpolate(elemFunction,xl);
      }
    };

    namespace {

      template<typename IB, typename LF, typename XG>
      struct InterpolateLeafFromScalarVisitor
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename LFS, typename TreePath>
        void leaf(const LFS& lfs, TreePath treePath) const
        {
          std::vector<typename XG::ElementType> xl(lfs.size());

          // call interpolate for the basis
          ib.interpolate(lfs.finiteElement(), lf, xl);

          // write coefficients into local vector
          xg.write_sub_container(lfs,xl);
        }

        InterpolateLeafFromScalarVisitor(const IB& ib_, const LF& lf_, XG& xg_)
          : ib(ib_)
          , lf(lf_)
          , xg(xg_)
        {}

        const IB& ib;
        const LF& lf;
        XG& xg;

      };


      template<typename IB, typename LF, typename XG>
      struct InterpolateLeafFromVectorVisitor
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename LFS, typename TreePath>
        void leaf(const LFS& lfs, TreePath treePath) const
        {
          std::vector<typename XG::ElementType> xl(lfs.size());

          // call interpolate for the basis

          typedef SelectComponentAdapter<LF> LFCOMP;

          LFCOMP localfcomp(lf,treePath.back());
          ib.interpolate(lfs.finiteElement(), localfcomp, xl);

          // write coefficients into local vector
          xg.write_sub_container(lfs,xl);
        }

        InterpolateLeafFromVectorVisitor(const IB& ib_, const LF& lf_, XG& xg_)
          : ib(ib_)
          , lf(lf_)
          , xg(xg_)
        {}

        const IB& ib;
        const LF& lf;
        XG& xg;

      };


      template<typename IB, typename E, typename XG>
      struct InterpolateVisitor
        : public TypeTree::TreePairVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename F, typename LFS, typename TreePath>
        typename enable_if<F::isLeaf && LFS::isLeaf>::type
        leaf(const F& f, const LFS& lfs, TreePath treePath) const
        {
          std::vector<typename XG::ElementType> xl(lfs.size());

          // call interpolate for the basis
          ib.interpolate(lfs.finiteElement(),
                         GridFunctionToLocalFunctionAdapter<F>(f,e), xl);

          // write coefficients into local vector
          xg.write_sub_container(lfs,xl);
        }

        // interpolate PowerLFS from vector-valued function
        template<typename F, typename LFS, typename TreePath>
        typename enable_if<F::isLeaf && F::Traits::dimRange == 1
                           && (!LFS::isLeaf)>::type
        leaf(const F& f, const LFS& lfs, TreePath treePath) const
        {
          static_assert((TypeTree::TreeInfo<LFS>::depth == 2),
                        "Automatic interpolation of vector-valued function " \
                        "is restricted to trees of depth 1");

          typedef GridFunctionToLocalFunctionAdapter<F> LF;
          LF localf(f,e);

          TypeTree::applyToTree(lfs,InterpolateLeafFromScalarVisitor<IB,LF,XG>(ib,localf,xg));

        }

        // interpolate PowerLFS from vector-valued function
        template<typename F, typename LFS, typename TreePath>
        typename enable_if<F::isLeaf && (F::Traits::dimRange > 1) &&
                           (!LFS::isLeaf)>::type
        leaf(const F& f, const LFS& lfs, TreePath treePath) const
        {
          static_assert((TypeTree::TreeInfo<LFS>::depth == 2),
                        "Automatic interpolation of vector-valued function " \
                        "is restricted to trees of depth 1");
          static_assert(LFS::CHILDREN == F::Traits::dimRange,
                        "Number of children and dimension of range type " \
                        "must match for automatic interpolation of "    \
                        "vector-valued function");

          typedef GridFunctionToLocalFunctionAdapter<F> LF;
          LF localf(f,e);

          TypeTree::applyToTree(lfs,InterpolateLeafFromVectorVisitor<IB,LF,XG>(ib,localf,xg));
        }

        InterpolateVisitor(IB ib_, const E& e_, XG& xg_)
          : ib(ib_)
          , e(e_)
          , xg(xg_)
        {}

      private:
        IB ib;
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

      // make local function space
      typedef LocalFunctionSpace<GFS> LFS;
      LFS lfs(gfs);
      typedef LFSIndexCache<LFS,EmptyTransformation,false> LFSCache;
      LFSCache lfs_cache(lfs);
      typedef typename XG::template LocalView<LFSCache> XView;

      XView x_view(xg);

      // loop once over the grid
      for (const auto& element : elements(entity_set))
        {
          // bind local function space to element
          lfs.bind(element,std::false_type{});
          lfs_cache.update();
          x_view.bind(lfs_cache);

          // call interpolate
          TypeTree::applyToTreePair(f,lfs,InterpolateVisitor<InterpolateBackendStandard,Element,XView>(InterpolateBackendStandard(),element,x_view));

          x_view.unbind();
        }

      x_view.detach();
    }

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
