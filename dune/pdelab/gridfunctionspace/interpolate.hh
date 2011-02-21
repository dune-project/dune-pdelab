// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_INTERPOLATE_HH
#define DUNE_PDELAB_INTERPOLATE_HH

#include<vector>

#include<dune/common/exceptions.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include "../common/typetree.hh"
#include "../common/typetree/pairtraversal.hh"
#include "../common/function.hh"

#include "gridfunctionspace.hh"

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
          lfs.vwrite(xl,xg);
        }

        // interpolate PowerLFS from vector-valued function
        template<typename F, typename LFS, typename TreePath>
        typename enable_if<F::isLeaf && F::Traits::dimRange == 1
                           && (!LFS::isLeaf)>::type
        leaf(const F& f, const LFS& lfs, TreePath treePath) const
        {
#warning TODO: Implement interpolate from scalar function onto composite GFS
          dune_static_assert(LFS::isPower,
                             "Automatic interpolation of vector-valued function " \
                             "only works for PowerGridFunctionSpace");
          dune_static_assert((LFS::template Child<0>::Type::isLeaf),
                             "Automatic interpolation of vector-valued function " \
                             "is restricted to trees of depth 1");

          typedef GridFunctionToLocalFunctionAdapter<F> LF;
          LF localf(f,e);

          for (std::size_t k=0; k<LFS::CHILDREN; ++k)
            {
              // allocate vector where to store coefficients from basis
              std::vector<typename XG::ElementType> xl(lfs.child(k).size());

              // call interpolate for the basis
              ib.interpolate(lfs.child(k).finiteElement(), lf, xl);

              // write coefficients into local vector
              lfs.child(k).vwrite(xl,xg);
            }
        }

        // interpolate PowerLFS from vector-valued function
        template<typename F, typename LFS, typename TreePath>
        typename enable_if<F::isLeaf && (F::Traits::dimRange > 1) &&
                           (!LFS::isLeaf)>::type
        leaf(const F& f, const LFS& lfs, TreePath treePath) const
        {
          dune_static_assert(LFS::isPower,
                             "Automatic interpolation of vector-valued function " \
                             "only works for PowerGridFunctionSpace");
          dune_static_assert((LFS::template Child<0>::Type::isLeaf),
                             "Automatic interpolation of vector-valued function " \
                             "is restricted to trees of depth 1");
          dune_static_assert(LFS::CHILDREN == F::Traits::dimRange,
                             "Number of children and dimension of range type " \
                             "must match for automatic interpolation of " \
                             "vector-valued function");
          for (std::size_t k=0; k<LFS::CHILDREN; ++k)
            {
              // allocate vector where to store coefficients from basis
              std::vector<typename XG::ElementType> xl(lfs.child(k).size());

              // call interpolate for the basis
              typedef GridFunctionToLocalFunctionAdapter<F> LF;
              LF localf(f,e);
              typedef SelectComponentAdapter<LF> LFCOMP;
              LFCOMP localfcomp(localf,k);
              ib.interpolate(lfs.child(k).finiteElement(), localfcomp, xl);

              // write coefficients into local vector
              lfs.child(k).vwrite(xl,xg);
            }
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
      typedef typename GFS::Traits::GridViewType GV;
      typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
      typedef typename GV::Traits::template Codim<0>::Entity Element;

      // make local function space
      typedef LocalFunctionSpace<GFS> LFS;
      LFS lfs(gfs);

      // loop once over the grid
      for (ElementIterator it = gfs.gridview().template begin<0>();
           it!=gfs.gridview().template end<0>(); ++it)
        {
          // bind local function space to element
          lfs.bind(*it);

          // call interpolate
          TypeTree::applyToTreePair(f,lfs,InterpolateVisitor<InterpolateBackendStandard,Element,XG>(InterpolateBackendStandard(),*it,xg));
        }
    }

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
