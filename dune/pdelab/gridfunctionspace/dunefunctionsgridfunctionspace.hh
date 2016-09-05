// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_DUNEFUNCTIONSGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_DUNEFUNCTIONSGRIDFUNCTIONSPACE_HH

#include <cstddef>
#include <map>
#include <ostream>
#include <set>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/dunefunctionslocalfunctionspace.hh>

#include <dune/pdelab/backend/istl/dunefunctions.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

    namespace Experimental {

      /** \brief A pdelab grid function space implemented by a dune-functions function space basis
       *
       *  \tparam DFBasis A dune-functions function space basis
       *  \tparam V       The type of the underlying ISTL vector
       *  \tparam CE      Type for constraints assembler
       */
      template<typename DFBasis, typename V, typename CE>
      class GridFunctionSpace
        : public TypeTree::LeafNode
        , public GridFunctionOutputParameters
        //        , public DataHandleProvider<DuneFunctionsGridFunctionSpace<DFBasis,CE,B,P> >
      {
        using GV = typename DFBasis::GridView;

        template<typename,typename>
        friend class GridFunctionSpaceBase;

      public:
        //! export Traits class

        struct Traits {

          using GridView  = Dune::PDELab::impl::GridView<typename DFBasis::GridView>;
          using EntitySet = Dune::PDELab::impl::EntitySet<typename DFBasis::GridView>;

          using size_type = std::size_t;
          using ConstraintsType = CE;

          using Basis = DFBasis;

          using Backend = istl::SimpleVectorBackend<V>;

          struct FEM
          {
            struct Traits
            {
              using FiniteElement = typename DFBasis::LocalView::Tree::FiniteElement;
              using FiniteElementType = FiniteElement;
            };
          };

          using FiniteElementMap = FEM;
          using FiniteElementMapType = FEM;

        };

        using Basis          = DFBasis;

        struct Ordering {

          struct Traits {

            using DOFIndex       = typename DFBasis::MultiIndex;
            using ContainerIndex = DOFIndex;
            using SizeType       = std::size_t;

          };

          using DOFIndex       = typename DFBasis::MultiIndex;
          using ContainerIndex = DOFIndex;
          using size_type      = std::size_t;

          using CacheTag       = DuneFunctionsCacheTag;
          using ContainerAllocationTag = FlatContainerAllocationTag;

          Ordering(const GridFunctionSpace& gfs)
            : _gfs(gfs)
          {}

          size_type size() const
          {
            return _gfs.basis().size();
          }

          size_type blockCount() const
          {
            return size() / V::block_type::dimension;
          }

          size_type maxLocalSize() const
          {
            return _gfs.basis().maxLocalSize();
          }

          ContainerIndex mapIndex(const DOFIndex& di) const
          {
            return di;
          }

          void mapIndex(const DOFIndex& di, ContainerIndex& ci) const
          {
            ci = di;
          }

          void update()
          {}

        private:

          const GridFunctionSpace& _gfs;

        };


        //! extract type for storing constraints
        template<typename E>
        struct ConstraintsContainer
        {

          //! \brief define Type as the Type of a container of E's
          using Type = std::conditional_t<
            std::is_same<
              CE,
              NoConstraints
              >::value,
            EmptyTransformation,
            ConstraintsTransformation<typename Ordering::Traits::DOFIndex,typename Ordering::Traits::ContainerIndex,E>
            >;

        private:
          ConstraintsContainer () {}
        };

        // ****************************************************************************************************
        // Construct from a dune-functions basis
        // ****************************************************************************************************

        //! constructor
        GridFunctionSpace (std::shared_ptr<DFBasis> df_basis, std::shared_ptr<CE> ce)
          : _es(df_basis->gridView())
          , _df_basis(std::move(df_basis))
          , _pce(std::move(ce))
          , _ordering(*this)
        {}

        GridFunctionSpace (std::shared_ptr<DFBasis> df_basis)
          : _es(df_basis->gridView())
          , _df_basis(std::move(df_basis))
          , _pce(std::make_shared<CE>())
          , _ordering(*this)
        {}

        //! get grid view
        const typename Traits::GridView& gridView () const
        {
          return _es.gridView();
        }

        //! get EntitySet
        const typename Traits::EntitySet& entitySet () const
        {
          return _es;
        }

        //! return constraints engine
        const typename Traits::ConstraintsType& constraints () const
        {
          return *_pce;
        }

        //! return storage of constraints engine
        std::shared_ptr<const CE> constraintsStorage () const
        {
          return _pce;
        }

        //! Direct access to the DOF ordering.
        const Ordering& ordering() const
        {
          return _ordering;
        }

        const std::string& name() const
        {
          return _name;
        }

        void name(const std::string& name)
        {
          _name = name;
        }

        bool isRootSpace() const
        {
          return true;
        }

        const Basis& basis() const
        {
          return *_df_basis;
        }

      private:

        typename Traits::EntitySet _es;
        std::shared_ptr<DFBasis> _df_basis;
        std::shared_ptr<CE const> _pce;
        Ordering _ordering;
        std::string _name;
      };

    } // namespace Experimental
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_DUNEFUNCTIONSGRIDFUNCTIONSPACE_HH
