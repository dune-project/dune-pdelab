// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_DUNEFUNCTIONSGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_DUNEFUNCTIONSGRIDFUNCTIONSPACE_HH

#include <cstddef>
#include <map>
#include <bitset>

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>

#include <dune/typetree/leafnode.hh>
#include <dune/typetree/compositenode.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/dunefunctionslocalfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/dunefunctionslfsindexcache.hh>

#include <dune/pdelab/backend/istl/dunefunctions.hh>
#include <dune/pdelab/backend/istl.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

    namespace Experimental {

      // The following code recognizes whether the given VectorBackend (VBE) is an ISTL backend.
      // If this is the case, then we need to replace it by ISTL::SimpleVectorBackend,
      // because we cannot handle anything more complicated at the moment.
      template<typename VBE>
      VBE* registerDuneFunctionsCompatibleVBE(VBE*);

      template<std::size_t block_size>
      ISTL::SimpleVectorBackend<(block_size > 0 ? block_size : 1)>* registerDuneFunctionsCompatibleVBE(ISTL::VectorBackend<ISTL::Blocking::none,block_size>*);

      template<typename VBE>
      using DuneFunctionsCompatibleVBE = std::decay_t<decltype(*registerDuneFunctionsCompatibleVBE(std::declval<VBE*>()))>;

      /** \brief A pdelab grid function space implemented by a dune-functions function space basis
       *
       * \warning This class works only under quite restrictive assumptions:
       *  - The dune-functions basis has to be scalar-valued
       *  - The basis has to be such that the local finite element type for a given element
       *    can be infered from the GeometryType of the element alone.
       *    (Due to a restriction in the current implementation of the FiniteElementMap)
       *
       *  \tparam DFBasis A dune-functions function space basis
       *  \tparam VBE     The type of the underlying vector backend
       *  \tparam CE      Type for constraints assembler
       */
      template<typename DFBasis, typename VBE, typename CE>
      class GridFunctionSpace
        : public TypeTree::LeafNode
        , public GridFunctionOutputParameters
        , public DataHandleProvider<GridFunctionSpace<DFBasis,VBE,CE> >
      {
        using GV = typename DFBasis::GridView;

        template<typename,typename>
        friend class GridFunctionSpaceBase;

      public:
        //! export Traits class

        struct Traits {

          using GridView  = Dune::PDELab::impl::GridView<typename DFBasis::GridView>;
          using GridViewType = GridView;    // DiscreteGridFunction wants this
          using EntitySet = Dune::PDELab::impl::EntitySet<typename DFBasis::GridView>;

          using size_type = std::size_t;
          using SizeType  = size_type;
          using ConstraintsType = CE;

          using FiniteElementType = typename DFBasis::LocalView::Tree::FiniteElement;  // DiscreteGridFunction wants this

          using Basis = DFBasis;

          // The following code recognizes whether the given VectorBackend (VBE) is an ISTL backend.
          // If this is the case, then we replace it by ISTL::SimpleVectorBackend,
          // because we cannot handle anything more complicated at the moment.
          using Backend = DuneFunctionsCompatibleVBE<VBE>;

          /** \brief Rudimentary internal implementation of a FiniteElementMap */
          struct FEM
          {
            struct Traits
            {
              using FiniteElement = typename DFBasis::LocalView::Tree::FiniteElement;
              using FiniteElementType = FiniteElement;
            };

            FEM(const std::shared_ptr<DFBasis>& basis)
            : _basis(basis)
            {}

            /** \brief Get local basis functions for entity
             *
             * This method makes a few short-cuts.  The problem is that dune-functions bases return LocalFiniteElement objects
             * by value, but the method 'find' here hands them out by reference.  Therefore, the FiniteElementMap class
             * has to assume ownership of these objects.  In principle they can be different for each element.  However,
             * storing a LocalFiniteElement for each element can be very expensive.  Therefore, we make the simplifying
             * assumption that the basis is such that the LocalFiniteElement type can be inferred from the GeometryType
             * of the element alone.  This will work for many interesting case, but it will fail, for example, for
             * p-adaptive or XFEM-type bases.
             */
            const typename Traits::FiniteElementType& find (const typename GridView::template Codim<0>::Entity& element) const
            {
              auto type = element.type();
              auto mapEntry = geometryTypeToLocalView_.find(type);
              if (mapEntry == geometryTypeToLocalView_.end())
              {
                auto newLocalView = std::make_shared<typename DFBasis::LocalView>(_basis->localView());
                newLocalView->bind(element);
                auto insertedLocalView = geometryTypeToLocalView_.insert(std::make_pair(type, newLocalView));
                return insertedLocalView.first->second->tree().finiteElement();
              }
              else
              {
                return mapEntry->second->tree().finiteElement();
              }
            }

            void update()
            {
              geometryTypeToLocalView_.clear();
            }

          private:
            const std::shared_ptr<DFBasis> _basis;

            mutable std::map<GeometryType, std::shared_ptr<typename DFBasis::LocalView> > geometryTypeToLocalView_;
          };

          using FiniteElementMap = FEM;
          using FiniteElementMapType = FEM;

        };

        using Basis          = DFBasis;

        /** \brief The actual Ordering object of the grid function space
         *
         * This class is the leaf in the ordering tree of the dune-functions grid function space.
         * It is the class that implements that actual ordering.  PDELab requires orderings to be
         * trees with at least two nodes even if the basis itself is represented by a tree with
         * a single node only.
         */
        struct LeafOrdering
          : public TypeTree::LeafNode
        {

          struct Traits {

            /** \brief A DOF index that is independent of any ordering */
            using DOFIndex       = PDELab::DOFIndex<std::size_t,1,2>;

            /** \brief The index to access containers with
             *
             * The implementation does not support blocking or other forms of multi-indices yet.
             * Therefore the container index type is always a multi-index with one digit,
             * i.e., an integer.
             */
            using ContainerIndex = PDELab::MultiIndex<std::size_t,1>;
            using size_type      = std::size_t;
            using SizeType       = size_type;

            using DOFIndexAccessor = Dune::PDELab::DefaultDOFIndexAccessor;
          };

          using DOFIndex       = typename Traits::DOFIndex;
          using ContainerIndex = typename Traits::ContainerIndex;
          using size_type      = std::size_t;

          LeafOrdering(const GridFunctionSpace& gfs)
            : _gfs(gfs)
          {
            update();
          }

          size_type size() const
          {
            return _gfs.basis().size();
          }

          /** \brief Number of degrees of freedom per entity */
          size_type size(const typename DOFIndex::EntityIndex& entity) const
          {
            return _containerIndices[entity[0]][entity[1]].size();
          }

          size_type maxLocalSize() const
          {
            return _gfs.basis().localView().maxSize();
          }

          /** \brief True if there is at least one entity of the given codim that has a dof
           */
          bool contains(typename Traits::SizeType codim) const
          {
            return _contains[codim];
          }

          /** \brief True if all entities of the given codimension have the same number of dofs
           */
          bool fixedSize(typename Traits::SizeType codim) const
          {
            return _fixedSize[codim];
          }

          // child_index: Steffen sagt: unklar, im Zweifel einfach ignorieren
          template<typename CIOutIterator, typename DIOutIterator = DummyDOFIndexIterator>
          typename Traits::SizeType
          extract_entity_indices(const typename Traits::DOFIndex::EntityIndex& entityIndex,
                                 typename Traits::SizeType child_index,
                                 CIOutIterator ci_out, const CIOutIterator ci_end,
                                 DIOutIterator dummy) const
          {
            for (size_type i=0; i<_containerIndices[entityIndex[0]][entityIndex[1]].size(); i++)
            {
              *ci_out = _containerIndices[entityIndex[0]][entityIndex[1]][i];
              ci_out++;
            }

            return _containerIndices[entityIndex[0]][entityIndex[1]].size();
          }

          ContainerIndex containerIndex(const DOFIndex& i) const
          {
            return _containerIndices[i.entityIndex()[0]][i.entityIndex()[1]][i.treeIndex()[0]];
          }

          void update()
          {
            _containerIndices.clear();

            constexpr auto dim = GV::dimension;
            const auto  gridView = _gfs.gridView();
            const auto& indexSet = gridView.indexSet();

            // Count how many dofs there are for each individual entity
            std::vector<std::vector<size_type> > dofsPerEntity(GlobalGeometryTypeIndex::size(dim));
            for (size_type codim=0; codim<=dim; codim++)
              for (auto&& type : indexSet.types(codim))
              {
                dofsPerEntity[GlobalGeometryTypeIndex::index(type)].resize(gridView.size(type));
                std::fill(dofsPerEntity[GlobalGeometryTypeIndex::index(type)].begin(),
                          dofsPerEntity[GlobalGeometryTypeIndex::index(type)].end(), 0);
              }

            typename DFBasis::LocalView localView = _gfs.basis().localView();
            for (auto&& element : elements(gridView))
            {
              localView.bind(element);
              const auto refElement = ReferenceElements<typename GV::ctype,dim>::general(element.type());

              const auto& localFiniteElement = localView.tree().finiteElement();

              for (size_type i=0; i<localFiniteElement.size(); i++)
              {
                const auto& localKey = localFiniteElement.localCoefficients().localKey(i);

                // Type of the entity that the current local dof is attached to
                auto subentityTypeIndex = GlobalGeometryTypeIndex::index(refElement.type(localKey.subEntity(), localKey.codim()));

                // Global index of the entity that the current local dof is attached to
                auto subentityIndex = indexSet.subIndex(element, localKey.subEntity(), localKey.codim());

                dofsPerEntity[subentityTypeIndex][subentityIndex]
                  = std::max(dofsPerEntity[subentityTypeIndex][subentityIndex], (size_type)localKey.index()+1);
              }
            }

            // Set up the nested container for the container indices
            _containerIndices.resize(GlobalGeometryTypeIndex::size(dim));

            for (size_type codim=0; codim<=dim; codim++)
              for (auto&& type : indexSet.types(codim))
              {
                _containerIndices[GlobalGeometryTypeIndex::index(type)].resize(gridView.size(type));
                for (size_type i=0; i<_containerIndices[GlobalGeometryTypeIndex::index(type)].size(); i++)
                  _containerIndices[GlobalGeometryTypeIndex::index(type)][i].resize(dofsPerEntity[GlobalGeometryTypeIndex::index(type)][i]);
              }

            // Actually set the container indices for all dofs on all entities
            for (auto&& element : elements(gridView))
            {
              localView.bind(element);
              const auto refElement = ReferenceElements<typename GV::ctype,dim>::general(element.type());

              const auto& localFiniteElement = localView.tree().finiteElement();

              for (size_type i=0; i<localFiniteElement.size(); i++)
              {
                const auto& localKey = localFiniteElement.localCoefficients().localKey(i);

                // Type of the entity that the current local dof is attached to
                GeometryType subentityType = refElement.type(localKey.subEntity(), localKey.codim());

                // Global index of the entity that the current local dof is attached to
                auto subentityIndex = indexSet.subIndex(element, localKey.subEntity(), localKey.codim());

                _containerIndices[GlobalGeometryTypeIndex::index(subentityType)][subentityIndex][localKey.index()].set({localView.index(i)});
              }

            }

            // Precompute for each codimension whether there is at least one entity with degrees of freedom,
            // and whether all entities have the same number of dofs
            for (size_type codim=0; codim<=dim; codim++)
            {
              _contains[codim] = false;
              _fixedSize[codim] = true;
              for (auto&& type : indexSet.types(codim))
              {
                const auto& dofs = dofsPerEntity[GlobalGeometryTypeIndex::index(type)];

                if (dofs[0] > 0)
                  _contains[codim] = true;

                for (size_type i=1; i<dofs.size(); i++)
                {
                  if (dofs[i] > 0)
                    _contains[codim] = true;

                  if (dofs[i-1] != dofs[i])
                    _fixedSize[codim] = false;
                }
              }
            }
          }

        private:

          const GridFunctionSpace& _gfs;

          // Container that contains the ContainerIndices for all dofs, accessible by entities
          std::vector<std::vector<std::vector<ContainerIndex> > > _containerIndices;

          /** \brief True if there is at least one entity of the given codim
           *         that has degrees of freedom
           */
          std::bitset<GV::dimension+1> _contains;

          /** \brief True if all entities of the given codim have the same number of dofs
           */
          std::bitset<GV::dimension+1> _fixedSize;
        };

        /** \brief Root of the ordering tree
         *
         * PDELab requires ordering trees to have at least two nodes even if the corresponding
         * basis tree consists of a single node only.  So here is an artificial root node to
         * please PDELab.  All it does is forward all method calls to its single child.
         */
        struct Ordering
          : public TypeTree::CompositeNode<LeafOrdering>
        {
          friend class LocalFunctionSpace<GridFunctionSpace>;

          using Traits = typename LeafOrdering::Traits;

          static const bool consume_tree_index = false;

          using DOFIndex       = typename Traits::DOFIndex;
          using ContainerIndex = typename Traits::ContainerIndex;
          using size_type      = std::size_t;

          using CacheTag       = DuneFunctionsCacheTag;
          using ContainerAllocationTag = FlatContainerAllocationTag;

          Ordering(const GridFunctionSpace& gfs)
            : TypeTree::CompositeNode<LeafOrdering>(LeafOrdering(gfs))
          {}

          size_type size() const
          {
            return this->child(Indices::_0).size();
          }

          /** \brief Same as size(), because block size is always 1
           */
          size_type blockCount() const
          {
            return this->child(Indices::_0).size();
          }

          size_type maxLocalSize() const
          {
            return this->child(Indices::_0).maxLocalSize();
          }

          /** \brief Returns true if there is at least one entity of the given codim
           *         for which data needs to be communicated.
           */
          bool contains(typename Traits::SizeType codim) const
          {
            return this->child(Indices::_0).contains(codim);
          }

          /** \brief True if for all entities of the given codim the same number of data items has to be communicated
           */
          bool fixedSize(typename Traits::SizeType codim) const
          {
            return this->child(Indices::_0).fixedSize(codim);
          }

          template<typename CIOutIterator, typename DIOutIterator = DummyDOFIndexIterator>
          typename Traits::SizeType
          extract_entity_indices(const typename Traits::DOFIndex::EntityIndex& ei,
                                 typename Traits::SizeType child_index,
                                 CIOutIterator ci_out, const CIOutIterator ci_end) const
          {
            return 0;
          }

          void update()
          {
            this->child(Indices::_0).update();
          }

        private:

          ContainerIndex containerIndex(const DOFIndex& i) const
          {
            return this->child(Indices::_0).containerIndex(i);
          }

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
          : _es(df_basis->gridView(), Traits::EntitySet::allCodims())
          , _df_basis(std::move(df_basis))
          , _finiteElementMap(_df_basis)
          , _pce(std::move(ce))
          , _ordering(*this)
        {}

        GridFunctionSpace (std::shared_ptr<DFBasis> df_basis)
          : _es(df_basis->gridView(), Traits::EntitySet::allCodims())
          , _df_basis(std::move(df_basis))
          , _finiteElementMap(_df_basis)
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

        //! get finite element map
        const auto& finiteElementMap () const
        {
          return _finiteElementMap;
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

        typename Traits::SizeType size() const
        {
          return _ordering.size();
        }

        typename Traits::SizeType blockCount() const
        {
          return _ordering.blockCount();
        }

        typename Traits::SizeType globalSize() const
        {
          return _ordering.size();
        }

        typename Traits::SizeType maxLocalSize () const
        {
          return _ordering.maxLocalSize();
        }

        /** \brief Update the indexing information of the GridFunctionSpace.
         *
         * \ param force   Set to true if the underlying grid has changed (e.g. due to adaptivity)
         *                 to force an update of the embedded EntitySet.
         */
        void update(bool force = false)
        {
          _es.update(force);
          _df_basis->update(_es.gridView());
          _finiteElementMap.update();
          // Apparently there is no need to update the constraints assembler '_pce';
          _ordering.update();
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
        typename Traits::FiniteElementMap _finiteElementMap;
        std::shared_ptr<CE const> _pce;
        Ordering _ordering;
        std::string _name;
      };

    } // namespace Experimental

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_DUNEFUNCTIONSGRIDFUNCTIONSPACE_HH
