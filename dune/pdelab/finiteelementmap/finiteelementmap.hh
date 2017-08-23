// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENTMAP_FINITEELEMENTMAP_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_FINITEELEMENTMAP_HH

#include <dune/common/deprecated.hh>
#include <dune/pdelab/common/exceptions.hh>

#include <dune/geometry/referenceelements.hh>

namespace Dune {
  namespace PDELab {

    //! \ingroup FiniteElementMap
    //! \{

    //! \brief general FiniteElementMap exception
    class FiniteElementMapError : public Exception {};
    //! \brief FiniteElementMap exception concerning the computation of the FiniteElementMap size
    class VariableElementSize : public FiniteElementMapError {};
    //! \brief FiniteElementMap exception raised when trying to obtain a finite element for an unsupported GeometryType.
    class InvalidGeometryType : public FiniteElementMapError {};

    //! \brief collect types exported by a finite element map
    template<class T>
    struct FiniteElementMapTraits
    {
      //! Type of finite element from local functions
      typedef T FiniteElementType;

      //! Type of finite element from local functions
      typedef T FiniteElement;
    };

    //! \brief collect types exported by a finite element map
    template<class T>
    struct LocalFiniteElementMapTraits : FiniteElementMapTraits<T> {};

    //! \brief interface for a finite element map
    template<class T, class Imp>
    class LocalFiniteElementMapInterface
    {
    public:
      //! \brief Export traits
      typedef T Traits;

      /** \brief Return local basis for the given entity.

          The return value is a reference to Traits::LocalBasisType. If
          there is a different local basis for two elements then this
          type must be polymorphic.
      */
      template<class EntityType>
      const typename Traits::FiniteElementType&
      find (const EntityType& e) const = delete;

      /** @name Size calculation
       *  The FiniteElementMap provides different methods to compute
       *  the size of the GridFunctionSpace (if possible) without
       *  iterating the grid. The approach is as follows (pseudo code):
       *
       *  \code
       *  computeNumberOfDofs(GridView, FEM):
       *    if(FEM.fixedSize()):
       *      sum(FEM.size(gt)*GridView.size(gt) for gt in GeometryTypes)
       *    else
       *      sum(FEM.find(E).basis().size() for E in GridView.entities<0>())
       *  \endcode
       */
      /** @{ */
      /** \brief a FiniteElementMap is fixedSize iif the size of the local
       * functions space for each GeometryType is fixed.
       */
      bool fixedSize() const = delete;
      /** \brief if the FiniteElementMap is fixedSize, the size
       * methods computes the number of DOFs for given GeometryType.
       */
      std::size_t size(GeometryType gt) const = delete;
      /** @} */
      /** \brief return if FiniteElementMap has degrees of freedom for
       * given codimension
       */
      bool hasDOFs(int codim) const = delete;
      /** \brief compute an upper bound for the local number of DOFs.
       *
       * this upper bound is used to avoid reallocations in
       * std::vectors used during the assembly.
       */
      std::size_t maxLocalSize() const = delete;
    };

    namespace Impl {

      template<int dim>
      struct inject_dimension
      {
        static constexpr int dimension = dim;
      };

      template<>
      struct inject_dimension<-1>
      {
        DUNE_DEPRECATED_MSG("Your FEM implementation does not export the dimension. This will cause compilation errors after PDELab 2.6.") inject_dimension() {};
      };

    }

    //! simple implementation where all entities have the same finite element
    template<class Imp, int dim_ = -1>
    class SimpleLocalFiniteElementMap :
      public LocalFiniteElementMapInterface<LocalFiniteElementMapTraits<Imp>,
                                            SimpleLocalFiniteElementMap<Imp> >,
      public Impl::inject_dimension<dim_>
    {
    public:
      //! \brief export type of the signature
      typedef LocalFiniteElementMapTraits<Imp> Traits;

      //! \brief Use when Imp has a standard constructor
      SimpleLocalFiniteElementMap ()
      {}

      //! \brief Constructor where an instance of Imp can be provided
      SimpleLocalFiniteElementMap (const Imp& imp_) : imp(imp_)
      {}

      //! \brief get local basis functions for entity
      template<class EntityType>
      const typename Traits::FiniteElementType& find (const EntityType& e) const
      {
        return imp;
      }

    private:
      Imp imp; // create once
    };

    /** \brief implementation for finite elements requiring oriented edges
     *
     * This is for edge elements. It works for one type of Geometry only, and
     * the requirements for the local finite element are:
     *
     *  - The default orientation of the shape functions on the edges must be
     *    from the vertex with the lower index within the reference element to
     *    the vertex with the higher index.
     *  - The local finite element must allow assignment.
     *  - The local finite element must have a default constructor.
     *  - the local finite element hust have a constructor of the form
     *    FE(unsigned int orientation), where orientation is a bitfield.  If
     *    bit i in orientation is set it means that the orientation of the
     *    shape function corresponding to the edge with id i in the reference
     *    element is inverted.
     *
     * \tparam GV  Type of gridview to work with
     * \tparam FE  Type of local finite element
     * \tparam Imp Type of the final LocalFiniteElementMap implementation
     */
    template<typename GV, typename FE, typename Imp>
    class EdgeS0LocalFiniteElementMap
      : public LocalFiniteElementMapInterface<
          LocalFiniteElementMapTraits<FE>,
          Imp
        >
    {
      typedef typename GV::IndexSet IndexSet;
      static const int dim = GV::dimension;

    public:
      //! \brief export type of the signature
      typedef LocalFiniteElementMapTraits<FE> Traits;

      //! The dimension of the finite elements returned by this map.
      static constexpr int dimension = GV::dimension;

      //! \brief construct EdgeSLocalFiniteElementMap
      EdgeS0LocalFiniteElementMap (const GV& gv_)
        : gv(gv_), orient(gv_.size(0))
      {
        using ct = typename GV::Grid::ctype;
        auto& refElem =
          ReferenceElements<ct, dim>::general(FE().type());

        auto &idSet = gv.grid().globalIdSet();

        // create all variants
        variant.resize(1 << refElem.size(dim-1));
        for (std::size_t i=0; i<variant.size(); i++)
          variant[i] = FE(i);

        // compute orientation for all elements
        auto& indexSet = gv.indexSet();

        // loop once over the grid
        for (const auto& element : elements(gv))
          {
            auto elemid = indexSet.index(element);
            orient[elemid] = 0;

            std::vector<typename GV::Grid::GlobalIdSet::IdType> vid(refElem.size(dim));
            for(std::size_t i = 0; i < vid.size(); ++i)
              vid[i] = idSet.subId(element, i, dim);

            // loop over all edges of the element
            for(int i = 0; i < refElem.size(dim-1); ++i) {
              auto v0 = refElem.subEntity(i, dim-1, 0, dim);
              auto v1 = refElem.subEntity(i, dim-1, 1, dim);
              // if (edge orientation in refelement) != (edge orientation in indexset)
              if((v0 > v1) != (vid[v0] > vid[v1]))
                orient[elemid] |= 1 << i;
            }
          }
      }

      //! \brief get local basis functions for entity
      template<class EntityType>
      const typename Traits::FiniteElementType& find (const EntityType& e) const
      {
        return variant[orient[gv.indexSet().index(e)]];
      }

    private:
      GV gv;
      std::vector<FE> variant;
      std::vector<unsigned char> orient;
    };

    //! wrap up element from local functions
    //! \ingroup FiniteElementMap
    template<typename GV, typename FE, typename Imp, std::size_t Variants>
    class RTLocalFiniteElementMap :
      public LocalFiniteElementMapInterface<
        LocalFiniteElementMapTraits<FE>,
        Imp >
    {
      typedef FE FiniteElement;
      typedef typename GV::IndexSet IndexSet;

    public:
      //! \brief export type of the signature
      typedef LocalFiniteElementMapTraits<FE> Traits;

      //! The dimension of the finite elements returned by this map.
      static constexpr int dimension = GV::dimension;

      //! \brief Use when Imp has a standard constructor
      RTLocalFiniteElementMap(const GV& gv_)
        : gv(gv_), is(gv_.indexSet()), orient(gv_.size(0))
      {
        // create all variants
        for (std::size_t i = 0; i < Variants; i++)
        {
          variant[i] = FiniteElement(i);
        }

        // compute orientation for all elements
        // loop once over the grid
        for(const auto& cell : elements(gv))
        {
          auto myId = is.index(cell);
          orient[myId] = 0;

          for (const auto& intersection : intersections(gv,cell))
          {
            if (intersection.neighbor()
                && is.index(intersection.outside()) > myId)
            {
              orient[myId] |= 1 << intersection.indexInInside();
            }
          }
        }
      }

      //! \brief get local basis functions for entity
      template<class EntityType>
      const typename Traits::FiniteElementType& find(const EntityType& e) const
      {
        return variant[orient[is.index(e)]];
      }

    private:
      GV gv;
      std::array<FiniteElement,Variants> variant;
      const IndexSet& is;
      std::vector<unsigned char> orient;
    };

    //! \} group FiniteElementMap

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_FINITEELEMENTMAP_HH
