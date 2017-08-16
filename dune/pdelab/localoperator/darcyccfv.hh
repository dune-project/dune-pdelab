// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_DARCYCCFV_HH
#define DUNE_PDELAB_LOCALOPERATOR_DARCYCCFV_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/typetraits.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/localfunctions/raviartthomas/raviartthomascube.hh>

#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

// A DataHandle class to exchange RT0 coefficients in the overlap
template<class GV, class V> // mapper type and vector type
class VectorExchange
  : public Dune::CommDataHandleIF<VectorExchange<GV,V>,
                                  typename V::value_type>
{
  using IndexSet = typename GV::IndexSet;
  using IndexType = typename IndexSet::IndexType;

  GV gv;
  V& c;
  const IndexSet& indexSet;

public:
  //! export type of data for message buffer
  using DataType = typename V::value_type;

  //! constructor
  VectorExchange (const GV& gv_, V& c_)
    : gv(gv_), c(c_), indexSet(gv.indexSet())
  {}

  //! returns true if data for this codim should be communicated
  bool contains (int dim, int codim) const
  {
    return (codim==0);
  }

  //! returns true if size per entity of given dim and codim is a constant
  bool fixedsize (int dim, int codim) const
  {
    return true;
  }

  /*! how many objects of type DataType have to be sent for a given entity

    Note: Only the sender side needs to know this size.
  */
  template<class EntityType>
  size_t size (EntityType& e) const
  {
    return 1;
  }

  //! pack data from user to message buffer
  template<class MessageBuffer, class EntityType>
  void gather (MessageBuffer& buff, const EntityType& e) const
  {
    buff.write(c[indexSet.index(e)]);
  }

  /*! unpack data from message buffer to user

    n is the number of objects sent by the sender
  */
  template<class MessageBuffer, class EntityType>
  void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
  {
    DataType x;
    buff.read(x);
    c[indexSet.index(e)]=x;
  }
};

/** \brief Provide velocity field for liquid phase

    Uses RT0 interpolation on a cell.

    - T  : provides TwoPhaseParameterInterface
    - PL : P0 function for liquid phase pressure
*/
template<typename  T, typename PL>
class DarcyVelocityFromHeadCCFV
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                                                           typename PL::Traits::RangeFieldType,
                                                                           PL::Traits::GridViewType::dimension,
                                                                           Dune::FieldVector<typename PL::Traits::RangeFieldType,PL::Traits::GridViewType::dimension> >,
                                          DarcyVelocityFromHeadCCFV<T,PL> >
{
  // extract useful types
  using GV = typename PL::Traits::GridViewType;
  using IndexSet = typename GV::IndexSet;
  using DF = typename GV::Grid::ctype;
  using RF = typename PL::Traits::RangeFieldType;
  using RangeType = typename PL::Traits::RangeType;
  enum { dim = PL::Traits::GridViewType::dimension };
  using Element = typename GV::Traits::template Codim<0>::Entity;
  using IntersectionIterator = typename GV::IntersectionIterator;
  using Intersection = typename IntersectionIterator::Intersection;
  using RT0RangeType = typename Dune::RaviartThomasCubeLocalFiniteElement<DF,RF,dim,0>::Traits::LocalBasisType::Traits::RangeType;
  using BCType = typename Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type;

  const T& t;
  const PL& pl;
  Dune::RaviartThomasCubeLocalFiniteElement<DF,RF,dim,0> rt0fe;

  mutable Dune::FieldMatrix<DF,dim,dim> B;
  mutable RF determinant;
  mutable int cachedindex;
  typename T::Traits::RangeFieldType time;

  using RT0Coeffs = Dune::FieldVector<RF,2*dim>;
  GV gv;
  const IndexSet& is;
  std::vector<RT0Coeffs> storedcoeffs;
  mutable std::vector<RT0RangeType> rt0vectors;

public:
  using Traits = Dune::PDELab::GridFunctionTraits<GV,RF,dim,Dune::FieldVector<RF,dim> >;
  using BaseT = Dune::PDELab::GridFunctionBase<Traits,DarcyVelocityFromHeadCCFV<T,PL> >;

  DarcyVelocityFromHeadCCFV (const T& t_, const PL& pl_)
    : t(t_), pl(pl_), cachedindex(-1), time(0), gv(pl_.getGridView()), is(gv.indexSet()), storedcoeffs(is.size(0)),
      rt0vectors(rt0fe.localBasis().size())
  {
    // compute RT0 coefficients for all interior cells
    for (const auto& element : elements(gv,Dune::Partitions::interior))
      {
        // get local cell number
        int index = is.index(element);

        // get geometry
        auto geo = element.geometry();

        // cell geometry
        auto ref_el = referenceElement(geo);
        auto inside_cell_center_local = ref_el.position(0,0);
        auto inside_cell_center_global = geo.global(inside_cell_center_local);

        // absolute permeability in primary cell
        auto tensor_inside = t.A(element,inside_cell_center_local);

        // pressure evaluation
        typename PL::Traits::RangeType pl_inside;
        pl.evaluate(element,inside_cell_center_local,pl_inside);

        // for coefficient computation
        RF vn[2*dim];    // normal velocities
        auto B = geo.jacobianInverseTransposed(inside_cell_center_local); // the transformation. Assume it is linear
        auto determinant = B.determinant();

        // loop over cell neighbors
        for (const auto& intersection : intersections(pl.getGridView(),element))
          {
            // set to zero for processor boundary
            vn[intersection.indexInInside()] = 0.0;

            // face geometry
            auto face_local = referenceElement(intersection.geometry()).position(0,0);

            // interior face
            if (intersection.neighbor())
              {
                auto outside_cell = intersection.outside();
                auto outside_cell_center_local = referenceElement(outside_cell.geometry()).position(0,0);
                auto outside_cell_center_global = outside_cell.geometry().global(outside_cell_center_local);

                // distance of cell centers
                auto d(outside_cell_center_global);
                d -= inside_cell_center_global;
                auto distance = d.two_norm();

                // absolute permeability
                auto tensor_outside = t.A(outside_cell,outside_cell_center_local);
                auto n_F = intersection.centerUnitOuterNormal();
                Dune::FieldVector<RF,dim> An_F;
                tensor_inside.mv(n_F,An_F);
                auto k_inside = n_F*An_F;
                tensor_outside.mv(n_F,An_F);
                auto k_outside = n_F*An_F;
                auto k_avg = 2.0/(1.0/(k_inside+1E-30) + 1.0/(k_outside+1E-30));

                // pressure evaluation
                typename PL::Traits::RangeType pl_outside;
                pl.evaluate(outside_cell,outside_cell_center_local,pl_outside);

                // set coefficient
                vn[intersection.indexInInside()] = k_avg*(pl_inside-pl_outside)/distance;
              }

            // boundary face
            if (intersection.boundary())
              {
                // distance of cell center to boundary
                auto d = intersection.geometry().global(face_local);
                d -= inside_cell_center_global;
                auto distance = d.two_norm();

                // evaluate boundary condition type
                auto bctype = t.bctype(intersection,face_local);

                // liquid phase Dirichlet boundary
                if (bctype==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet)
                  {
                    auto iplocal_s = intersection.geometryInInside().global(face_local);
                    auto g_l = t.g(intersection.inside(),iplocal_s);
                    auto n_F = intersection.centerUnitOuterNormal();
                    Dune::FieldVector<RF,dim> An_F;
                    tensor_inside.mv(n_F,An_F);
                    auto k_inside = n_F*An_F;
                    vn[intersection.indexInInside()] = k_inside * (pl_inside-g_l)/distance;
                  }

                // liquid phase Neumann boundary
                if (bctype==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann)
                  {
                    auto j = t.j(intersection,face_local);
                    vn[intersection.indexInInside()] = j;
                  }
              }

            // compute coefficient
            auto vstar=intersection.centerUnitOuterNormal(); // normal on tranformef element
            vstar *= vn[intersection.indexInInside()];
            Dune::FieldVector<RF,dim> normalhat(0); // normal on reference element
            if (intersection.indexInInside()%2==0)
              normalhat[intersection.indexInInside()/2] = -1.0;
            else
              normalhat[intersection.indexInInside()/2] =  1.0;
            Dune::FieldVector<DF,dim> vstarhat(0);
            B.umtv(vstar,vstarhat); // Piola backward transformation
            vstarhat *= determinant;
            storedcoeffs[index][intersection.indexInInside()] = vstarhat*normalhat;
          }
      }

    // communicate coefficients in overlap
    VectorExchange<GV,std::vector<RT0Coeffs> > dh(gv,storedcoeffs);
    if (gv.grid().comm().size()>1)
      gv.grid().communicate(dh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
  }

  // set time where operator is to be evaluated (i.e. end of the time intervall)
  void set_time (typename T::Traits::RangeFieldType time_)
  {
    time = time_;
  }

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    // local cell number
    int index = is.index(e);

    // compute velocity on reference element
    rt0fe.localBasis().evaluateFunction(x,rt0vectors);
    typename Traits::RangeType yhat(0);
    for (unsigned int i=0; i<rt0fe.localBasis().size(); i++)
      yhat.axpy(storedcoeffs[index][i],rt0vectors[i]);

    // apply Piola transformation
    if (index != cachedindex)
      {
        B = e.geometry().jacobianTransposed(x); // the transformation. Assume it is linear
        determinant = B.determinant();
        cachedindex = index;
      }
    y = 0;
    B.umtv(yhat,y);
    y *= determinant;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return pl.getGridView();
  }
};

#endif // DUNE_PDELAB_LOCALOPERATOR_DARCYCCFV_HH
