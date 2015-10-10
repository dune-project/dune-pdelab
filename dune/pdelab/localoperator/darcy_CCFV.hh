// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DarcyVelocityFromHeadCCFV_HH
#define DarcyVelocityFromHeadCCFV_HH

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
  typedef typename GV::IndexSet IndexSet;
  typedef typename IndexSet::IndexType IndexType;

  GV gv;
  V& c;
  const IndexSet& indexSet;

public:
  //! export type of data for message buffer
  typedef typename V::value_type DataType;

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
  typedef typename PL::Traits::GridViewType GV;
  typedef typename GV::IndexSet IndexSet;
  typedef typename GV::Grid::ctype DF;
  typedef typename PL::Traits::RangeFieldType RF;
  typedef typename PL::Traits::RangeType RangeType;
  enum { dim = PL::Traits::GridViewType::dimension };
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;
  typedef typename Dune::RaviartThomasCubeLocalFiniteElement<DF,RF,dim,0>::Traits::LocalBasisType::Traits::RangeType RT0RangeType;
  typedef typename Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

  const T& t;
  const PL& pl;
  Dune::RaviartThomasCubeLocalFiniteElement<DF,RF,dim,0> rt0fe;

  mutable Dune::FieldMatrix<DF,dim,dim> B;
  mutable RF determinant;
  mutable int cachedindex;
  typename T::Traits::RangeFieldType time;

  typedef Dune::FieldVector<RF,2*dim> RT0Coeffs;
  GV gv;
  const IndexSet& is;
  std::vector<RT0Coeffs> storedcoeffs;
  mutable std::vector<RT0RangeType> rt0vectors;

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,dim,Dune::FieldVector<RF,dim> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Traits,DarcyVelocityFromHeadCCFV<T,PL> > BaseT;

  DarcyVelocityFromHeadCCFV (const T& t_, const PL& pl_)
    : t(t_), pl(pl_), cachedindex(-1), time(0), gv(pl_.getGridView()), is(gv.indexSet()), storedcoeffs(is.size(0)),
      rt0vectors(rt0fe.localBasis().size())
  {
    // iterate over grid and store values
    typedef typename GV::Traits::template Codim<0>::template Partition<Dune::Interior_Partition>::Iterator ElementIterator;
    //std::cout << "Allocated std::vector with size " << storedcoeffs.size() << std::endl;

    // compute RT0 coefficients for all interior cells
    for (ElementIterator it = gv.template begin<0,Dune::Interior_Partition>();
         it!=gv.template end<0,Dune::Interior_Partition>(); ++it)
      {
        // get local cell number
        int index = is.index(*it);

        // cell geometry
        const Dune::FieldVector<DF,dim>
          inside_cell_center_local = Dune::ReferenceElements<DF,dim>::
          general(it->type()).position(0,0);
        Dune::FieldVector<DF,dim>
          inside_cell_center_global = it->geometry().global(inside_cell_center_local);

        // absolute permeability in primary cell
        typename T::Traits::PermTensorType tensor_inside;
        tensor_inside = t.A(*it,inside_cell_center_local);

        // pressure evaluation
        typename PL::Traits::RangeType pl_inside;
        pl.evaluate(*it,inside_cell_center_local,pl_inside);

        // for coefficient computation
        RF vn[2*dim];    // normal velocities
        Dune::FieldMatrix<typename Traits::DomainFieldType,dim,dim>
          B = it->geometry().jacobianInverseTransposed(inside_cell_center_local); // the transformation. Assume it is linear
        RF determinant = B.determinant();

        // loop over cell neighbors
        IntersectionIterator endit = pl.getGridView().iend(*it);
        for (IntersectionIterator iit = pl.getGridView().ibegin(*it); iit!=endit; ++iit)
          {
            // set to zero for processor boundary
            vn[iit->indexInInside()] = 0.0;

            // face geometry
            const Dune::FieldVector<DF,dim-1>&
              face_local = Dune::ReferenceElements<DF,dim-1>::general(iit->geometry().type()).position(0,0);

            // interior face
            if (iit->neighbor())
              {
                auto outside_cell = iit->outside();
                const Dune::FieldVector<DF,dim>
                  outside_cell_center_local = Dune::ReferenceElements<DF,dim>::
                  general(outside_cell.type()).position(0,0);
                Dune::FieldVector<DF,dim>
                  outside_cell_center_global = outside_cell.geometry().global(outside_cell_center_local);

                // distance of cell centers
                Dune::FieldVector<DF,dim> d(outside_cell_center_global);
                d -= inside_cell_center_global;
                RF distance = d.two_norm();

                // absolute permeability
                typename T::Traits::PermTensorType tensor_outside;
                tensor_outside = t.A(outside_cell,outside_cell_center_local);
                const Dune::FieldVector<DF,dim> n_F = iit->centerUnitOuterNormal();
                Dune::FieldVector<RF,dim> An_F;
                tensor_inside.mv(n_F,An_F);
                RF k_inside = n_F*An_F;
                tensor_outside.mv(n_F,An_F);
                RF k_outside = n_F*An_F;
                RF k_avg = 2.0/(1.0/(k_inside+1E-30) + 1.0/(k_outside+1E-30));

                // pressure evaluation
                typename PL::Traits::RangeType pl_outside;
                pl.evaluate(outside_cell,outside_cell_center_local,pl_outside);

                // set coefficient
                vn[iit->indexInInside()] = k_avg*(pl_inside-pl_outside)/distance;
              }

            // boundary face
            if (iit->boundary())
              {
                // distance of cell center to boundary
                Dune::FieldVector<DF,dim> d = iit->geometry().global(face_local);
                d -= inside_cell_center_global;
                RF distance = d.two_norm();

                // evaluate boundary condition type
                BCType bctype;
                bctype = t.bctype(*iit,face_local);

                // liquid phase Dirichlet boundary
                if (bctype==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet)
                  {
                    Dune::FieldVector<DF,dim> iplocal_s = iit->geometryInInside().global(face_local);
                    RF g_l = t.g(iit->inside(),iplocal_s);
                    const Dune::FieldVector<DF,dim> n_F = iit->centerUnitOuterNormal();
                    Dune::FieldVector<RF,dim> An_F;
                    tensor_inside.mv(n_F,An_F);
                    RF k_inside = n_F*An_F;
                    vn[iit->indexInInside()] = k_inside * (pl_inside-g_l)/distance;
                  }

                // liquid phase Neumann boundary
                if (bctype==Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann)
                  {
                    typename T::Traits::RangeFieldType j = t.j(*iit,face_local);
                    vn[iit->indexInInside()] = j;
                  }
              }

            // compute coefficient
            Dune::FieldVector<DF,dim> vstar=iit->centerUnitOuterNormal(); // normal on tranformef element
            vstar *= vn[iit->indexInInside()];
            Dune::FieldVector<RF,dim> normalhat(0); // normal on reference element
            if (iit->indexInInside()%2==0)
              normalhat[iit->indexInInside()/2] = -1.0;
            else
              normalhat[iit->indexInInside()/2] =  1.0;
            Dune::FieldVector<DF,dim> vstarhat(0);
            B.umtv(vstar,vstarhat); // Piola backward transformation
            vstarhat *= determinant;
            storedcoeffs[index][iit->indexInInside()] = vstarhat*normalhat;
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

#endif
