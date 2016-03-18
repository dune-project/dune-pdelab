#ifndef DUNE_PDELAB_LOCALOPERATOR_PERMEABILITY_ADAPTER_HH
#define DUNE_PDELAB_LOCALOPERATOR_PERMEABILITY_ADAPTER_HH

#include <dune/pdelab/common/function.hh>

/*! Adapter that extracts permeability from parameter class

  \tparam T  model of ConvectionDiffusionParameterInterface
*/
template<typename T>
class PermeabilityAdapter
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                       typename T::Traits::RangeFieldType,
                                       1,Dune::FieldVector<typename T::Traits::RangeFieldType,1> >
                      ,PermeabilityAdapter<T> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                       typename T::Traits::RangeFieldType,
                       1,Dune::FieldVector<typename T::Traits::RangeFieldType,1> > Traits;

  //! constructor
  PermeabilityAdapter (const typename Traits::GridViewType& g_, T& t_)
    : g(g_), t(t_)
  {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
            const typename Traits::DomainType& x,
            typename Traits::RangeType& y) const
  {
    y = log(t.A(e,x)[0][0]);
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return g;
  }

  inline void setTime(double time_)
  {
    t.setTime(time_);
  }

private:
  const typename Traits::GridViewType& g;
  T& t;
};

/*! Adapter that extracts diagonal of permeability tensor from parameter class

  \tparam T  model of ConvectionDiffusionParameterInterface
*/
template<typename T>
class DiagonalPermeabilityAdapter
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                       typename T::Traits::RangeFieldType,
                                       T::Traits::dimDomain,Dune::FieldVector<typename T::Traits::RangeFieldType,T::Traits::dimDomain> >
                      ,DiagonalPermeabilityAdapter<T> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                       typename T::Traits::RangeFieldType,
                       T::Traits::dimDomain,Dune::FieldVector<typename T::Traits::RangeFieldType,T::Traits::dimDomain> > Traits;

  //! constructor
  DiagonalPermeabilityAdapter (const typename Traits::GridViewType& g_, T& t_)
    : g(g_), t(t_)
  {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
            const typename Traits::DomainType& x,
            typename Traits::RangeType& y) const
  {
    for (int i=0; i<T::Traits::dimDomain; i++)
      y[i] = log10(t.A(e,x)[i][i]);
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return g;
  }

  inline void setTime(double time_)
  {
    t.setTime(time_);
  }

private:
  const typename Traits::GridViewType& g;
  T& t;
};

#endif // DUNE_PDELAB_LOCALOPERATOR_PERMEABILITY_ADAPTER_HH
