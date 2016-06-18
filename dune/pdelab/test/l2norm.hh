#ifndef DUNE_PDELAB_TEST_L2NORM_HH
#define DUNE_PDELAB_TEST_L2NORM_HH

#include <cmath>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>


// Calculate the squared L2 norm of a function
template<typename U>
typename U::Traits::RangeFieldType
l2norm2 (const U& u, int qorder=1)
{
  // constants and types
  typedef typename U::Traits::GridViewType GV;
  static const int dimDomain = U::Traits::dimDomain;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename U::Traits::DomainFieldType ct;
  typedef Dune::QuadratureRule<ct,dimDomain> QR;
  
  // loop over grid view
  typename U::Traits::RangeFieldType sum = 0.0;
  for (ElementIterator eit = u.getGridView().template begin<0>();
       eit!=u.getGridView().template end<0>(); ++eit)
  {

    Dune::GeometryType gt = eit->geometry().type();
    const QR& rule = Dune::QuadratureRules<ct,dimDomain>::rule(gt,qorder);

    for (typename QR::const_iterator qit=rule.begin(); qit!=rule.end(); ++qit)
    {
      // evaluate the given grid functions at integration point
      typename U::Traits::RangeType u_val;
      u.evaluate(*eit,qit->position(),u_val);

      // accumulate error
      sum += u_val.two_norm2()*qit->weight()*
        eit->geometry().integrationElement(qit->position());
    }
  }
  return sum;
}

// Calculate the L2 norm of a function
template<typename U> 
typename U::Traits::RangeFieldType
l2norm (const U& u, int qorder=1)
{
  return std::sqrt(l2norm2(u, qorder));
}

#endif // DUNE_PDELAB_TEST_L2NORM_HH
