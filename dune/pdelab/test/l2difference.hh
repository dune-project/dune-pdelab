#ifndef DUNE_PDELAB_TEST_L2DIFFERENCE_HH
#define DUNE_PDELAB_TEST_L2DIFFERENCE_HH

#include <cmath>

#include <dune/common/geometrytype.hh>

#include <dune/grid/common/quadraturerules.hh>


// Calculate the L2 differerence of two functions
template<typename GV, typename U, typename V> 
double l2difference (const GV & gv, const U& u, const V &v, int qorder=1)
{
  // constants and types
  const int dim = GV::Grid::dimension;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Grid::ctype ct;
  
  // loop over grid view
  double sum = 0.0;
  for (ElementIterator eit = gv.template begin<0>();
       eit!=gv.template end<0>(); ++eit)
  {

    Dune::GeometryType gt = eit->geometry().type();
    const Dune::QuadratureRule<ct,dim>& 
      rule = Dune::QuadratureRules<ct,dim>::rule(gt,qorder);

    for (typename Dune::QuadratureRule<ct,dim>::const_iterator qit=rule.begin();
         qit!=rule.end(); ++qit)
    {
      // evaluate the given grid functions at integration point
      typename U::Traits::RangeType u_val;
      u.evaluate(*eit,qit->position(),u_val);

      typename V::Traits::RangeType v_val;
      v.evaluate(*eit,qit->position(),v_val);

      // accumulate error
      v_val -= u_val;
      sum += v_val.two_norm2()*qit->weight()*
        eit->geometry().integrationElement(qit->position());
    }
  }
  return std::sqrt(sum);
}

#endif // DUNE_PDELAB_TEST_L2DIFFERENCE_HH
