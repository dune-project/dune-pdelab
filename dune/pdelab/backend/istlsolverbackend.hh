#ifndef DUNE_ISTLSOLVERBACKEND_HH
#define DUNE_ISTLSOLVERBACKEND_HH

#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>

namespace Dune {
  namespace PDELab {

	template<typename X, typename Y, typename GOS>
	class OnTheFlyOperator : public Dune::LinearOperator<X,Y>
	{
	public:
	  typedef X domain_type;
	  typedef Y range_type;
	  typedef typename X::field_type field_type;

	  enum {category=Dune::SolverCategory::sequential};

	  OnTheFlyOperator (GOS& gos_)
		: gos(gos_)
	  {}

	  virtual void apply (const X& x, Y& y) const
	  {
		y = 0.0;
		gos.jacobian_apply(x,y);
	  }

	  virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
	  {
		Y temp(y);
		temp = 0.0;
		gos.jacobian_apply(x,temp);
		y.axpy(alpha,temp);
	  }

	private:
	  GOS& gos;
	};

  } // namespace PDELab
} // namespace Dune

#endif
