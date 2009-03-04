// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_GRIDOPERATORSPACEUTILITIES_HH
#define DUNE_PDELAB_GRIDOPERATORSPACEUTILITIES_HH

#include<dune/common/exceptions.hh>

#include"../common/geometrywrapper.hh"

namespace Dune {
  namespace PDELab {

	// compile time switching of function call
    template<typename LA, bool doIt>
    struct LocalAssemblerCallSwitch
    {
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
	  {
	  }
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_volume_apply (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y)
	  {
	  }
    };
    template<typename LA>
    struct LocalAssemblerCallSwitch<LA,true>
    {
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
	  {
		la.alpha_volume(eg,lfsu,x,lfsv,r);
	  }
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_volume_apply (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y)
	  {
		la.jacobian_volume_apply(eg,lfsu,x,lfsv,y);
	  }
    };


	// derive from this class to add numerical evaluation of jacobian apply
	template<typename Imp>
	class NumericalJacobianVolumeApply
	{
	public:

	  // provide numerical evaluation of jacobian
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
	  void jacobian_volume_apply (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y) const
	  {
		typedef typename X::value_type R;
		const R epsilon=1E-8; // problem: this depends on data type R!
		const int m=lfsv.size();
		const int n=lfsu.size();

		X u(x);
		std::vector<R> down(m,0.0),up(m);

		y.resize(m);
		asImp().alpha_volume(eg,lfsu,u,lfsv,down);
		for (int j=0; j<n; j++) // loop over columns
		  {
			for (int k=0; k<m; k++) up[k]=0.0;
			R delta = epsilon*(1.0+std::abs(u[j]));
			u[j] += delta;
			asImp().alpha_volume(eg,lfsu,u,lfsv,up);
			for (int i=0; i<m; i++)
			  y[i] += ((up[i]-down[i])/delta)*x[j];
			u[j] = x[j];
		  }
	  }

	private:
	  Imp& asImp () {return static_cast<Imp &> (*this);}
	  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
	};

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
