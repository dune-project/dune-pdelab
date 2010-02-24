#ifndef DUNE_PDELAB_STATIONARYLINEARPROBLEM_HH
#define DUNE_PDELAB_STATIONARYLINEARPROBLEM_HH

namespace Dune {
  namespace PDELab {

    //===============================================================
    // A class for solving linear stationary problems.
    // It assembles the matrix, computes the right hand side and
    // solves the problem.
    // This is only a first vanilla implementation which has to be improved.
    //===============================================================

    template<class GOS, class LS, class V> 
    class StationaryLinearProblemSolver
    {
      typedef typename V::ElementType Real;
      typedef typename GOS::template MatrixContainer<Real>::Type M;
      typedef typename GOS::Traits::TrialGridFunctionSpace::template VectorContainer<Real>::Type W;

    public:

      StationaryLinearProblemSolver (const GOS& gos_, V& x_, LS& ls_, typename V::ElementType reduction_)
	: gos(gos_), ls(ls_), x(&x_), reduction(reduction_)
      {
      }

      StationaryLinearProblemSolver (const GOS& gos_, LS& ls_, typename V::ElementType reduction_)
	: gos(gos_), ls(ls_), x(0), reduction(reduction_)
      {
      }

      void apply (V& x_) {
        x = &x_;
        apply();
      }
      void apply ()
      {
	// assemble matrix; optional: assemble only on demand!
	M m(gos); 
	m = 0.0;
	gos.jacobian(*x,m);

	// assemble residual
	W r(gos.testGridFunctionSpace(),0.0);
        gos.residual(*x,r);  // residual is additive

	// compute correction
	V z(gos.trialGridFunctionSpace(),0.0);
	ls.apply(m,z,r,reduction); // solver makes right hand side consistent

	// and update
        *x -= z;
      }

    private:
      const GOS& gos;
      LS& ls;
      V* x;
      typename V::ElementType reduction;
    };

  } // namespace PDELab
} // namespace Dune

#endif
