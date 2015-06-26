#include<dune/common/fvector.hh>

#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>

//! \brief Parameter class selecting boundary conditions
class BCTypeParam
  : public Dune::PDELab::DirichletConstraintsParameters
{
public:
  //! Test whether boundary is Dirichlet-constrained
  template<typename I>
  bool isDirichlet(const I & intersection,
                   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                   ) const
  {
    Dune::FieldVector<typename I::ctype, I::dimension>
      xg = intersection.geometry().global( coord );

    if( xg[0]>1.0-1E-6 )
      return false; // no Dirichlet b.c. on the eastern boundary

    return true;  // Dirichlet b.c. on all other boundaries
  }

};
