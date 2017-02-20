// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_TEST_TESTCOMPLEXNUMBERS_PROBLEM_HH
#define DUNE_PDELAB_TEST_TESTCOMPLEXNUMBERS_PROBLEM_HH

#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

/** \file

    \brief Parameters for spherical wave problem.
    \author Philipp Stekl, Marian Piatkowski
*/
template<typename GV, typename RF, typename DF>
class ParametersSphericalWave  : public Dune::PDELab::DirichletConstraintsParameters
{
public:

  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  ParametersSphericalWave (double omega_)
    : omega(omega_)
  {
  }


  //! Dirichlet boundary condition type function
  template<typename IG>
  inline bool isDirichlet(const IG& intersection, const Dune::FieldVector<typename IG::ctype, GV::dimension-1>& xlocal) const
  {
    return false;
  }

  //! Neumann boundary condition type function
  template<typename IG>
  inline bool isNeumann(const IG& intersection, const Dune::FieldVector<typename IG::ctype, GV::dimension-1>& xlocal) const
  {
    return true;
  }

  //! Dirichlet boundary condition value
  RF g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal ) const
  {
    return RF(0.,0.);
  }

  //! Neumann boundary condition value
  RF j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& xlocal, RF u) const
  {
    RF i(0., 1.);
    return -i*omega*u;
  }



  //! source term
  RF f (const typename Traits::ElementType& e, const typename Traits::DomainType& position) const
  {
    auto xg = e.geometry().global(position);
    Dune::FieldVector<DF,Traits::dimDomain> midpoint(0.5);

    xg -= midpoint;

    RF res(0.,0.);

    auto r = 0.2;
    auto tmp = xg.two_norm();
    if(tmp < r)
      res = RF(25,0.0);

    return res;
  }


#if 0
  //! exact solution / analytic solution
  RF ua (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal ) const
  {

    // const int dim = Traits::GridViewType::Grid::dimension;
    // typedef typename Traits::GridViewType::Grid::ctype ctype;
    // Dune::FieldVector<ctype,dim> xglobal = e.geometry().global(xlocal);

    std::cout<<"Error: Analytical Solution has not been implemented"<<std::endl;

    return RF(0.,0.);

  }

  const static int dim = Traits::GridViewType::Grid::dimension;

  //! exact grad of solution / analytic grad of solution
  Dune::FieldVector<RF,dim> ugrada (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal ) const
  {


    // typedef typename Traits::GridViewType::Grid::ctype ctype;
    // Dune::FieldVector<ctype,dim> xglobal = e.geometry().global(xlocal);

    std::cout<<"Error: Analytical Solution has not been implemented"<<std::endl;
    Dune::FieldVector<RF,dim> res;
    res[0] = 0.;
    res[1] = 0.;
    return res;

  }
#endif

  const double omega;

};

#endif // DUNE_PDELAB_TEST_TESTCOMPLEXNUMBERS_PROBLEM_HH
