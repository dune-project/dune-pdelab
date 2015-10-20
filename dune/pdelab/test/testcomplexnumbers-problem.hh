// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_TEST_TESTCOMPLEXNUMBERS_PROBLEM_HH
#define DUNE_PDELAB_TEST_TESTCOMPLEXNUMBERS_PROBLEM_HH


/** \file

    \brief Parameters for plane wave problem.
    \author Philipp Stekl, Marian Piatkowski.
*/

template<typename GV, typename RF, typename DF>
class ParametersPlaneWave  : public Dune::PDELab::DirichletConstraintsParameters
{
public:

  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;
  // typedef DiffusionParameterTraits<GV,RF> Traits;

  ParametersPlaneWave (double omega_, double theta_)
    : theta(theta_), omega(omega_)
  {
  }


  //! Dirichlet boundary condition type function
  template<typename IG>
  inline bool isDirichlet(const IG& intersection, const Dune::FieldVector<typename IG::ctype, IG::dimension-1>& xlocal) const
  {
    return false;
  }

  //! Neumann boundary condition type function
  template<typename IG>
  inline bool isNeumann(const IG& intersection, const Dune::FieldVector<typename IG::ctype, IG::dimension-1>& xlocal) const
  {
    return true;
  }

  //! Dirichlet boundary condition value
  RF g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal ) const
  {
      return RF(0.,0.);
  }

  //! Neumann boundary condition value
  template <typename IG>
  RF j (const IG& ig, const Dune::FieldVector<typename IG::ctype, IG::dimension-1>& xlocal, RF u) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> xglobal = ig.geometry().global(xlocal);
    RF i(0., 1.);
    RF x,y;
    x = xglobal[0];
    y = xglobal[1];
    // get the outer normal vector at the face center
    const Dune::FieldVector<RF, dim> normal = ig.centerUnitOuterNormal();
    Dune::FieldVector<RF, dim> gradu;
    gradu[0] = i*omega*cos(theta)*std::exp( i*omega*(x*cos(theta) + y *sin(theta)) );
    gradu[1] = i*omega*sin(theta)*std::exp( i*omega*(x*cos(theta) + y *sin(theta)) );
    RF res = gradu * normal;
    res *= -1;
    return res;
  }



  //! source term
  template <typename EG>
  RF f (const EG& eg,  const Dune::FieldVector<typename EG::Geometry::ctype, EG::Geometry::mydimension>& position) const
  {
    RF res(0.,0.);
    return res;
  }


  //! exact solution / analytic solution
  RF ua (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal ) const
  {

    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> xglobal = e.geometry().global(xlocal);

    RF i(0., 1.);
    RF x,y;
    x = xglobal[0];
    y = xglobal[1];

    return std::exp( i*omega*(x*cos(theta) + y *sin(theta)) );

  }

#if 0
  static const int dim = Traits::GridViewType::Grid::dimension;

  //! exact grad of solution / analytic grad of solution
  Dune::FieldVector<RF,dim> ugrada (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal ) const
  {

    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> xglobal = e.geometry().global(xlocal);
    RF i(0., 1.);
    RF x,y;
    x = xglobal[0];
    y = xglobal[1];
    Dune::FieldVector<RF,dim> res;
    res[0] = i*omega*cos(theta)*std::exp( i*omega*(x*cos(theta) + y *sin(theta)) );
    res[1] = i*omega*sin(theta)*std::exp( i*omega*(x*cos(theta) + y *sin(theta)) );

    return res;

  }
#endif


  const double theta;
  const double omega;

};
#endif
