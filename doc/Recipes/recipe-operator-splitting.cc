// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
// C and C++ includes
#include <iostream>
#include <sys/stat.h> // subfolder for VTK files
#include <cmath> // e.g. usage of pow(e,x)
#include <math.h>
#include <vector>
// dune-common includes
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/densematrix.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/timer.hh>
#include <dune/common/fvector.hh>
// dune-geometry includes
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>
// dune-grid includes
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/typetree/treepath.hh>

#include <dune/pdelab.hh>

/**
 * \page recipe-operator-splitting Operator splitting
 *
 * This recipe shows two implementation issues arising when operator
 * splitting approach is used: Communicating the data between two
 * systems of PDE to another, and calculating operator splitting errors.
 *
 * \section splitting-access-data Accessing data of another system from inside the LocalOperator
 *
 * We need a reference to the data vector and its GridFunctionSpace (GFS).
 * GFS tells us which parts of the vector to extract.
 * We store shared pointers to the data vector (type Data) and GFS
 * \snippet recipe-operator-splitting.cc Declare objects for data communication
 * The mutable token is necessary, because functions in the local oparator are const.
 *
 * They are instantiated in the constructor:
 * \snippet recipe-operator-splitting.cc Initialize objects for data communication
 *
 * We load data with respect to the local Entity. If we are accumulating Skeleton terms,
 * we need to repeat the process to get data from both inside and outside Entity.
 * \snippet recipe-operator-splitting.cc Read data
 * We got a vector with local data (from the other system). If the vector is a collection
 * of data from multiple variables, it places values of the same variable together. The data to
 * the first variable are at indices from 0 to lfs.size()-1, the second
 * variable occupies lfs.size() to 2*lfs.size()-1, etc.
 *
 * \section splitting-errors Accumulate splitting errors
 *
 * This example calculates a correction calculated by one operator splitting iteration.
 * Firstly, we need variables storing data (previous and current iteration results).
 * To extract a single data type from a coupled system use GridFunctionSubSpace.
 * The path to the data is defined through TreePath.
 * Create two DiscreteGridFunction-s which are fed to an Adapter (analytic function can be used in Adapter too),
 * and choose the Adapter -we chose DifferenceSquaredAdapter.
 * \snippet recipe-operator-splitting.cc Define objects for measuring splitting correction
 *
 * Then we use integrateGridFunction to get a vector with (local) results from Adapter,
 * \snippet recipe-operator-splitting.cc integrateGridFunction
 *
 * and sum the corrections. If the program is parallel, communicate.
 * \snippet recipe-operator-splitting.cc Global error
 * Note that communication passes a (const) reference and returns a value.
 * Some collective communication is necessary. Otherwise we risk one rank
 * continuing operator splitting iteration and other going to the next step
 * producing a deadlock.
 *
 * For a short summary of the communication, check @ref recipe-communication
 *
 * Full example code: @ref recipe-operator-splitting.cc
 * \example recipe-operator-splitting.cc
 * See explanation at @ref recipe-operator-splitting
 */


/**
 * Program solving water flow and contaminant transport in porous medium.
 * One phase unsaturated water flow is driven by Richards equation,
 * variable is hydraulic head and uses van Genuchten-Mualem model for
 * water retention curve and hydraulic permeability.
 * Contaminant part has two components C0, C1 representing concentrations.
 * They are transported via advection and diffusion (no dispersion).
 *
 * Whole system is split into two parts: water flow and contaminant transport.
 * Water flow does not depend on contaminants and can be solved separately.
 * The contaminant system is split by operator splitting technique into
 * two parts: advection and diffusion. Second order (in time) symmetric
 * Strang splitting is implemented. Firstly we take half step of convection,
 * then full diffusion step follower by second convection half step.
 * Operator splitting is iterative and we stop iterating when the correction
 * size is small.
 *
 * All three parts use finite volume elements on a regular 2D axiparallel
 * grid. Time stepping scheme is implicit Euler.
 */

template <typename Number, typename LType>
class Parameters
{
  Number time = 0.;
  Number finaltime = 1.; // reduced to get a reasonable testing run time
  Number dt = 0.2;
  Number othertime = 0.; // used in contaminant temporal part to find out which flow data should be loaded
  Number grav = -9.81; // negative, acts in a direction (0,-1)
  Number rho = 1.; // rescaled water density to get hydraulic head instead of pressure -> better conditioned system
  Number viscosity = 1.002e-3;
  Number pw_ini = -1.; // hydraulic head, pressure of 1 m tall water column
  // hydraulic properties of sand in van Genuchten-Mualem model, but not reflecting the rescaling to hydraulic head
  Number K_intristic = 1.76e-4;
  Number phi = 0.32; // porosity, usually 0.37 with 0.05 residual saturation, this simplified model uses zero residual saturations
  Number alpha = 3.5;
  Number n = 3.19;
  Number m; // = 1-1/n
  // some random coefficients for contaminant parts:
  Number C0_ini = 0.;
  Number C1_ini = 1.;
  Number D0 = 2e-6;
  Number D1 = 5e-7;
  const LType& L; // default: Dune::FieldVector<Number,2> L{0.1,0.1};
  Number inflow_ = 0.001;
  Number inflow_C0 = 1.;
  Number inflow_C1 = 0.;
public:
  using value_type = Number;

  Parameters (const LType& L_)
  : L(L_)
  {
    m = 1.-1./n;
  }

  void setTime (Number t_)
  {
    time = t_;
  }

  Number getTime () const
  {
    return time;
  }

  Number getT () const
  {
    return finaltime;
  }

  Number getTimeStep () const
  {
    return dt;
  }

  Number getWidth() const
  {
    return L[0];
  }

  Number getHeight() const
  {
    return L[1];
  }

  void setTimeStep (Number dt_)
  {
    dt = dt_;
  }

  void adjustTimeStep (Number coef)
  {
    dt *= coef;
  }

  void setOtherTime (Number t_)
  {
    othertime = t_;
  }

  Number getOtherTime () const
  {
    return othertime;
  }

  Number Sw (const Number& pw_) const
  {
    return (pw_ < 0.) ? 1./pow(1.+pow(-alpha*pw_,n),m) : 1.;
  }

  Number K (const Number& Sw_) const
  {
    // Sw_ is effective saturation, residual saturations are zero
    return K_intristic/viscosity*k_relative(Sw_);
  }

  Number k_relative (const Number& eSw) const
  {
    // uses effective saturation eSw \in (0,1)
    Number pom = 1. - pow( 1.-pow(eSw,1./m) ,m);
    Number result=sqrt(eSw)*pom*pom;
    return result;
  }

  Number getphi () const
  {
    return phi;
  }

  Number getrho () const
  {
    return rho;
  }

  Number getgrav () const
  {
    return grav;
  }

  Number getD0 (const Number& Sw_) const
  {
    // return D0*pow(phi,4./3.)*pow(Sw_,10./3.)
    return D0*phi*Sw_*Sw_*Sw_*cbrt(phi*Sw_);
  }

  Number getD1 (const Number& Sw_) const
  {
    // return D1*pow(phi,4./3.)*pow(Sw_,10./3.)
    return D1*phi*Sw_*Sw_*Sw_*cbrt(phi*Sw_);
  }

  template <typename E, typename X>
  Number g (const E& e, const X& x) const
  {
    auto global = e.geometry().global(x);
    // starting from a steady state
    return pw_ini+global[1]*grav*rho;
  }

  template <typename E, typename X>
  Number g_C0 (const E& e, const X& x) const
  {
    return C0_ini;
  }

  template <typename E, typename X>
  Number g_C1 (const E& e, const X& x) const
  {
    return C1_ini;
  }

  Number harg(Number a, Number b)
  {
    if (a<1e-20 || b<1e-20)
      return 0.;
    return 2.*a*b/(a+b);
  }

  template <typename I, typename X>
  bool b (const I& i, const X& x)
  {
    return false;
  }

  template <typename I, typename X>
  bool bC0 (const I& i, const X& x)
  {
    return false;
  }

  template <typename I, typename X>
  bool bC1 (const I& i, const X& x)
  {
    return false;
  }

  template <typename I, typename X>
  Number inflow (const I& i, const X& x)
  {
    return inflow_;
  }

  template <typename I, typename X>
  Number inflowC0 (const I& i, const X& x)
  {
    return inflow_C0;
  }

  template <typename I, typename X>
  Number inflowC1 (const I& i, const X& x)
  {
    return inflow_C1;
  }
};

class OpSplitException : public Dune::Exception
{};

/** a local operator for solving one-phase flow in porous medium
 * (Richards equation) using cell-centered finite volume method.
 * Van Genuchten-Mualem model is used for water retention curve and
 * hydraulic permeability. Variable is water pressure p_w, \phi is
 * porosity, \rho density, S_w saturation, \vec{g} gravity vector,
 * and K permeability.
 *
 * \f{align*}{
 *   \partial_t (\phi S_w) + \mathrm{div} \left( -K(S_w) (\nabla p_w - \vec{g} \rho) \right) &=& 0 \ x \in \Omega \\
 *   pw                                    &=& pw_ini     \ t=0 \\
 *   pw \cdot \nu                          &=& 0          \ x \in \Gamma_N \\
 *   K (\nabla pw - \vec{g}\rho) \cdot \nu &=& q_{inflow} \ x \in \Gamma_{inflow} \\
 *   pw \leq 0, K(\nabla pw -\vec{g}\rho) \cdot \nu &\leq& 0, pw \left(K(\nabla pw -\vec{g}\rho) \cdot \nu\right) = 0 \ x \in \Gamma_{free outflow}
 * \f}
 *
 */
template<typename Param>
class LOP_spatial_flow
  : public Dune::PDELab::FullVolumePattern
  , public Dune::PDELab::FullSkeletonPattern
  , public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::NumericalJacobianVolume<LOP_spatial_flow<Param> >
  , public Dune::PDELab::NumericalJacobianSkeleton<LOP_spatial_flow<Param> >
  , public Dune::PDELab::NumericalJacobianBoundary<LOP_spatial_flow<Param> >
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  Param& param;
  using RF = typename Param::value_type;

public:
  // pattern assembly flags
  enum { doPatternVolume = true };
  enum { doPatternSkeleton = true };

  // residual assembly flags
  // enum { doLambdaBoundary = true };// side boundaries in alpha
  // enum { doAlphaVolume   = true }; // terms from integrals over cells
  enum { doAlphaSkeleton = true }; // terms from integrals over edges between cells
  enum { doAlphaBoundary = true }; // terms from integrals over boundaries

  //! constructor stores a copy of the parameter object
  LOP_spatial_flow (Param& param_)
    : param(param_)
  {}

  void setTime(RF t)
  {
    param.setTime(t);
  }

  //! residual contribution from skeleton terms
  template<typename IG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_skeleton (const IG& ig,
         const LFSU& lfsu_i, const X& x_i, const LFSV& lfsv_i,
         const LFSU& lfsu_o, const X& x_o, const LFSV& lfsv_o,
         R& r_i, R& r_o) const
  {
    // inside and outside cells
    auto cell_inside = ig.inside();
    auto cell_outside = ig.outside();

    // inside and outside cell geometries
    auto insidegeo = cell_inside.geometry();
    auto outsidegeo = cell_outside.geometry();

    // cell centers in global coordinates
    auto inside_global = insidegeo.center();
    auto outside_global = outsidegeo.center();

    // distance between the two cell centers
    inside_global -= outside_global;
    auto distance = inside_global.two_norm();

    // face volume for integration
    auto facegeo = ig.geometry();
    auto volume = facegeo.volume();

    // evaluations for gradient, gravity, and permeability
    auto pw_i = x_i(lfsu_i,0);
    auto Sw_i = param.Sw(pw_i);
    auto pw_o = x_o(lfsu_o,0);
    auto Sw_o = param.Sw(pw_o);
    auto dpwdn = (pw_o-pw_i)/distance; // \nabla pw \cdot\nu
    auto g = param.getgrav()*param.getrho()*inside_global[1]/distance; // -\vec{g}\cdot\vec{\nu}
    auto K_i = param.K(Sw_i);
    auto K_o = param.K(Sw_o);
    auto K = param.harg(K_i,K_o);

    // contribution to residual on inside and outside elements
    auto q = -K*(dpwdn+g);
    r_i.accumulate(lfsv_i,0, q*volume);
    r_o.accumulate(lfsv_o,0,-q*volume);
  }

  template<typename IG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_boundary (const IG& ig,
                       const LFSU& lfsu_i, const X& x_i,
                       const LFSV& lfsv_i, R& r_i) const
  {
    const auto& volume = ig.geometry().volume();
    const auto& global = ig.geometry().center();

    if (global[0]<=1e-6)
    {
      if ( (global[1]<(0.9*param.getHeight()) ) && (global[1]>(0.5*param.getHeight()) ) )
      {
        r_i.accumulate(lfsv_i,0,-param.inflow(ig,x_i)*volume);
      }
    }
    else if (global[0]>=param.getWidth()-1e-6)
    {
      auto pw = x_i(lfsu_i,0);
      auto K = param.K(param.Sw(pw));
      auto outflow = -K*(0.-pw);
      r_i.accumulate(lfsv_i,0,std::max(0.,outflow)*volume);
    }

  } // end alpha_boundary

}; // end LOP_spatial_flow

template<typename Param>
class LOP_time_flow
  : public Dune::PDELab::FullVolumePattern
  , public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::NumericalJacobianVolume<LOP_time_flow<Param> >
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  Param& param;
  using RF = typename Param::value_type;
public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };

  LOP_time_flow(Param& param_)
    : param(param_)
  {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
                     const LFSV& lfsv, R& r) const
  {
    // types & dimension
    // const int dim = EG::Entity::dimension;
    const auto volume = eg.geometry().volume();

    auto pw  = x(lfsu,0);
    auto Sw  = param.Sw(pw);
    auto phi = param.getphi();

    r.accumulate(lfsv,0, phi*Sw*volume );
  }
}; // end LOP_time_contaminant


/** a local operator for solving contaminant transport in porous
 *  medium using cell-centered finite volume method.
 *  Variables are molar concentrations in water C0, and C1.
 *  \phi is porosity, \vec{q_w} is water flux, D_* diffusion coefficient (no dispersion)
 *
 * \f{align*}{
 *   \partial_t(\phi S_w C0) + \mathrm{div}[ \frac{\vec{q_w}}{\rho_w} C0 - D_1\nabla C0 ] &=& 0 \ x \in \Omega,\\
 *   \partial_t(\phi S_w C1) + \mathrm{div}[ \frac{\vec{q_w}}{\rho_w} C1 - D_2\nabla C1 ] &=& 0 \ x \in \Omega,\\
 *   C0 = C1_{ini}, C1 &=& C2_{ini}                               \ t=0 \\
 *   (\frac{\vec{q_w}}{\rho_w} C0 - D_1\nabla C0) \cdot \nu &=& 0 \ x \in \Gamma_N \\
 *   (\frac{\vec{q_w}}{\rho_w} C1 - D_2\nabla C1) \cdot \nu &=& 0 \ x \in \Gamma_N \\
 *   C0 = C0_{inflow}, C1 &=& C1_{inflow}                         \ x \in \Gamma_{inflow} \\
 *   \nabla C0 \cdot \nu = 0, \nabla C1 \cdot \nu &=& 0,          \ x \in \Gamma_{free outflow}
 * \f}
 *
 */

template<typename Param, typename GFSF, typename ZF, typename GFSC, typename ZC, bool convection>
// GFSF,ZF are of water flow type
class LOP_spatial_contaminant
  : public Dune::PDELab::FullVolumePattern
  , public Dune::PDELab::FullSkeletonPattern
  , public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::NumericalJacobianVolume<LOP_spatial_contaminant<Param,GFSF,ZF,GFSC,ZC,convection> >
  , public Dune::PDELab::NumericalJacobianSkeleton<LOP_spatial_contaminant<Param,GFSF,ZF,GFSC,ZC,convection> >
  , public Dune::PDELab::NumericalJacobianBoundary<LOP_spatial_contaminant<Param,GFSF,ZF,GFSC,ZC,convection> >
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  Param& param;
  using RF = typename Param::value_type;
  // [Declare objects for data communication]
  using LFSF      = Dune::PDELab::LocalFunctionSpace<GFSF>;
  using LFSFCache = Dune::PDELab::LFSIndexCache<LFSF>;
  std::shared_ptr<ZF> dataf;
  std::shared_ptr<const GFSF> pgfsf;
  mutable LFSF lfsf;
  mutable LFSFCache lfsf_cache;
  //! [Declare objects for data communication]
  using LFSC      = Dune::PDELab::LocalFunctionSpace<GFSC>;
  using LFSCCache = Dune::PDELab::LFSIndexCache<LFSC>;
  std::shared_ptr<ZC> datac;
  std::shared_ptr<const GFSC> pgfsc;
  mutable LFSC lfsc;
  mutable LFSCCache lfsc_cache;

public:
  // pattern assembly flags
  enum { doPatternVolume = true };
  enum { doPatternSkeleton = true };

  // residual assembly flags
  // enum { doLambdaBoundary = true }; // side boundaries in alpha
  // enum { doAlphaVolume   = true }; // reaction and source terms
  enum { doAlphaSkeleton = true }; // flow between cells
  enum { doAlphaBoundary = true }; // boundaries

  // [Initialize objects for data communication]
  LOP_spatial_contaminant (Param&  param_, const GFSF& gfsf_, ZF& zf_, const GFSC& gfsc_, ZC& zc_)
    : param(param_)
    , dataf(stackobject_to_shared_ptr(zf_))
    , pgfsf(stackobject_to_shared_ptr(gfsf_))
    , lfsf(pgfsf)
    , lfsf_cache(lfsf)
  //! [Initialize objects for data communication]
    , datac(stackobject_to_shared_ptr(zc_))
    , pgfsc(stackobject_to_shared_ptr(gfsc_))
    , lfsc(pgfsc)
    , lfsc_cache(lfsc)
  {}

  void setTime(double t)
  {
    param.setTime(t);
  }

  void adjustTimeStep(double coef)
  {
    param.adjustTimeStep(coef);
  }

  //! residual contribution from skeleton terms
  template<typename IG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_skeleton (const IG& ig,
         const LFSU& lfsu_i, const X& x_i, const LFSV& lfsv_i,
         const LFSU& lfsu_o, const X& x_o, const LFSV& lfsv_o,
         R& r_i, R& r_o) const
  {
    using namespace Dune::TypeTree::Indices;
    auto lfsv_i0 = lfsv_i.child(_0);
    auto lfsv_i1 = lfsv_i.child(_1);
    auto lfsu_i0 = lfsu_i.child(_0);
    auto lfsu_i1 = lfsu_i.child(_1);
    auto lfsv_o0 = lfsv_o.child(_0);
    auto lfsv_o1 = lfsv_o.child(_1);
    auto lfsu_o0 = lfsu_o.child(_0);
    auto lfsu_o1 = lfsu_o.child(_1);
    // inside and outside cells
    // [Read data]
    auto cell_inside = ig.inside(); // entity object
    typename ZF::template LocalView<LFSFCache> p_view(*dataf);
    lfsf.bind(cell_inside);
    lfsf_cache.update();
    std::vector<RF> pw(lfsf.size());
    p_view.bind(lfsf_cache);
    p_view.read(pw);
    p_view.unbind();
    //! [Read data]
    auto cell_outside = ig.outside();

    // inside and outside cell geometries
    auto insidegeo = cell_inside.geometry();
    auto outsidegeo = cell_outside.geometry();

    // cell centers in global coordinates
    auto inside_global = insidegeo.center();
    auto outside_global = outsidegeo.center();

    // distance between the two cell centers
    inside_global -= outside_global;
    auto distance = inside_global.two_norm();

    // face volume for integration
    auto facegeo = ig.geometry();
    auto volume = facegeo.volume();
    auto pw_i = pw[0];
    auto Sw_i = param.Sw(pw_i);
    lfsf.bind(cell_outside);
    lfsf_cache.update();
    p_view.bind(lfsf_cache);
    p_view.read(pw);
    p_view.unbind();
    auto pw_o = pw[0];
    auto Sw_o = param.Sw(pw_o);
    typename ZC::template LocalView<LFSCCache> c_view(*datac);
    lfsc.bind(cell_inside);
    lfsc_cache.update();
    std::vector<RF> C_i(lfsc.size());
    c_view.bind(lfsc_cache);
    c_view.read(C_i);
    c_view.unbind();
    lfsc.bind(cell_outside);
    lfsc_cache.update();
    std::vector<RF> C_o(lfsc.size());
    c_view.bind(lfsc_cache);
    c_view.read(C_o);
    c_view.unbind();

    // data used for convection
    auto C0c_i = convection ? x_i(lfsu_i0,0) : C_i[0];
    auto C1c_i = convection ? x_i(lfsu_i1,0) : C_i[1];
    auto C0c_o = convection ? x_o(lfsu_o0,0) : C_o[0];
    auto C1c_o = convection ? x_o(lfsu_o1,0) : C_o[1];
    // data used for diffusion
    auto C0d_i = convection ? C_i[0] : x_i(lfsu_i0,0);
    auto C1d_i = convection ? C_i[1] : x_i(lfsu_i1,0);
    auto C0d_o = convection ? C_o[0] : x_o(lfsu_o0,0);
    auto C1d_o = convection ? C_o[1] : x_o(lfsu_o1,0);

    // water flow
    auto dpwdn = (pw_o-pw_i)/distance; // \nabla pw \cdot\nu
    auto g = param.getgrav()*param.getrho()*inside_global[1]/distance; // -\vec{g}\cdot\vec{\nu}
    auto K_i = param.K(Sw_i);
    auto K_o = param.K(Sw_o);
    auto K = param.harg(K_i,K_o);
    auto q = -K*(dpwdn+g);

    auto dC0dn = (C0d_o-C0d_i)/distance;
    auto dC1dn = (C1d_o-C1d_i)/distance;

    // convection
    // simple upwind
    decltype(C0c_i) C0;
    decltype(C1c_i) C1;
    if ( q > 0 )
    {
      C0 = C0c_i;
      C1 = C1c_i;
    }
    else
    {
      C0 = C0c_o;
      C1 = C1c_o;
    }
    r_i.accumulate(lfsv_i0,0, q*C0*volume);
    r_o.accumulate(lfsv_o0,0,-q*C0*volume);
    r_i.accumulate(lfsv_i1,0, q*C1*volume);
    r_o.accumulate(lfsv_o1,0,-q*C1*volume);

    // diffusion
    auto D0_i = param.getD0(Sw_i);
    auto D0_o = param.getD0(Sw_o);
    auto D0 = param.harg(D0_i,D0_o);
    auto D1_i = param.getD1(Sw_i);
    auto D1_o = param.getD1(Sw_o);
    auto D1 = param.harg(D1_i,D1_o);
    r_i.accumulate(lfsv_i0,0,-D0*dC0dn*volume);
    r_o.accumulate(lfsv_o0,0, D0*dC0dn*volume);
    r_i.accumulate(lfsv_i1,0,-D1*dC1dn*volume);
    r_o.accumulate(lfsv_o1,0, D1*dC1dn*volume);
  }

  template<typename IG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_boundary (const IG& ig,
                       const LFSU& lfsu_i, const X& x_i,
                       const LFSV& lfsv_i, R& r_i) const
  {
    using namespace Dune::TypeTree::Indices;
    auto lfsv_i0 = lfsv_i.child(_0);
    auto lfsv_i1 = lfsv_i.child(_1);
    auto lfsu_i0 = lfsu_i.child(_0);
    auto lfsu_i1 = lfsu_i.child(_1);
    const auto& volume = ig.geometry().volume();
    const auto& global = ig.geometry().center();
    const auto& inside_cell = ig.inside();

    if (global[0]<=1e-6)
    {
      if (global[1]<0.9*param.getHeight() && global[1]>0.5*param.getHeight())
      {
        // only convective term present at the inflow boundary
        r_i.accumulate(lfsv_i0,0,-param.inflow(ig,x_i)*param.inflowC0(ig,x_i)*volume);
        r_i.accumulate(lfsv_i1,0,-param.inflow(ig,x_i)*param.inflowC1(ig,x_i)*volume);
      }
    }
    else if (global[0]>=param.getWidth()-1e-6)
    {
      typename ZF::template LocalView<LFSFCache> p_view(*dataf);
      lfsf.bind(inside_cell);
      lfsf_cache.update();
      std::vector<RF> pw(lfsf.size());
      p_view.bind(lfsf_cache);
      p_view.read(pw);
      p_view.unbind();
      auto K = param.K(param.Sw(pw[0]));
      auto outflow = -K*(0.-pw[0]);
      auto C0 = x_i(lfsv_i0,0);
      auto C1 = x_i(lfsv_i1,0);
      r_i.accumulate(lfsv_i0,0,std::max(0.,outflow)*C0*volume);
      r_i.accumulate(lfsv_i1,0,std::max(0.,outflow)*C1*volume);
    }

  } // end alpha_boundary

}; // end LOP_spatial_contaminant

template<typename Param, typename GFSF, typename ZF, bool convection>
class LOP_time_contaminant
  : public Dune::PDELab::FullVolumePattern
  , public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::NumericalJacobianVolume<LOP_time_contaminant<Param,GFSF,ZF,convection> >
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  Param& param;
  using RF = typename Param::value_type;
  using LFSF = Dune::PDELab::LocalFunctionSpace<GFSF>;
  using LFSFCache = Dune::PDELab::LFSIndexCache<LFSF>;
  std::shared_ptr<ZF> dataf;
  std::shared_ptr<ZF> datafold;
  std::shared_ptr<const GFSF> pgfsf;
  mutable LFSF lfsf;
  mutable LFSFCache lfsf_cache;
public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };

  LOP_time_contaminant(Param& param_, const GFSF& gfsf_, ZF& zf_, ZF& zfold_)
    : param(param_)
    , dataf(stackobject_to_shared_ptr(zf_))
    , datafold(stackobject_to_shared_ptr(zfold_))
    , pgfsf(stackobject_to_shared_ptr(gfsf_))
    , lfsf(pgfsf)
    , lfsf_cache(lfsf)
  {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
                     const LFSV& lfsv, R& r) const
  {
    // get access to each vector space
    using namespace Dune::TypeTree::Indices;
    const auto& lfsu_0 = lfsu.child(_0);
    const auto& lfsu_1 = lfsu.child(_1);
    const auto& lfsv_0 = lfsv.child(_0);
    const auto& lfsv_1 = lfsv.child(_1);
    const auto& inside_cell = eg.entity();

    // determine which term of the time derivative we are evaluating
    bool new_or_old = param.getTime() > param.getOtherTime();

    // load flow data
    typename ZF::template LocalView<LFSFCache> p_view(new_or_old ? *dataf : *datafold);
    lfsf.bind(inside_cell);
    lfsf_cache.update();
    std::vector<RF> pw(lfsf.size());
    p_view.bind(lfsf_cache);
    p_view.read(pw);
    p_view.unbind();

    auto Sw = param.Sw(pw[0]);
    auto C0 = x(lfsu_0,0);
    auto C1 = x(lfsu_1,0);
    auto phi = param.getphi();

    const auto& volume = eg.geometry().volume();
    r.accumulate(lfsv_0,0, phi*Sw*C0*volume);
    r.accumulate(lfsv_1,0, phi*Sw*C1*volume);
  }
}; // end LOP_time_contaminant


template<typename GV, typename FEM, typename Param>
void driver (const GV& gv, const FEM& fem, Param& param)
{
  // type for computations
  using RF = typename Param::value_type;

  // Make grid function space
  using CON = Dune::PDELab::P0ParallelConstraints;
  using VBE = Dune::PDELab::ISTL::VectorBackend<>;
  using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE>;
  GFS gfs(gv,fem); // GridFuntionSpace for flow
  gfs.name("pressure");

  using VBES = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
  using OrderingTag = Dune::PDELab::EntityBlockedOrderingTag;
  using GFSC = Dune::PDELab::PowerGridFunctionSpace<GFS,2,VBES,OrderingTag>;
  GFSC gfsc(gfs); // GridFunctionSpace for contaminant
  // gfsc names for VTK output
  using namespace Dune::TypeTree::Indices;
  gfsc.child(_0).name("C0");
  gfsc.child(_1).name("C1");

  // Coefficient vectors
  using Z = Dune::PDELab::Backend::Vector<GFS, RF>;
  Z z(gfs);
  using ZC = Dune::PDELab::Backend::Vector<GFSC, RF>;
  ZC zc(gfsc);

  // vectors for data manipulation
  Z  zfnew(z);  // result from apply() for flow
  ZC zcnew(zc); // result from apply() for contaminant
  Z  zfdata(z); // flow data accessed from the contaminant systems
  ZC zcdata(zc); // contaminant data accessed from the other contaminant system
  ZC zcstep(zc); // previous splitting iteration data, used to compute corrections
  ZC zchalfstep(zc); // contaminant data from convection part after first dt/2 step in op-split iteration

  // initial conditions
  RF time = 0.0;
  auto glambdaf = [&](const auto& e, const auto& x)
  {
    RF gf;
    gf = param.g(e,x);
    return gf;
  };
  auto glambdac = [&] (const auto& e, const auto& x)
  {
    Dune::FieldVector<RF,2> gc;
    gc[0] = param.g_C0(e,x);
    gc[1] = param.g_C1(e,x);
    return gc;
  };
  auto gf = Dune::PDELab::
    makeInstationaryGridFunctionFromCallable(gv,glambdaf,param);
  auto gc = Dune::PDELab::
    makeInstationaryGridFunctionFromCallable(gv,glambdac,param);

  // initialize coefficient vectors
  Dune::PDELab::interpolate(gf,gfs,z);
  Dune::PDELab::interpolate(gc,gfsc,zc);

  // boundary conditions: finite volume does not use any, but parallel solvers need them for overlap
  using CF = typename GFS::template ConstraintsContainer<RF>::Type;
  using CC = typename GFSC::template ConstraintsContainer<RF>::Type;
  CF cf;
  CC cc;
  Dune::PDELab::constraints(gfs,cf);
  Dune::PDELab::constraints(gfsc,cc);

  // Make instationary grid operator for flow
  using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  MBE mbe((int)10);
  using LOPF = LOP_spatial_flow<Param>;
  LOPF lopf(param);
  using GOF0 = Dune::PDELab::GridOperator<GFS,GFS,LOPF,MBE,RF,RF,RF,CF,CF>;
  GOF0 gof0(gfs,cf,gfs,cf,lopf,mbe);
  using TLOPF = LOP_time_flow<Param>;
  TLOPF tlopf(param);
  using GOF1 = Dune::PDELab::GridOperator<GFS,GFS,TLOPF,MBE,RF,RF,RF,CF,CF>;
  GOF1 gof1(gfs,cf,gfs,cf,tlopf,mbe);
  using IGOF = Dune::PDELab::OneStepGridOperator<GOF0,GOF1>;
  IGOF igof(gof0,gof1);
  // same for contaminant transport
  constexpr bool convection{true};
  using LOPC = LOP_spatial_contaminant<Param,GFS,Z,GFSC,ZC,convection>;
  LOPC lopc(param,gfs,zfdata,gfsc,zcdata);
  using GOC0 = Dune::PDELab::GridOperator<GFSC,GFSC,LOPC,MBE,RF,RF,RF,CC,CC>;
  GOC0 goc0(gfsc,cc,gfsc,cc,lopc,mbe);
  using TLOPC = LOP_time_contaminant<Param,GFS,Z,convection>;
  TLOPC tlopc(param,gfs,zfdata,z);
  using GOC1 = Dune::PDELab::GridOperator<GFSC,GFSC,TLOPC,MBE,RF,RF,RF,CC,CC>;
  GOC1 goc1(gfsc,cc,gfsc,cc,tlopc,mbe);
  using IGOC = Dune::PDELab::OneStepGridOperator<GOC0,GOC1>;
  IGOC igoc(goc0,goc1);
  // same for contaminant transport
  constexpr bool diffusion{false};
  using LOPD = LOP_spatial_contaminant<Param,GFS,Z,GFSC,ZC,diffusion>;
  LOPD lopd(param,gfs,zfdata,gfsc,zcdata);
  using GOD0 = Dune::PDELab::GridOperator<GFSC,GFSC,LOPD,MBE,RF,RF,RF,CC,CC>;
  GOD0 god0(gfsc,cc,gfsc,cc,lopd,mbe);
  using TLOPD = LOP_time_contaminant<Param,GFS,Z,diffusion>;
  TLOPD tlopd(param,gfs,zfdata,z);
  using GOD1 = Dune::PDELab::GridOperator<GFSC,GFSC,TLOPD,MBE,RF,RF,RF,CC,CC>;
  GOD1 god1(gfsc,cc,gfsc,cc,tlopd,mbe);
  using IGOD = Dune::PDELab::OneStepGridOperator<GOD0,GOD1>;
  IGOD igod(god0,god1);

  // Select a parallel linear solver backend
  using LSF = Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,CF>;
  using LSC = Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFSC,CC>;
  int verbose=0;
  if (gfs.gridView().comm().rank()==0) verbose=1;
  LSF lsf(gfs,cf,100,5,verbose);
  LSC lsc(gfsc,cc,100,5,verbose);

  // select a nonlinear solver for flow
  constexpr int newtonMaxIt{10};
  using PDESOLVERF = Dune::PDELab::NewtonMethod<IGOF,LSF>;
  PDESOLVERF pdesolverf(igof,lsf);
  Dune::ParameterTree newtonParam;
  newtonParam["ReassembleThreshold"] = "0.0";
  newtonParam["VerbosityLevel"] = "2";
  newtonParam["Reduction"] = "1e-8";
  newtonParam["MinLinearReduction"] = "1e-4";
  newtonParam["MaxIterations"] = std::to_string(newtonMaxIt);
  newtonParam["LineSearchMaxIterations"] = "10";
  pdesolverf.setParameters(newtonParam);
  // select a solver for contaminant
  using PDESOLVERC = Dune::PDELab::NewtonMethod<IGOC,LSC>;
  PDESOLVERC pdesolverc(igoc,lsc);
  pdesolverc.setParameters(newtonParam);
  using PDESOLVERD = Dune::PDELab::NewtonMethod<IGOD,LSC>;
  PDESOLVERD pdesolverd(igod,lsc);
  pdesolverd.setParameters(newtonParam);

  // select and prepare time-stepping scheme
  Dune::PDELab::OneStepThetaParameter<RF> method(1.0);
  Dune::PDELab::TimeSteppingParameterInterface<RF>* pmethod=&method;
  using OSMF = Dune::PDELab::OneStepMethod<RF,IGOF,PDESOLVERF,Z,Z>;
  OSMF  osmf(*pmethod,igof,pdesolverf);
  osmf.setVerbosityLevel(2);
  using OSMC = Dune::PDELab::OneStepMethod<RF,IGOC,PDESOLVERC,ZC,ZC>;
  OSMC  osmc(*pmethod,igoc,pdesolverc);
  osmc.setVerbosityLevel(2);
  using OSMD = Dune::PDELab::OneStepMethod<RF,IGOD,PDESOLVERD,ZC,ZC>;
  OSMD  osmd(*pmethod,igod,pdesolverd);
  osmd.setVerbosityLevel(2);

  // [Define objects for measuring splitting correction]
  using GFSC_Sub0 = Dune::PDELab::GridFunctionSubSpace<GFSC,Dune::TypeTree::StaticTreePath<0>>;
  GFSC_Sub0 gfsc_sub0(gfsc);
  using SubC0 = Dune::PDELab::DiscreteGridFunction<GFSC_Sub0,ZC>;
  SubC0 subC0(gfsc_sub0,zcstep);
  SubC0 subC0new(gfsc_sub0,zcnew);
  using ErrC0 = Dune::PDELab::DifferenceSquaredAdapter<SubC0,SubC0>;
  ErrC0 errC0(subC0,subC0new);
  //! [Define objects for measuring splitting correction]
  using GFSC_Sub1 = Dune::PDELab::GridFunctionSubSpace<GFSC,Dune::TypeTree::StaticTreePath<1>>;
  GFSC_Sub1 gfsc_sub1(gfsc);
  using SubC1 = Dune::PDELab::DiscreteGridFunction<GFSC_Sub1,ZC>;
  SubC1 subC1(gfsc_sub1,zcstep);
  SubC1 subC1new(gfsc_sub1,zcnew);
  using ErrC1 = Dune::PDELab::DifferenceSquaredAdapter<SubC1,SubC1>;
  ErrC1 errC1(subC1,subC1new);

  // prepare VTK writer and write first file
  int subsampling=1;
  using VTKWRITER = Dune::SubsamplingVTKWriter<GV>;
  VTKWRITER vtkwriter(gv,Dune::refinementIntervals(subsampling));
  std::string filename="recipe-operator-splitting_parallel";
  struct stat stf;
  if( stat( filename.c_str(), &stf ) != 0 )
  {
    int statf = 0;
    statf = mkdir(filename.c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
    if( statf != 0 && statf != -1)
      std::cout << "Error: Cannot create directory "
                << filename << std::endl;
  }
  using VTKSEQUENCEWRITER = Dune::VTKSequenceWriter<GV>;
  VTKSEQUENCEWRITER vtkSequenceWriter(
    std::make_shared<VTKWRITER>(vtkwriter),filename,filename,"");
  // add data field(s) for all components of the space to the VTK writer
  Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter,gfs,z);
  Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter,gfsc,zc);
  vtkSequenceWriter.write(0.0,Dune::VTK::appendedraw);

  // time loop --------------------------------------------------------------------
  RF T = param.getT();
  param.setTime(time);
  param.setOtherTime(time); // used in contaminant temporal term to load the correct flow data
  RF lastvtk = 0.; // time of the last vtk file, used to avoid creating too many frames

  while (time<T-1e-8)
  {
    try
    {
      unsigned int iterations{0};
      unsigned int split_iter{0};
      bool continue_opsplit{true};
      if (verbose>=1)
        std::cout << std::endl << std::endl << "flow part:\n";
      osmf.apply(time,param.getTimeStep(),z,gf,zfnew);
      zfdata = zfnew;
      iterations=pdesolverf.result().iterations;

      // iterative operator splitting for contaminant
      do
      {
        if (verbose>=1)
          std::cout << std::endl << std::endl << "contaminant convection part 1/2:\n";
        osmc.apply(time,param.getTimeStep()/2,zc,gc,zcnew);
        zcdata = zcnew;
        zchalfstep = zcnew;
        iterations=std::max(iterations,pdesolverc.result().iterations);

        if (verbose>=1)
          std::cout << std::endl << "contaminant diffusion part:\n";
        osmd.apply(time,param.getTimeStep(),zc,gc,zcnew);
        zcdata = zcnew;
        iterations=std::max(iterations,pdesolverd.result().iterations);

        if (verbose>=1)
          std::cout << std::endl << "contaminant convection part 2/2:\n";
        osmc.apply(time+param.getTimeStep()/2,param.getTimeStep()/2,zchalfstep,gc,zcnew);
        zcdata = zcnew;
        iterations=std::max(iterations,pdesolverc.result().iterations);

        // [integrateGridFunction]
        typename ErrC0::Traits::RangeType spliterrorcont0;
        Dune::PDELab::integrateGridFunction(errC0,spliterrorcont0,1);
        //! [integrateGridFunction]
        typename ErrC1::Traits::RangeType spliterrorcont1;
        Dune::PDELab::integrateGridFunction(errC1,spliterrorcont1,1);
        // [Global error]
        RF sperrC0 = spliterrorcont0.one_norm();
        sperrC0 = gv.comm().sum(sperrC0);
        //! [Global error]
        RF sperrC1 = spliterrorcont1.one_norm();
        sperrC1 = gv.comm().sum(sperrC1);

        if (sperrC0 < 1e-17 && sperrC1 < 1e-19 && split_iter > 0)
        {
          if (verbose>=1)
            std::cout << "low spliterrors (" << sperrC0 << ", " << sperrC1 << "), go to next timestep" << std::endl;
          continue_opsplit = false;
        }
        else
        {
          if (split_iter >= 5)
          {
            if (verbose>=1)
              std::cout << "max operator splitting iteration number reached, current errors: " << sperrC0 << ", " << sperrC1 << ", reseting time step" << std::endl;
            continue_opsplit = true;
            DUNE_THROW(OpSplitException,"Operator splitting did not converge");
          }
          else
          {
            if (verbose>=1)
              std::cout << "going to " << split_iter+1 << ". operator splitting iteration, errors: " << sperrC0 << ", " << sperrC1 << std::endl;
            continue_opsplit = true;
          }
        }
        zcstep = zcnew;
        zcdata = zcnew;
        ++split_iter;
      } while (continue_opsplit);
      zc=zcstep;
      z=zfnew;
      time+=param.getTimeStep();
      param.setOtherTime(time);

      // write solution
      if (time-lastvtk > time*0.01)
      {
        vtkSequenceWriter.write(time,Dune::VTK::appendedraw);
        lastvtk = time;
      }

      if (iterations<(newtonMaxIt*4)/10)
      {
        param.adjustTimeStep(1.25);
      }
      else
        if ((iterations>(newtonMaxIt*7)/10) && param.getTimeStep()>1e-4)
        {
          iterations = 0;
          param.adjustTimeStep(0.75);
        }
    } // end try
    catch(Dune::Exception& e)
    {
      if (verbose==1)
      {
        std::cout << e.what() << std::endl;
      }
      if (param.getTimeStep()<1e-4)
      {
        if (verbose==1)
        {
          std::cout << "too little time step" << std::endl;
        }
        vtkSequenceWriter.write(time+T,Dune::VTK::appendedraw);
        throw;
      }
      else
      {
        param.adjustTimeStep(0.75);
        zfnew = z;
        zcnew = zc;
        if (verbose==1)
        {
          std::cout << "Reducing TimeStep size to " << param.getTimeStep() << std::endl;
        }
      }
    }
    catch(std::exception& e)
    {
      if (verbose==1)
      {
        std::cout << e.what() << std::endl;
      }
      if (param.getTimeStep()<1e-4)
      {
        if (verbose==1)
        {
          std::cout << "too little time step" << std::endl;
        }
        vtkSequenceWriter.write(time+T,Dune::VTK::appendedraw);
        throw;
      }
      else
      {
        param.adjustTimeStep(0.75);
        zfnew = z;
        zcnew = zc;
        if (verbose==1)
        {
          std::cout << "Reducing TimeStep size to " << param.getTimeStep() << std::endl;
        }
      }
    }
    catch(...) // in case of unexpected exception terminate
    {
      throw;
    }
  } // end while
} // end driver


template <int dim>
class YaspPartition : public Dune::YLoadBalance<dim>
{
public:
  using iTuple = std::array<int,dim>;
  void loadbalance (const iTuple& size, int P, iTuple& dims) const
  {
    // greedy algorithm for splitting P into product of dim numbers, closest to \sqrt[dim]{P}
    for (int j=0; j<dim; ++j)
    {
      int nextP = P;
      for (int i=1; pow(i,dim-j)<=P; ++i)
      {
        if(P%i == 0)
        {
          dims[j] = i;
          nextP = P/i;
        }
      }
      P = nextP;
    }
    auto greater = [](int a,int b) -> bool {return a>b;};
    std::sort(dims.begin(),dims.end(),greater);
  }
};

int main(int argc, char** argv)
{
  try{
    // Maybe initialize MPI
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
      std::cout<<"I am rank "<<helper.rank()<<" of "<<helper.size()
        <<" processes!"<<std::endl;
    constexpr int dim = 2;
    using RF = double;
    using Grid = Dune::YaspGrid<dim>;
    using DF = Grid::ctype;
    using LType = Dune::FieldVector<DF,dim>;
    LType L;
    L[0] = 0.1;
    L[1] = 0.1;
    std::array<int,dim> N;
    N[0] = 16;
    N[1] = 16;
    std::bitset<dim> periodic(false);
    int overlap=1;
    int refinement = 1;
    YaspPartition<dim> yp;
    std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(new Grid(L,N,periodic,overlap,Dune::MPIHelper::getCollectiveCommunication(),&yp));
    gridp->refineOptions(false); // keep overlap in cells
    gridp->globalRefine(refinement);
    using GV = Grid::LeafGridView;
    GV gv=gridp->leafGridView();
    using FEM = Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim>;
    FEM fem(Dune::GeometryTypes::cube(dim));
    using Param = Parameters<RF,LType>;
    Param param(L);
    driver(gv,fem,param);
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
