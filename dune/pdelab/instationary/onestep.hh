// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ONESTEP_HH
#define DUNE_PDELAB_ONESTEP_HH

#include <iomanip>
#include <iostream>
#include <ostream>
#include <vector>

#include <stdio.h>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/ios_state.hh>

#include <dune/pdelab/common/logtag.hh>
#include <dune/pdelab/gridoperator/common/timesteppingparameterinterface.hh>

namespace Dune {
  namespace PDELab {

    /**
     *  @addtogroup OneStepMethod
     *  @{
     */
    /**
     * \brief Parameters to turn the ExplicitOneStepMethod into an
     * explicite Euler method.
     *
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class ExplicitEulerParameter : public TimeSteppingParameterInterface<R>
    {
    public:

      ExplicitEulerParameter ()
      {
        D[0] = 0.0;  D[1] = 1.0;
        A[0][0] = -1.0; A[0][1] = 1.0;
        B[0][0] = 1.0;  B[0][1] = 0.0;
      }

      /*! \brief Return true if method is implicit
       */
      virtual bool implicit () const
      {
        return false;
      }

      /*! \brief Return number of stages s of the method
       */
      virtual unsigned s () const
      {
        return 1;
      }

      /*! \brief Return entries of the A matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R a (int r, int i) const
      {
        return A[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R b (int r, int i) const
      {
        return B[r-1][i];
      }

      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R d (int i) const
      {
        return D[i];
      }

      /*! \brief Return name of the scheme
       */
      virtual std::string name () const
      {
        return std::string("explicit Euler");
      }

    private:
      Dune::FieldVector<R,2> D;
      Dune::FieldMatrix<R,1,2> A;
      Dune::FieldMatrix<R,1,2> B;
    };

    /**
     * \brief Parameters to turn the OneStepMethod into an
     * one step theta method.
     *
     * For theta=0 this parameter class can be used with the
     * ExplicitOneStepMethod
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class OneStepThetaParameter : public TimeSteppingParameterInterface<R>
    {
      //! hide default constructor, otherwise theta will be undefined
      OneStepThetaParameter();

    public:
      //! construct OneStepThetaParameter class
      OneStepThetaParameter (R theta_)
      : theta(theta_)
      {
        D[0] = 0.0;  D[1] = 1.0;
        A[0][0] = -1.0; A[0][1] = 1.0;
        B[0][0] = 1.0-theta;  B[0][1] = theta;
      }

      /*! \brief Return true if method is implicit
       */
      virtual bool implicit () const
      {
        if (theta>0.0)
          return true;
        else
          return false;
      }

      /*! \brief Return number of stages s of the method
       */
      virtual unsigned s () const
      {
        return 1;
      }

      /*! \brief Return entries of the A matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R a (int r, int i) const
      {
        return A[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R b (int r, int i) const
      {
        return B[r-1][i];
      }

      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R d (int i) const
      {
        return D[i];
      }

      /*! \brief Return name of the scheme
       */
      virtual std::string name () const
      {
        return std::string("one step theta");
      }

    private:
      R theta;
      Dune::FieldVector<R,2> D;
      Dune::FieldMatrix<R,1,2> A;
      Dune::FieldMatrix<R,1,2> B;
    };

    /**
     * \brief Parameters to turn the ExplicitOneStepMethod into a
     * Heun scheme.
     *
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class HeunParameter : public TimeSteppingParameterInterface<R>
    {
    public:

      HeunParameter ()
      {
        D[0] = 0.0;     D[1] = 1.0;     D[2] = 1.0;

        A[0][0] = -1.0; A[0][1] = 1.0;  A[0][2] = 0.0;
        A[1][0] = -0.5; A[1][1] = -0.5; A[1][2] = 1.0;

        B[0][0] =  1.0; B[0][1] = 0.0;  B[0][2] = 0.0;
        B[1][0] =  0.0; B[1][1] = 0.5;  B[1][2] = 0.0;
      }

      /*! \brief Return true if method is implicit
       */
      virtual bool implicit () const
      {
        return false;
      }

      /*! \brief Return number of stages s of the method
       */
      virtual unsigned s () const
      {
        return 2;
      }

      /*! \brief Return entries of the A matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R a (int r, int i) const
      {
        return A[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R b (int r, int i) const
      {
        return B[r-1][i];
      }

      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R d (int i) const
      {
        return D[i];
      }

      /*! \brief Return name of the scheme
       */
      virtual std::string name () const
      {
        return std::string("Heun");
      }

    private:
      Dune::FieldVector<R,3> D;
      Dune::FieldMatrix<R,2,3> A;
      Dune::FieldMatrix<R,2,3> B;
    };

    /**
     * \brief Parameters to turn the ExplicitOneStepMethod into a
     * third order strong stability preserving (SSP) scheme.
     *
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class Shu3Parameter : public TimeSteppingParameterInterface<R>
    {
    public:

      Shu3Parameter ()
      {
        D[0] = 0.0;     D[1] = 1.0;     D[2] = 0.5; D[3] = 1.0;

        A[0][0] = -1.0; A[0][1] = 1.0; A[0][2] = 0.0; A[0][3] = 0.0;
        A[1][0] = -0.75; A[1][1] = -0.25; A[1][2] = 1.0; A[1][3] = 0.0;
        A[2][0] = -1.0/3.0; A[2][1] = 0.0; A[2][2] = -2.0/3.0; A[2][3] = 1.0;

        B[0][0] =  1.0; B[0][1] = 0.0;  B[0][2] = 0.0;     B[0][3] = 0.0;
        B[1][0] =  0.0; B[1][1] = 0.25;  B[1][2] = 0.0;     B[1][3] = 0.0;
        B[2][0] =  0.0; B[2][1] = 0.0;  B[2][2] = 2.0/3.0; B[2][3] = 0.0;

      }

      /*! \brief Return true if method is implicit
       */
      virtual bool implicit () const
      {
        return false;
      }

      /*! \brief Return number of stages s of the method
       */
      virtual unsigned s () const
      {
        return 3;
      }

      /*! \brief Return entries of the A matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R a (int r, int i) const
      {
        return A[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R b (int r, int i) const
      {
        return B[r-1][i];
      }

      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R d (int i) const
      {
        return D[i];
      }

      /*! \brief Return name of the scheme
       */
      virtual std::string name () const
      {
        return std::string("Shu's third order method");
      }

    private:
      Dune::FieldVector<R,4> D;
      Dune::FieldMatrix<R,3,4> A;
      Dune::FieldMatrix<R,3,4> B;
    };


    /**
     * \brief Parameters to turn the ExplicitOneStepMethod into a
     * classical fourth order Runge-Kutta method
     *
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class RK4Parameter : public TimeSteppingParameterInterface<R>
    {
    public:

      RK4Parameter ()
      {
        D[0] = 0.0;     D[1] = 0.5;     D[2] = 0.5; D[3] = 1.0; D[4] = 1.0;

        A[0][0] = -1.0; A[0][1] = 1.0; A[0][2] = 0.0; A[0][3] = 0.0;  A[0][4] = 0.0;
        A[1][0] = -1.0; A[1][1] = 0.0; A[1][2] = 1.0; A[1][3] = 0.0;  A[1][4] = 0.0;
        A[2][0] = -1.0; A[2][1] = 0.0; A[2][2] = 0.0; A[2][3] = 1.0;  A[2][4] = 0.0;
        A[3][0] = -1.0; A[3][1] = 0.0; A[3][2] = 0.0; A[3][3] = 0.0;  A[3][4] = 1.0;

        B[0][0] =  0.5;     B[0][1] = 0.0;     B[0][2] = 0.0;     B[0][3] = 0.0;      B[0][4] = 0.0;
        B[1][0] =  0.0;     B[1][1] = 0.5;     B[1][2] = 0.0;     B[1][3] = 0.0;      B[1][4] = 0.0;
        B[2][0] =  0.0;     B[2][1] = 0.0;     B[2][2] = 1.0;     B[2][3] = 0.0;      B[2][4] = 0.0;
        B[3][0] =  1.0/6.0; B[3][1] = 1.0/3.0; B[3][2] = 1.0/3.0; B[3][3] = 1.0/6.0;  B[3][4] = 0.0;

      }

      /*! \brief Return true if method is implicit
       */
      virtual bool implicit () const
      {
        return false;
      }

      /*! \brief Return number of stages s of the method
       */
      virtual unsigned s () const
      {
        return 4;
      }

      /*! \brief Return entries of the A matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R a (int r, int i) const
      {
        return A[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R b (int r, int i) const
      {
        return B[r-1][i];
      }

      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R d (int i) const
      {
        return D[i];
      }

      /*! \brief Return name of the scheme
       */
      virtual std::string name () const
      {
        return std::string("RK4");
      }

    private:
      Dune::FieldVector<R,5> D;
      Dune::FieldMatrix<R,4,5> A;
      Dune::FieldMatrix<R,4,5> B;
    };




    /**
     * \brief Parameters to turn the OneStepMethod into an
     * Alexander scheme.
     *
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class Alexander2Parameter : public TimeSteppingParameterInterface<R>
    {
    public:

      Alexander2Parameter ()
      {
        alpha = 1.0 - 0.5*sqrt(2.0);

        D[0] = 0.0;     D[1] = alpha;     D[2] = 1.0;

        A[0][0] = -1.0; A[0][1] = 1.0; A[0][2] = 0.0;
        A[1][0] = -1.0; A[1][1] = 0.0; A[1][2] = 1.0;

        B[0][0] =  0.0; B[0][1] = alpha;  B[0][2] = 0.0;
        B[1][0] =  0.0; B[1][1] = 1.0-alpha;  B[1][2] = alpha;
      }

      /*! \brief Return true if method is implicit
       */
      virtual bool implicit () const
      {
        return true;
      }

      /*! \brief Return number of stages s of the method
       */
      virtual unsigned s () const
      {
        return 2;
      }

      /*! \brief Return entries of the A matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R a (int r, int i) const
      {
        return A[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R b (int r, int i) const
      {
        return B[r-1][i];
      }

      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R d (int i) const
      {
        return D[i];
      }

      /*! \brief Return name of the scheme
       */
      virtual std::string name () const
      {
        return std::string("Alexander (order 2)");
      }

    private:
      R alpha;
      Dune::FieldVector<R,3> D;
      Dune::FieldMatrix<R,2,3> A;
      Dune::FieldMatrix<R,2,3> B;
    };

    /**
     * \brief Parameters to turn the OneStepMethod into a
     * fractional step theta scheme.
     *
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class FractionalStepParameter : public TimeSteppingParameterInterface<R>
    {
    public:

      FractionalStepParameter ()
      {
        R alpha, theta, thetap, beta;
        theta = 1.0 - 0.5*sqrt(2.0);
        thetap = 1.0-2.0*theta;
        alpha = 2.0-sqrt(2.0);
        beta = 1.0-alpha;

        D[0] = 0.0;     D[1] = theta;     D[2] = 1.0-theta; D[3] = 1.0;

        A[0][0] = -1.0; A[0][1] = 1.0; A[0][2] = 0.0; A[0][3] = 0.0;
        A[1][0] = 0.0; A[1][1] = -1.0; A[1][2] = 1.0; A[1][3] = 0.0;
        A[2][0] = 0.0; A[2][1] = 0.0; A[2][2] = -1.0; A[2][3] = 1.0;

        B[0][0] =  beta*theta; B[0][1] = alpha*theta;  B[0][2] = 0.0; B[0][3] = 0.0;
        B[1][0] =  0.0; B[1][1] = alpha*thetap;  B[1][2] = alpha*theta; B[1][3] = 0.0;
        B[2][0] =  0.0; B[2][1] = 0.0; B[2][2] = beta*theta;  B[2][3] = alpha*theta;
      }

      /*! \brief Return true if method is implicit
       */
      virtual bool implicit () const
      {
        return true;
      }

      /*! \brief Return number of stages s of the method
       */
      virtual unsigned s () const
      {
        return 3;
      }

      /*! \brief Return entries of the A matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R a (int r, int i) const
      {
        return A[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R b (int r, int i) const
      {
        return B[r-1][i];
      }

      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R d (int i) const
      {
        return D[i];
      }

      /*! \brief Return name of the scheme
       */
      virtual std::string name () const
      {
        return std::string("Fractional step theta");
      }

    private:
      Dune::FieldVector<R,4> D;
      Dune::FieldMatrix<R,3,4> A;
      Dune::FieldMatrix<R,3,4> B;
    };

    /**

     * \brief Parameters to turn the OneStepMethod into an
     * Alexander3 scheme.
     *
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class Alexander3Parameter : public TimeSteppingParameterInterface<R>
    {
    public:

      Alexander3Parameter ()
      {
        R alpha = 0.4358665215;

        // Newton iteration for alpha
        for (int i=1; i<=10; i++)
          {
            alpha = alpha - (alpha*(alpha*alpha-3.0*(alpha-0.5))-1.0/6.0)/(3.0*alpha*(alpha-2.0)+1.5);
            //         std::cout.precision(16);
            //         std::cout << "alpha "
            //               << std::setw(8) << i << "  "
            //               << std::scientific << alpha << std::endl;
          }

        R tau2 = (1.0+alpha)*0.5;
        R b1 = -(6.0*alpha*alpha -16.0*alpha + 1.0)*0.25;
        R b2 = (6*alpha*alpha - 20.0*alpha + 5.0)*0.25;

        //         std::cout.precision(16);
        //         std::cout << "alpha " << std::scientific << alpha << std::endl;
        //         std::cout << "tau2  " << std::scientific << tau2 << std::endl;
        //         std::cout << "b1    " << std::scientific << b1 << std::endl;
        //         std::cout << "b2    " << std::scientific << b2 << std::endl;

        D[0] = 0.0;     D[1] = alpha;     D[2] = tau2; D[3] = 1.0;

        A[0][0] = -1.0; A[0][1] = 1.0; A[0][2] = 0.0; A[0][3] = 0.0;
        A[1][0] = -1.0; A[1][1] = 0.0; A[1][2] = 1.0; A[1][3] = 0.0;
        A[2][0] = -1.0; A[2][1] = 0.0; A[2][2] = 0.0; A[2][3] = 1.0;

        B[0][0] =  0.0; B[0][1] = alpha;      B[0][2] = 0.0;   B[0][3] = 0.0;
        B[1][0] =  0.0; B[1][1] = tau2-alpha; B[1][2] = alpha; B[1][3] = 0.0;
        B[2][0] =  0.0; B[2][1] = b1;         B[2][2] = b2;    B[2][3] = alpha;
      }

      /*! \brief Return true if method is implicit
       */
      virtual bool implicit () const
      {
        return true;
      }

      /*! \brief Return number of stages s of the method
       */
      virtual unsigned s () const
      {
        return 3;
      }

      /*! \brief Return entries of the A matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R a (int r, int i) const
      {
        return A[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R b (int r, int i) const
      {
        return B[r-1][i];
      }

      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R d (int i) const
      {
        return D[i];
      }

      /*! \brief Return name of the scheme
       */
      virtual std::string name () const
      {
        return std::string("Alexander (claims order 3)");
      }

    private:
      R alpha, theta, thetap, beta;
      Dune::FieldVector<R,4> D;
      Dune::FieldMatrix<R,3,4> A;
      Dune::FieldMatrix<R,3,4> B;
    };


    /**
     * \brief Controller interface for adaptive time stepping.
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class TimeControllerInterface
    {
    public:
      typedef R RealType;

      /*! \brief Return name of the scheme
       */
      virtual RealType suggestTimestep (RealType time, RealType givendt) = 0;

      //! every abstract base class has a virtual destructor
      virtual ~TimeControllerInterface () {}
    };

    //! Default time controller; just returns given dt
    /**
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class SimpleTimeController : public TimeControllerInterface<R>
    {
    public:
      typedef R RealType;

      /*! \brief Return name of the scheme
       */
      virtual RealType suggestTimestep (RealType time, RealType givendt)
      {
        return givendt;
      }
    };


    //! limit time step to maximum dt * CFL number
    /**
     * \tparam R     C++ type of the floating point parameters
     * \tparam IGOS  instationary grid operator space
     */
    template<class R, class IGOS>
    class CFLTimeController : public TimeControllerInterface<R>
    {
    public:
      typedef R RealType;

      CFLTimeController (R cfl_, const IGOS& igos_) : cfl(cfl_), target(1e100), igos(igos_)
      {}

      CFLTimeController (R cfl_, R target_, const IGOS& igos_) : cfl(cfl_), target(target_), igos(igos_)
      {}

      void setTarget (R target_)
      {
        target = target_;
      }

      /*! \brief Return name of the scheme
       */
      virtual RealType suggestTimestep (RealType time, RealType givendt)
      {
        RealType suggested = cfl*igos.suggestTimestep(givendt);
        if (time+2.0*suggested<target)
          return suggested;
        if (time+suggested<target)
          return 0.5*(target-time);
        return target-time;
      }

    private:
      R cfl;
      R target;
      const IGOS& igos;
    };


    // Status information of Newton's method
    struct OneStepMethodPartialResult
    {
      unsigned int timesteps;
      double assembler_time;
      double linear_solver_time;
      int linear_solver_iterations;
      int nonlinear_solver_iterations;

      OneStepMethodPartialResult() :
        timesteps(0),
        assembler_time(0.0),
        linear_solver_time(0.0),
        linear_solver_iterations(0),
        nonlinear_solver_iterations(0)
      {}
    };

    struct OneStepMethodResult
    {
      OneStepMethodPartialResult total;
      OneStepMethodPartialResult successful;
      OneStepMethodResult() : total(), successful()
      {}
    };

    //! Do one step of a time-stepping scheme
    /**
     * \tparam T          type to represent time values
     * \tparam IGOS       assembler for instationary problems
     * \tparam PDESOLVER  solver problem in each step (typically Newton)
     * \tparam TrlV       vector type to represent coefficients of solutions
     * \tparam TstV       vector type to represent residuals
     */
    template<class T, class IGOS, class PDESOLVER, class TrlV, class TstV = TrlV>
    class OneStepMethod
    {
      typedef typename PDESOLVER::Result PDESolverResult;

    public:
      typedef OneStepMethodResult Result;

      //! construct a new one step scheme
      /**
       * \param method_    Parameter object. This chooses the actual method
       *                   used.
       * \param igos_      Assembler object (instationary grid operator space).
       * \param pdesolver_ solver object (typically Newton).
       *
       * The contructed method object stores references to the object it is
       * constructed with, so these objects should be valid for as long as the
       * constructed object is used (or until setMethod() is called, see
       * there).
       */
      OneStepMethod(const TimeSteppingParameterInterface<T>& method_,
                    IGOS& igos_, PDESOLVER& pdesolver_)
        : method(&method_), igos(igos_), pdesolver(pdesolver_), verbosityLevel(1), step(1), res()
      {
        if (igos.trialGridFunctionSpace().gridView().comm().rank()>0)
          verbosityLevel = 0;
      }

      //! change verbosity level; 0 means completely quiet
      void setVerbosityLevel (int level)
      {
        if (igos.trialGridFunctionSpace().gridView().comm().rank()>0)
          verbosityLevel = 0;
        else
          verbosityLevel = level;
      }

      //! change number of current step
      void setStepNumber(int newstep) { step = newstep; }

      //! Access to the (non) linear solver
      const PDESOLVER & getPDESolver() const
      {
        return pdesolver;
      }

      //! Access to the (non) linear solver
      PDESOLVER & getPDESolver()
      {
        return pdesolver;
      }

      const Result& result() const
      {
        return res;
      }

      //! redefine the method to be used; can be done before every step
      /**
       * \param method_ Parameter object.
       *
       * The OneStepMethod object stores a reference to the method_ object.
       * The old method object is no longer referenced after this member
       * function returns.
       */
      void setMethod (const TimeSteppingParameterInterface<T>& method_)
      {
        method = &method_;
      }

      /*! \brief do one step;
       * \param[in]  time start of time step
       * \param[in]  dt suggested time step size
       * \param[in]  xold value at begin of time step
       * \param[in,out] xnew value at end of time step; contains initial guess for first substep on entry
       * \return selected time step size
       */
      T apply (T time, T dt, TrlV& xold, TrlV& xnew)
      {
        // save formatting attributes
        ios_base_all_saver format_attribute_saver(std::cout);

        // do statistics
        OneStepMethodPartialResult step_result;

        std::vector<TrlV*> x(1); // vector of pointers to all steps
        x[0] = &xold;            // initially we have only one

        if (verbosityLevel>=1){
          std::ios_base::fmtflags oldflags = std::cout.flags();
          std::cout << "TIME STEP [" << method->name() << "] "
                    << std::setw(6) << step
                    << " time (from): "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << time
                    << " dt: "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << dt
                    << " time (to): "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << time+dt
                    << std::endl;
          std::cout.flags(oldflags);
        }

        // prepare assembler
        igos.preStep(*method,time,dt);

        // loop over all stages
        for (unsigned r=1; r<=method->s(); ++r)
          {
            if (verbosityLevel>=2){
              std::ios_base::fmtflags oldflags = std::cout.flags();
              std::cout << "STAGE "
                        << r
                        << " time (to): "
                        << std::setw(12) << std::setprecision(4) << std::scientific
                        << time+method->d(r)*dt
                        << "." << std::endl;
              std::cout.flags(oldflags);
            }

            // prepare stage
            igos.preStage(r,x);

            // get vector for current stage
            if (r==method->s())
              {
                // last stage
                x.push_back(&xnew);
                if (r>1) xnew = *(x[r-1]); // if r=1 then xnew has already initial guess
              }
            else
              {
                // intermediate step
                x.push_back(new TrlV(igos.trialGridFunctionSpace()));
                if (r>1)
                  *(x[r]) = *(x[r-1]); // use result of last stage as initial guess
                else
                  *(x[r]) = xnew;
              }

            // solve stage
            try {
              pdesolver.apply(*x[r]);
            }
            catch (...)
              {
                // time step failed -> accumulate to total only
                PDESolverResult pderes = pdesolver.result();
                step_result.assembler_time += pderes.assembler_time;
                step_result.linear_solver_time += pderes.linear_solver_time;
                step_result.linear_solver_iterations += pderes.linear_solver_iterations;
                step_result.nonlinear_solver_iterations += pderes.iterations;
                res.total.assembler_time += step_result.assembler_time;
                res.total.linear_solver_time += step_result.linear_solver_time;
                res.total.linear_solver_iterations += step_result.linear_solver_iterations;
                res.total.nonlinear_solver_iterations += step_result.nonlinear_solver_iterations;
                res.total.timesteps += 1;
                throw;
              }
            PDESolverResult pderes = pdesolver.result();
            step_result.assembler_time += pderes.assembler_time;
            step_result.linear_solver_time += pderes.linear_solver_time;
            step_result.linear_solver_iterations += pderes.linear_solver_iterations;
            step_result.nonlinear_solver_iterations += pderes.iterations;

            // stage cleanup
            igos.postStage();
          }

        // delete intermediate steps
        for (unsigned i=1; i<method->s(); ++i) delete x[i];

        // step cleanup
        igos.postStep();

        // update statistics
        res.total.assembler_time += step_result.assembler_time;
        res.total.linear_solver_time += step_result.linear_solver_time;
        res.total.linear_solver_iterations += step_result.linear_solver_iterations;
        res.total.nonlinear_solver_iterations += step_result.nonlinear_solver_iterations;
        res.total.timesteps += 1;
        res.successful.assembler_time += step_result.assembler_time;
        res.successful.linear_solver_time += step_result.linear_solver_time;
        res.successful.linear_solver_iterations += step_result.linear_solver_iterations;
        res.successful.nonlinear_solver_iterations += step_result.nonlinear_solver_iterations;
        res.successful.timesteps += 1;
        if (verbosityLevel>=1){
          std::ios_base::fmtflags oldflags = std::cout.flags();
          std::cout << "::: timesteps      " << std::setw(6) << res.successful.timesteps
                    << " (" << res.total.timesteps << ")" << std::endl;
          std::cout << "::: nl iterations  " << std::setw(6) << res.successful.nonlinear_solver_iterations
                    << " (" << res.total.nonlinear_solver_iterations << ")" << std::endl;
          std::cout << "::: lin iterations " << std::setw(6) << res.successful.linear_solver_iterations
                    << " (" << res.total.linear_solver_iterations << ")" << std::endl;
          std::cout << "::: assemble time  " << std::setw(12) << std::setprecision(4) << std::scientific
                    << res.successful.assembler_time << " (" << res.total.assembler_time << ")" << std::endl;
          std::cout << "::: lin solve time " << std::setw(12) << std::setprecision(4) << std::scientific
                    << res.successful.linear_solver_time << " (" << res.total.linear_solver_time << ")" << std::endl;
          std::cout.flags(oldflags);
        }

        step++;
        return dt;
      }

      /*! \brief do one step;
       * This is a version which interpolates constraints at the start of each stage
       *
       * \param[in]  time start of time step
       * \param[in]  dt suggested time step size
       * \param[in]  xold value at begin of time step
       * \param[in]  f function to interpolate boundary conditions from
       * \param[in,out] xnew value at end of time step; contains initial guess for first substep on entry
       * \return selected time step size
       */
      template<typename F>
      T apply (T time, T dt, TrlV& xold, F& f, TrlV& xnew)
      {
        // do statistics
        OneStepMethodPartialResult step_result;

        // save formatting attributes
        ios_base_all_saver format_attribute_saver(std::cout);

        std::vector<TrlV*> x(1); // vector of pointers to all steps
        x[0] = &xold;            // initially we have only one

        if (verbosityLevel>=1){
          std::ios_base::fmtflags oldflags = std::cout.flags();
          std::cout << "TIME STEP [" << method->name() << "] "
                    << std::setw(6) << step
                    << " time (from): "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << time
                    << " dt: "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << dt
                    << " time (to): "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << time+dt
                    << std::endl;
          std::cout.flags(oldflags);
        }

        // prepare assembler
        igos.preStep(*method,time,dt);

        // loop over all stages
        for (unsigned r=1; r<=method->s(); ++r)
          {
            if (verbosityLevel>=2){
              std::ios_base::fmtflags oldflags = std::cout.flags();
              std::cout << "STAGE "
                        << r
                        << " time (to): "
                        << std::setw(12) << std::setprecision(4) << std::scientific
                        << time+method->d(r)*dt
                        << "." << std::endl;
              std::cout.flags(oldflags);
            }

            // prepare stage
            igos.preStage(r,x);

            // get vector for current stage
            if (r==method->s())
              {
                // last stage
                x.push_back(&xnew);
              }
            else
              {
                // intermediate step
                x.push_back(new TrlV(igos.trialGridFunctionSpace()));
              }

            // set boundary conditions and initial value
            igos.interpolate(r,*x[r-1],f,*x[r]);

            // solve stage
            try {
              pdesolver.apply(*x[r]);
            }
            catch (...)
              {
                // time step failed -> accumulate to total only
                PDESolverResult pderes = pdesolver.result();
                step_result.assembler_time += pderes.assembler_time;
                step_result.linear_solver_time += pderes.linear_solver_time;
                step_result.linear_solver_iterations += pderes.linear_solver_iterations;
                step_result.nonlinear_solver_iterations += pderes.iterations;
                res.total.assembler_time += step_result.assembler_time;
                res.total.linear_solver_time += step_result.linear_solver_time;
                res.total.linear_solver_iterations += step_result.linear_solver_iterations;
                res.total.nonlinear_solver_iterations += step_result.nonlinear_solver_iterations;
                res.total.timesteps += 1;
                throw;
              }
            PDESolverResult pderes = pdesolver.result();
            step_result.assembler_time += pderes.assembler_time;
            step_result.linear_solver_time += pderes.linear_solver_time;
            step_result.linear_solver_iterations += pderes.linear_solver_iterations;
            step_result.nonlinear_solver_iterations += pderes.iterations;

            // stage cleanup
            igos.postStage();
          }

        // delete intermediate steps
        for (unsigned i=1; i<method->s(); ++i) delete x[i];

        // step cleanup
        igos.postStep();

        // update statistics
        res.total.assembler_time += step_result.assembler_time;
        res.total.linear_solver_time += step_result.linear_solver_time;
        res.total.linear_solver_iterations += step_result.linear_solver_iterations;
        res.total.nonlinear_solver_iterations += step_result.nonlinear_solver_iterations;
        res.total.timesteps += 1;
        res.successful.assembler_time += step_result.assembler_time;
        res.successful.linear_solver_time += step_result.linear_solver_time;
        res.successful.linear_solver_iterations += step_result.linear_solver_iterations;
        res.successful.nonlinear_solver_iterations += step_result.nonlinear_solver_iterations;
        res.successful.timesteps += 1;
        if (verbosityLevel>=1){
          std::ios_base::fmtflags oldflags = std::cout.flags();
          std::cout << "::: timesteps      " << std::setw(6) << res.successful.timesteps
                    << " (" << res.total.timesteps << ")" << std::endl;
          std::cout << "::: nl iterations  " << std::setw(6) << res.successful.nonlinear_solver_iterations
                    << " (" << res.total.nonlinear_solver_iterations << ")" << std::endl;
          std::cout << "::: lin iterations " << std::setw(6) << res.successful.linear_solver_iterations
                    << " (" << res.total.linear_solver_iterations << ")" << std::endl;
          std::cout << "::: assemble time  " << std::setw(12) << std::setprecision(4) << std::scientific
                    << res.successful.assembler_time << " (" << res.total.assembler_time << ")" << std::endl;
          std::cout << "::: lin solve time " << std::setw(12) << std::setprecision(4) << std::scientific
                    << res.successful.linear_solver_time << " (" << res.total.linear_solver_time << ")" << std::endl;
          std::cout.flags(oldflags);
        }

        step++;
        return dt;
      }

    private:
      const TimeSteppingParameterInterface<T> *method;
      IGOS& igos;
      PDESOLVER& pdesolver;
      int verbosityLevel;
      int step;
      Result res;
    };

    //! Do one step of an explicit time-stepping scheme
    /**
     * \tparam T          type to represent time values
     * \tparam IGOS       assembler for instationary problems
     * \tparam LS         backend to solve diagonal linear system
     * \tparam TrlV       vector type to represent coefficients of solutions
     * \tparam TstV       vector type to represent residuals
     * \tparam TC         time controller class
     */
    template<class T, class IGOS, class LS, class TrlV, class TstV = TrlV, class TC = SimpleTimeController<T> >
    class ExplicitOneStepMethod
    {
      typedef typename TrlV::ElementType Real;
      typedef typename IGOS::template MatrixContainer<Real>::Type M;

    public:
      //! construct a new one step scheme
      /**
       * \param method_    Parameter object.
       * \param igos_      Assembler object (instationary grid operator space).
       * \param pdesolver_ solver object (typically Newton).
       *
       * The contructed method object stores references to the object it is
       * constructed with, so these objects should be valid for as long as the
       * constructed object is used (or until setMethod() is called, see
       * there).
       * Use SimpleTimeController that does not control the time step.
       */
      ExplicitOneStepMethod(const TimeSteppingParameterInterface<T>& method_, IGOS& igos_, LS& ls_)
        : method(&method_), igos(igos_), ls(ls_), verbosityLevel(1), step(1), D(igos),
          tc(new SimpleTimeController<T>()), allocated(true)
      {
        if (method->implicit())
          DUNE_THROW(Exception,"explicit one step method called with implicit scheme");
        if (igos.trialGridFunctionSpace().gridView().comm().rank()>0)
          verbosityLevel = 0;
      }

      //! construct a new one step scheme
      /**
       * \param method_    Parameter object.
       * \param igos_      Assembler object (instationary grid operator space).
       * \param pdesolver_ solver object (typically Newton).
       * \param tc_        a time controller object
       *
       * The contructed method object stores references to the object it is
       * constructed with, so these objects should be valid for as long as the
       * constructed object is used (or until setMethod() is called, see
       * there).
       */
      ExplicitOneStepMethod(const TimeSteppingParameterInterface<T>& method_, IGOS& igos_, LS& ls_, TC& tc_)
        : method(&method_), igos(igos_), ls(ls_), verbosityLevel(1), step(1), D(igos),
          tc(&tc_), allocated(false)
      {
        if (method->implicit())
          DUNE_THROW(Exception,"explicit one step method called with implicit scheme");
      }

      ~ExplicitOneStepMethod ()
      {
        if (allocated) delete tc;
      }

      //! change verbosity level; 0 means completely quiet
      void setVerbosityLevel (int level)
      {
        if (igos.trialGridFunctionSpace().gridView().comm().rank()>0)
          verbosityLevel = 0;
        else
          verbosityLevel = level;
      }

      //! change number of current step
      void setStepNumber(int newstep) { step = newstep; }

      //! redefine the method to be used; can be done before every step
      /**
       * \param method_ Parameter object.
       *
       * The OneStepMethod object stores a reference to the method_ object.
       * The old method object is no longer referenced after this member
       * function returns.
       */
      void setMethod (const TimeSteppingParameterInterface<T>& method_)
      {
        method = &method_;
        if (method->implicit())
          DUNE_THROW(Exception,"explicit one step method called with implicit scheme");
      }

      /*! \brief do one step;
       * \param[in]  time start of time step
       * \param[in]  dt suggested time step size
       * \param[in]  xold value at begin of time step
       * \param[in]  limiter
       * \param[in,out] xnew value at end of time step; contains initial guess for first substep on entry
       * \return time step size
       */
      T apply (T time, T dt, TrlV& xold, TrlV& xnew)
      {
        DefaultLimiter limiter;
        return apply(time,dt,xold,xnew,limiter);
      }

      template<typename Limiter>
      T apply (T time, T dt, TrlV& xold, TrlV& xnew, Limiter& limiter)
      {
        // save formatting attributes
        ios_base_all_saver format_attribute_saver(std::cout);
        LocalTag mytag;
        mytag << "ExplicitOneStepMethod::apply(): ";

        std::vector<TrlV*> x(1); // vector of pointers to all steps
        x[0] = &xold;         // initially we have only one
        if(verbosityLevel>=4)
          std::cout << mytag << "Creating residual vectors alpha and beta..."
                    << std::endl;
        TstV alpha(igos.testGridFunctionSpace()), beta(igos.testGridFunctionSpace()); // split residual vectors
        if(verbosityLevel>=4)
          std::cout << mytag
                    << "Creating residual vectors alpha and beta... done."
                    << std::endl;

        if (verbosityLevel>=1){
          std::ios_base::fmtflags oldflags = std::cout.flags();
          std::cout << "TIME STEP [" << method->name() << "] "
                    << std::setw(6) << step
                    << " time (from): "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << time
                    << " dt: "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << dt
                    << " time (to): "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << time+dt
                    << std::endl;
          std::cout.flags(oldflags);
        }

        // prepare assembler
        if(verbosityLevel>=4)
          std::cout << mytag << "Preparing assembler..." << std::endl;
        igos.preStep(*method,time,dt);
        if(verbosityLevel>=4)
          std::cout << mytag << "Preparing assembler... done." << std::endl;

        // loop over all stages
        for(unsigned r=1; r<=method->s(); ++r)
          {
            LocalTag stagetag(mytag);
            stagetag << "stage " << r << ": ";
            if (verbosityLevel>=4)
              std::cout << stagetag << "Start." << std::endl;

            if (verbosityLevel>=2){
              std::ios_base::fmtflags oldflags = std::cout.flags();
              std::cout << "STAGE "
                        << r
                        << " time (to): "
                        << std::setw(12) << std::setprecision(4) << std::scientific
                        << time+method->d(r)*dt
                        << "." << std::endl;
              std::cout.flags(oldflags);
            }

            // get vector for current stage
            if (r==method->s())
              {
                // last stage
                x.push_back(&xnew);
                if (r>1) xnew = *(x[r-1]); // if r=1 then xnew has already initial guess
              }
            else
              {
                // intermediate step
                x.push_back(new TrlV(igos.trialGridFunctionSpace()));
                if (r>1)
                  *(x[r]) = *(x[r-1]); // use result of last stage as initial guess
                else
                  *(x[r]) = xnew;
              }

            // compute residuals and jacobian
            if (verbosityLevel>=4) std::cout << "assembling D, alpha, beta ..." << std::endl;
            D = 0.0;
            alpha = 0.0;
            beta = 0.0;

            //apply slope limiter to old solution (e.g for finite volume reconstruction scheme)
            limiter.prestage(*x[r-1]);

            if(verbosityLevel>=4)
              std::cout << stagetag << "Assembling residual..." << std::endl;
            igos.explicit_jacobian_residual(r,x,D,alpha,beta);
            if(verbosityLevel>=4)
              std::cout << stagetag << "Assembling residual... done."
                        << std::endl;

            // let time controller compute the optimal dt in first stage
            if (r==1)
              {
                T newdt = tc->suggestTimestep(time,dt);
                newdt = std::min(newdt, dt);

                if (verbosityLevel>=4){
                  std::ios_base::fmtflags oldflags = std::cout.flags();
                  std::cout << stagetag
                            << "current dt: "
                            << std::setw(12) << std::setprecision(4) << std::scientific
                            << dt
                            << " suggested dt: "
                            << std::setw(12) << std::setprecision(4) << std::scientific
                            << newdt
                            << std::endl;
                  std::cout.flags(oldflags);
                }

                if (verbosityLevel>=2 && newdt!=dt)
                  {
                    std::ios_base::fmtflags oldflags = std::cout.flags();
                    std::cout << "changed dt to "
                              << std::setw(12) << std::setprecision(4) << std::scientific
                              << newdt
                              << std::endl;
                    std::cout.flags(oldflags);
                  }
                dt = newdt;
              }

            // combine residual with selected dt
            if (verbosityLevel>=4)
              std::cout << stagetag
                        << "Combining residuals with selected dt..."
                        << std::endl;
            alpha.axpy(dt,beta);
            if (verbosityLevel>=4)
              std::cout << stagetag
                        << "Combining residuals with selected dt... done."
                        << std::endl;

            // solve diagonal system
            if (verbosityLevel>=4)
              std::cout << stagetag << "Solving diagonal system..."
                        << std::endl;
            ls.apply(D,*x[r],alpha,0.99); // dummy reduction
            if (verbosityLevel>=4)
              std::cout << stagetag << "Solving diagonal system... done."
                        << std::endl;

            // apply slope limiter to new solution (e.g DG scheme)
            limiter.poststage(*x[r]);

            // stage cleanup
            if (verbosityLevel>=4)
              std::cout << stagetag << "Cleanup..." << std::endl;
            igos.postStage();
            if (verbosityLevel>=4)
              std::cout << stagetag << "Cleanup... done" << std::endl;

            if (verbosityLevel>=4)
              std::cout << stagetag << "Finished." << std::endl;
          }

        // delete intermediate steps
        for(unsigned i=1; i<method->s(); ++i) delete x[i];

        // step cleanup
        if (verbosityLevel>=4)
          std::cout << mytag << "Cleanup..." << std::endl;
        igos.postStep();
        if (verbosityLevel>=4)
          std::cout << mytag << "Cleanup... done." << std::endl;

        step++;
        return dt;
      }

    private:

      //! dummy default limiter
      class DefaultLimiter
      {
      public:
        template<typename V>
        void prestage(V& v)
        {}

        template<typename V>
        void poststage(V& v)
        {}
      };

      const TimeSteppingParameterInterface<T> *method;
      IGOS& igos;
      LS& ls;
      int verbosityLevel;
      int step;
      M D;
      TimeControllerInterface<T> *tc;
      bool allocated;
    };

    class FilenameHelper
    {
    public:
      FilenameHelper(const char *basename_, int i_=0)
        : i(i_)
      {
        sprintf(basename,"%s",basename_);
      }

      FilenameHelper(const std::string & basename_, int i_=0)
        : i(i_)
      {
        sprintf(basename,"%s",basename_.c_str());
      }

      const char *getName (int i_)
      {
        sprintf(fname,"%s-%05d",basename,i_);
        return fname;
      }

      const char *getName ()
      {
        sprintf(fname,"%s-%05d",basename,i);
        return fname;
      }

      void increment ()
      {
        i++;
      }

    private:
      char fname[255];
      char basename[255];
      int i;
    };

    /** @} */
  } // namespace PDELab
} // namespace Dune

#endif
