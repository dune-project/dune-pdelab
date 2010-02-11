// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ONESTEP_HH
#define DUNE_PDELAB_ONESTEP_HH

#include <iostream>
#include <iomanip>
#include <vector>

#include <stdio.h>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>

#include"../gridoperatorspace/instationarygridoperatorspace.hh"

namespace Dune {
  namespace PDELab {

    //===============================================================
    // Interface class for parameters of one step methods
    //===============================================================

    // some implementations of this interface

    //! Parameters specifying implicit euler
    /**
     * \tparam R C++ type of the floating point parameters
     */
    template<class R> 
    class ImplicitEulerParameter : public TimeSteppingParameterInterface<R>
    {
    public:
      
      ImplicitEulerParameter ()
      {
	D[0] = 0.0;  D[1] = 1.0;
	A[0][0] = -1.0; A[0][1] = 1.0;
	B[0][0] = 0.0;  B[0][1] = 1.0;
      }

      /*! \brief Return true if method is implicit
      */
      virtual bool implicit () const
      {
	return true;
      }
      
      /*! \brief Return number of stages s of the method
      */
      virtual int s () const
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
        return std::string("implicit Euler");
      }

    private:
      Dune::FieldVector<R,2> D;
      Dune::FieldMatrix<R,1,2> A;
      Dune::FieldMatrix<R,1,2> B;
    };

    //! Parameters specifying explicit euler
    /**
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
      virtual int s () const
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

    //! Parameters specifying the one step theta scheme
    /**
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
      virtual int s () const
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

    //! Parameters specifying the Heun scheme
    /**
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
      virtual int s () const
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

    //! Parameters specifying the Alexander scheme
    /**
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
      virtual int s () const
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

    //! Parameters specifying the fractional step theta scheme
    /**
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
      virtual int s () const
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

    //! Parameters specifying the Alexander3 scheme
    /**
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
// 	    std::cout.precision(16);
// 	    std::cout << "alpha " 
// 		      << std::setw(8) << i << "  " 
// 		      << std::scientific << alpha << std::endl;
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
      virtual int s () const
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
       */
      OneStepMethod(const TimeSteppingParameterInterface<T>& method_, IGOS& igos_, PDESOLVER& pdesolver_)
	: method(&method_), igos(igos_), pdesolver(pdesolver_), verbosityLevel(1), step(1)
      {
      }

      //! change verbosity level; 0 means completely quiet
      void setVerbosityLevel (int level)
      {
	verbosityLevel = level;
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

      //! do one step;
      /*
       * \param[in]  xold value at begin of time step
       * \param[out] xnew value at end of time step; contains initial guess for first substep on entry
       */
      void apply (T time, T dt, TrlV& xold, TrlV& xnew)
      {
	std::vector<TrlV*> x(1); // vector of pointers to all steps
	x[0] = &xold;            // initially we have only one
	TstV residual0(igos.testGridFunctionSpace()); // stores constant part of residual

	if (verbosityLevel>=1)
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

	// prepare assembler
	igos.preStep(*method,time,dt);

	// loop over all stages
	for (int r=1; r<=method->s(); ++r)
	  {
	    if (verbosityLevel>=2)
	      std::cout << "STAGE " 
                        << r 
                        << " time (to): "
                        << std::setw(12) << std::setprecision(4) << std::scientific
                        << time+method->d(r)*dt
                        << "." << std::endl;
	      
	    // prepare stage
	    residual0 = 0.0;
	    igos.preStage(r,x,residual0);

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
	    pdesolver.apply(*x[r]);
	  }

	// delete intermediate steps
	for (int i=1; i<method->s(); ++i) delete x[i];

	step++;
      }

    private:
      const TimeSteppingParameterInterface<T> *method;
      IGOS& igos;
      PDESOLVER pdesolver;
      int verbosityLevel;
      int step;
    };

    class FilenameHelper 
    {
    public:
      FilenameHelper(const char *basename_, int i_=0)
	: i(i_)
      {
	sprintf(basename,"%s",basename_);
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

  } // namespace PDELab
} // namespace Dune

#endif
