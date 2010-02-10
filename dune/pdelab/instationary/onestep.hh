#ifndef DUNE_PDELAB_ONESTEP_HH
#define DUNE_PDELAB_ONESTEP_HH

#include <iostream>
#include <iomanip>
#include <vector>

#include <stdio.h>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>

namespace Dune {
  namespace PDELab {

    //===============================================================
    // Interface class for parameters of one step methods
    //===============================================================

    template<class R> 
    class OneStepParameterInterface
    {
    public:
      typedef R RealType;
      
      /*! \brief Return true if method is implicit
      */
      virtual bool implicit () const = 0;
      
      /*! \brief Return number of stages of the method
      */
      virtual int s () const = 0;
      
      /*! \brief Return entries of the A matrix
	Note that r \in 1,...,s and i \in 0,...,r
      */
      virtual R a (int r, int i) const = 0;
      
      /*! \brief Return entries of the B matrix
	Note that r \in 1,...,s and i \in 0,...,r
      */
      virtual R b (int r, int i) const = 0;
      
      /*! \brief Return entries of the d Vector
	Note that i \in 0,...,r
      */
      virtual R d (int r) const = 0;
      
      //! every abstract base class has a virtual destructor
      virtual ~OneStepParameterInterface () {}
    };

    // some implementations of this interface

    template<class R> 
    class ImplicitEulerParameter : public OneStepParameterInterface<R>
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
	Note that r \in 1,...,s and i \in 0,...,r
      */
      virtual R a (int r, int i) const
      {
	return A[r-1][i];
      }
      
      /*! \brief Return entries of the B matrix
	Note that r \in 1,...,s and i \in 0,...,r
      */
      virtual R b (int r, int i) const
      {
	return B[r-1][i];
      }
      
      /*! \brief Return entries of the d Vector
	i runs from 0,...,s
      */
      virtual R d (int i) const
      {
	return D[i];
      }
      
    private:
      Dune::FieldVector<R,2> D;
      Dune::FieldMatrix<R,1,2> A;
      Dune::FieldMatrix<R,1,2> B;
    };

    template<class R> 
    class ExplicitEulerParameter : public OneStepParameterInterface<R>
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
	Note that r \in 1,...,s and i \in 0,...,r
      */
      virtual R a (int r, int i) const
      {
	return A[r-1][i];
      }
      
      /*! \brief Return entries of the B matrix
	Note that r \in 1,...,s and i \in 0,...,r
      */
      virtual R b (int r, int i) const
      {
	return B[r-1][i];
      }
      
      /*! \brief Return entries of the d Vector
	i runs from 0,...,s
      */
      virtual R d (int i) const
      {
	return D[i];
      }
      
    private:
      Dune::FieldVector<R,2> D;
      Dune::FieldMatrix<R,1,2> A;
      Dune::FieldMatrix<R,1,2> B;
    };

    template<class R> 
    class OneStepThetaParameter : public OneStepParameterInterface<R>
    {
    public:
      
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
	Note that r \in 1,...,s and i \in 0,...,r
      */
      virtual R a (int r, int i) const
      {
	return A[r-1][i];
      }
      
      /*! \brief Return entries of the B matrix
	Note that r \in 1,...,s and i \in 0,...,r
      */
      virtual R b (int r, int i) const
      {
	return B[r-1][i];
      }
      
      /*! \brief Return entries of the d Vector
	i runs from 0,...,s
      */
      virtual R d (int i) const
      {
	return D[i];
      }
      
    private:
      R theta;
      Dune::FieldVector<R,2> D;
      Dune::FieldMatrix<R,1,2> A;
      Dune::FieldMatrix<R,1,2> B;
    };

    template<class R> 
    class HeunParameter : public OneStepParameterInterface<R>
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
	Note that r \in 1,...,s and i \in 0,...,r
      */
      virtual R a (int r, int i) const
      {
	return A[r-1][i];
      }
      
      /*! \brief Return entries of the B matrix
	Note that r \in 1,...,s and i \in 0,...,r
      */
      virtual R b (int r, int i) const
      {
	return B[r-1][i];
      }
      
      /*! \brief Return entries of the d Vector
	i runs from 0,...,s
      */
      virtual R d (int i) const
      {
	return D[i];
      }
      
    private:
      Dune::FieldVector<R,3> D;
      Dune::FieldMatrix<R,2,3> A;
      Dune::FieldMatrix<R,2,3> B;
    };

    template<class R> 
    class AlexanderParameter : public OneStepParameterInterface<R>
    {
    public:
      
      AlexanderParameter ()
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
	Note that r \in 1,...,s and i \in 0,...,r
      */
      virtual R a (int r, int i) const
      {
	return A[r-1][i];
      }
      
      /*! \brief Return entries of the B matrix
	Note that r \in 1,...,s and i \in 0,...,r
      */
      virtual R b (int r, int i) const
      {
	return B[r-1][i];
      }
      
      /*! \brief Return entries of the d Vector
	i runs from 0,...,s
      */
      virtual R d (int i) const
      {
	return D[i];
      }
      
    private:
      R alpha;
      Dune::FieldVector<R,3> D;
      Dune::FieldMatrix<R,2,3> A;
      Dune::FieldMatrix<R,2,3> B;
    };

    template<class R> 
    class FractionalStepParameter : public OneStepParameterInterface<R>
    {
    public:
      
      FractionalStepParameter ()
      {
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
	Note that r \in 1,...,s and i \in 0,...,r
      */
      virtual R a (int r, int i) const
      {
	return A[r-1][i];
      }
      
      /*! \brief Return entries of the B matrix
	Note that r \in 1,...,s and i \in 0,...,r
      */
      virtual R b (int r, int i) const
      {
	return B[r-1][i];
      }
      
      /*! \brief Return entries of the d Vector
	i runs from 0,...,s
      */
      virtual R d (int i) const
      {
	return D[i];
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
      OneStepMethod(const OneStepParameterInterface<T>& method_, IGOS& igos_, PDESOLVER& pdesolver_)
	: method(&method_), igos(igos_), pdesolver(pdesolver_), verbosityLevel(1), step(1)
      {
      }

      //! change verbosity level; 0 means completely quiet
      void setVerbosityLevel (int level)
      {
	verbosityLevel = level;
      }

      //! redefine the method to be used; can be done before every step
      void setMethod (const OneStepParameterInterface<T>& method_)
      {
	method = &method;
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
	  std::cout << "TIME STEP " << std::setw(6) << step
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
	      std::cout << "STAGE " << r << std::endl;
	      
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
      const OneStepParameterInterface<T> *method;
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
