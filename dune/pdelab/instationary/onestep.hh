#ifndef DUNE_PDELAB_ONESTEP_HH
#define DUNE_PDELAB_ONESTEP_HH

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

  } // namespace PDELab
} // namespace Dune

#endif
