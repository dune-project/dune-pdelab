// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_INSTATIONARY_ONESTEPPARAMETER_HH
#define DUNE_PDELAB_INSTATIONARY_ONESTEPPARAMETER_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

namespace Dune {
  namespace PDELab {

    /**
     *  @addtogroup OneStepMethod
     *  @{
     */
    //! Base parameter class for time stepping scheme parameters
    /**
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class TimeSteppingParameterInterface
    {
    public:
      typedef R RealType;

      /*! \brief Return true if method is implicit
       */
      virtual bool implicit () const = 0;

      /*! \brief Return number of stages of the method
       */
      virtual unsigned s () const = 0;

      /*! \brief Return entries of the A matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R a (int r, int i) const = 0;

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R b (int r, int i) const = 0;

      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,r
      */
      virtual R d (int r) const = 0;

      /*! \brief Return name of the scheme
       */
      virtual std::string name () const = 0;

      //! every abstract base class has a virtual destructor
      virtual ~TimeSteppingParameterInterface () {}
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
        return (theta > 0.0);
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
     * \brief Parameters to turn the ExplicitOneStepMethod into an
     * explicit Euler method.
     *
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class ExplicitEulerParameter : public OneStepThetaParameter<R>
    {
    public:

      ExplicitEulerParameter ()
        : OneStepThetaParameter<R>(0.0)
      {}

      /*! \brief Return name of the scheme
       */
      virtual std::string name () const
      {
        return std::string("explicit Euler");
      }

    };


   /*!
    * \brief Parameters to turn the OneStepMethod into an
    * implicit Euler method.
    *
    * \tparam R C++ type of the floating point parameters
    */
   template<class R>
   class ImplicitEulerParameter : public OneStepThetaParameter<R>
   {
   public:

     ImplicitEulerParameter ()
       : OneStepThetaParameter<R>(1.0)
     {}

     /*! \brief Return name of the scheme
      */
     virtual std::string name () const
     {
       return std::string("implicit Euler");
     }

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

  } // end namespace PDELab
} // end namespace Dune
#endif // DUNE_PDELAB_INSTATIONARY_ONESTEPPARAMETER_HH
