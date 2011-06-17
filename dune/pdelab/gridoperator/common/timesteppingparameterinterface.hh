// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_TIMESTEPPINGPARAMETERINTERFACE_HH
#define DUNE_PDELAB_TIMESTEPPINGPARAMETERINTERFACE_HH

#include <string>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dune {
  namespace PDELab {

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
        return std::string("implicit Euler");
      }

    private:
      Dune::FieldVector<R,2> D;
      Dune::FieldMatrix<R,1,2> A;
      Dune::FieldMatrix<R,1,2> B;
    };

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
