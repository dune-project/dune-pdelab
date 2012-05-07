#ifndef DUNE_PDELAB_BACKEND_SOLVER_HH
#define DUNE_PDELAB_BACKEND_SOLVER_HH

namespace Dune {
  namespace PDELab {

    //! \addtogroup Backend
    //! \ingroup PDELab
    //! \{

    struct SequentialNorm
    {/*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      template<class V>
      typename Dune::template FieldTraits<typename V::ElementType >::real_type norm(const V& v) const
      {
        return v.base().two_norm();
      }
    };

    class LinearResultStorage
    {
    public:
      /*! \brief Return access to result data */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res;
      }

    protected:
      Dune::PDELab::LinearSolverResult<double> res;
    };
    
    //! \} group Backend
    
  } // end namespace PDELab
} // end namespace Dune


#endif // DUNE_PDELAB_BACKEND_SOLVER_HH
