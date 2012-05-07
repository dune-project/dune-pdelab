#ifndef DUNE_PDELAB_EIGENVECTORBACKEND_HH
#define DUNE_PDELAB_EIGENVECTORBACKEND_HH

#include<vector>

#if HAVE_EIGEN

#include <Eigen/Dense>

namespace Dune {
  namespace PDELab {

    //! ISTL backend for FunctionSpace
    class EigenVectorBackend
    {
    public:
      static const int BlockSize = 1;

      //! container construction
      template<typename T, typename E>
        class VectorContainer : public Eigen::VectorXd
      {
      public:
        typedef E ElementType;
        typedef Eigen::VectorXd BaseT;
        typedef EigenVectorBackend Backend;

        VectorContainer (const T& t_) : BaseT(t_.globalSize())
        {
        }
        VectorContainer (const T& t_, const E& e) : BaseT(t_.globalSize())
        {
          (*this).setConstant(t_.globalSize(), e);
        }
        VectorContainer& operator= (const E& e)
        {
          (*this).setConstant(this->rows(), e);
          return *this;
        }

        // for debugging and AMG access
        BaseT& base ()
        {
          return *this;
        }

        const BaseT& base () const
        {
          return *this;
        }

        size_t flatsize() const
        {
          return this->rows();
        }

        template<typename X>
          void std_copy_to (std::vector<X>& x) const
          {
            size_t n = flatsize();
            x.resize(n);
            for (Eigen::VectorXd::Index row=0; row<n; row++)
              x[row] = (*this)(row);
          }

        template<typename X>
          void std_copy_from (const std::vector<X>& x)
          {
            size_t n = flatsize();
            x.resize(n);
            for (Eigen::VectorXd::Index row=0; row<n; row++)
              (*this)(row) = x[row];
          }
      };

      // extract type of container element 
      template<class C>
        struct Value
        {
          typedef typename C::ElementType Type;
        };

      //! The size type
      typedef typename Eigen::VectorXd::Index size_type;

      // get const_reference to container element
      // note: this method does not depend on T!
      template<typename C>
        static const typename C::ElementType& access (const C& c, size_type row)
        {
          return c(row);
        }

      // get non const_reference to container element 
      // note: this method does not depend on T!
      template<typename C>
        static typename C::ElementType& access(C& c, size_type row)
        {
          return c(row);
        }
    };

    template<typename T, typename E>
    struct BackendVectorSelectorHelper<EigenVectorBackend,T,E>
    {
        typedef EigenVectorBackend::VectorContainer<T,E> Type;
    };

  } // namespace PDELab
} // namespace Dune

#endif

#endif // DUNE_PDELAB_EIGENVECTORBACKEND_HH

// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
