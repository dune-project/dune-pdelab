// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_BACKEND_SOLVERWRAPPER_HH
#define DUNE_PDELAB_BACKEND_SOLVERWRAPPER_HH

// This file implements a solver wrapper that does the following:
//
// It stores a linear solver backend
//
// The wrapper has a matrix-based and a matrix-free apply function. It uses
// SFINAE to detect if the linear solver backend implements those methods. If
// it does the call is forwarded to the actual backend. If it does not exist
// the wrapper throws an error (at runtime).
//
// The SFINAE check happens at compile time without throwing an error if the
// method in question does not exist.
//
// What was the purpose of this class: We want to have a
// StationaryLinearProblemSolver and a Newton implementation that works for
// both matrix-based and matrix-free solvers. This could also be solved in
// other ways (e.g. having a base class defining the solver backend interface)
// but this is the least invasive solution.

namespace Dune{
  namespace PDELab{

#ifndef DOXYGEN
    namespace impl{
      // SFINAE checks that test if the apply methods or the isMatrixFree
      // method exist. Note: In the apply you need to pass all the required
      // arguments. The 42 is passed as reduction. You could also pass
      // something like std::declval<int>() or something else.
      template<typename LS, typename M, typename V, typename W, typename = void>
      struct linearSolverHasMatrixBasedApply : std::false_type {};

      template<typename LS, typename M, typename V, typename W>
      struct linearSolverHasMatrixBasedApply<LS, M, V, W, std::void_t<decltype(std::declval<LS>().apply(std::declval<M&>(), std::declval<V&>(), std::declval<W&>(), 42))>> : std::true_type {};

      template<typename LS, typename V, typename W, typename = void>
      struct linearSolverHasMatrixFreeApply : std::false_type {};

      template<typename LS, typename V, typename W>
      struct linearSolverHasMatrixFreeApply<LS, V, W, std::void_t<decltype(std::declval<LS>().apply(std::declval<V&>(), std::declval<W&>(), 42))>> : std::true_type {};

      template<typename LS, typename = void>
      struct linearSolverHasIsMatrixFree : std::false_type {};

      template<typename LS>
      struct linearSolverHasIsMatrixFree<LS, std::void_t<decltype(std::declval<LS>().isMatrixFree())>> : std::true_type {};


      // A wrapper for linear solver backends, that checks if the linear solver
      // backend implements certain methods. If it does the call gets
      // forwarded, otherwise we throw an error.
      template <typename LS>
      struct LinearSolverWrapper{
        LS& _ls;

        LinearSolverWrapper (LS& ls) : _ls(ls) {}

        // Matrix based apply method
        template<class M, class V, class W>
        void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename W::ElementType >::real_type reduction)
        {
          if constexpr (linearSolverHasMatrixBasedApply<LS, M, V, W>::value){
            _ls.apply(A, z, r, reduction);
          }
          else{
            DUNE_THROW(Exception, "Your linear solver backend doesn't have a matrix-based apply method.");
          }
        }

        // Matrix free apply method
        template<class V, class W>
        void apply(V& z, W& r, typename Dune::template FieldTraits<typename W::ElementType >::real_type reduction)
        {
          if constexpr (linearSolverHasMatrixFreeApply<LS, V, W>::value){
            _ls.apply(z, r, reduction);
          }
          else{
            DUNE_THROW(Exception, "Your linear solver backend doesn't have a matrix-free apply method.");
          }
        }

        // Note: If the linear solver backend does not have an isMatrixFree
        // method we return false.
        bool isMatrixFree()
        {
          if constexpr (linearSolverHasIsMatrixFree<LS>::value){
            return _ls.isMatrixFree();
          }
          else{
            return false;
          }
        }
      };
    } // namespace impl
#endif
  } // namespace PDELab
}  // namespace Dune
#endif
