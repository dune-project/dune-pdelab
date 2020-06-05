/*
 * Modified version of the ISTL Arpack++ wrapper for supporting
 * generalized eigenproblems as required by the GenEO coarse space.
 */

#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_ARPACK_GENEO_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_ARPACK_GENEO_HH

#if HAVE_ARPACKPP

#include <cmath>  // provides std::abs, std::pow, std::sqrt

#include <iostream>  // provides std::cout, std::endl
#include <string>    // provides std::string

#include <algorithm>
#include <numeric>
#include <vector>

#include <dune/common/fvector.hh>     // provides Dune::FieldVector
#include <dune/common/exceptions.hh>  // provides DUNE_THROW(...)

#if DUNE_VERSION_GTE(ISTL,2,8)
#include <dune/istl/blocklevel.hh>
#endif
#include <dune/istl/bvector.hh>        // provides Dune::BlockVector
#include <dune/istl/istlexception.hh>  // provides Dune::ISTLError
#include <dune/istl/io.hh>             // provides Dune::printvector(...)

#include <dune/pdelab/backend/interface.hh>

#ifdef Status
#undef Status        // prevent preprocessor from damaging the ARPACK++
                     // code when "X11/Xlib.h" is included (the latter
                     // defines Status as "#define Status int" and
                     // ARPACK++ provides a class with a method called
                     // Status)
#endif
#include "arssym.h"  // provides ARSymStdEig
#include "argsym.h"  // provides ARSymGenEig
#include "argnsym.h"  // provides ARSymGenEig


namespace ArpackGeneo
{

  /**
   * \brief Wrapper for a DUNE-ISTL BCRSMatrix which can be used
   *        together with some algorithms of the ARPACK++ library
   *        which solve a generalized eigenproblems
   *
   *
   * \tparam BCRSMatrix Type of a DUNE-ISTL BCRSMatrix;
   *                    is assumed to have blocklevel 2.
   *
   */


  template <class BCRSMatrix>
  class ArPackPlusPlus_BCRSMatrixWrapperGen
  {
  public:
    //! Type of the underlying field of the matrix
    typedef typename BCRSMatrix::field_type Real;

    /*!
     * \brief Matrix wrapper DUNE to ARPACK++ for generalized eigenproblems
     */
  public:
    //! Construct from BCRSMatrix A
    ArPackPlusPlus_BCRSMatrixWrapperGen (const BCRSMatrix& A, const BCRSMatrix& B)
      : A_(A), B_(B),
        A_solver(A_),
        a_(A_.M() * mBlock), n_(A_.N() * nBlock)
    {
      // assert that BCRSMatrix type has blocklevel 2
#if DUNE_VERSION_LT(DUNE_ISTL,2,8)
      static_assert
        (BCRSMatrix::blocklevel == 2,
         "Only BCRSMatrices with blocklevel 2 are supported.");
#else
      static_assert
        (Dune::blockLevel<BCRSMatrix>() == 2,
         "Only BCRSMatrices with blocklevel 2 are supported.");
#endif
      // allocate memory for auxiliary block vector objects
      // which are compatible to matrix rows / columns
      domainBlockVector.resize(A_.N());
      rangeBlockVector.resize(A_.M());
    }

    //! Perform matrix-vector product w = A*v
    inline void multMv (Real* v, Real* w)
    {
      // get vector v as an object of appropriate type
      arrayToDomainBlockVector(v,domainBlockVector);

      // perform matrix-vector product
      B_.mv(domainBlockVector,rangeBlockVector);

      Dune::InverseOperatorResult result;
      auto rangeBlockVector2 = rangeBlockVector;
      A_solver.apply(rangeBlockVector2, rangeBlockVector, result);

      // get vector w from object of appropriate type
      rangeBlockVectorToArray(rangeBlockVector2,w);
    };

    inline void multMvB (Real* v, Real* w)
    {
      // get vector v as an object of appropriate type
      arrayToDomainBlockVector(v,domainBlockVector);

      // perform matrix-vector product
      B_.mv(domainBlockVector,rangeBlockVector);

      // get vector w from object of appropriate type
      rangeBlockVectorToArray(rangeBlockVector,w);
    };


    //! Return number of rows in the matrix
    inline int nrows () const { return a_; }

    //! Return number of columns in the matrix
    inline int ncols () const { return n_; }

  protected:
    // Number of rows and columns in each block of the matrix
    constexpr static int mBlock = BCRSMatrix::block_type::rows;
    constexpr static int nBlock = BCRSMatrix::block_type::cols;

    // Type of vectors in the domain of the linear map associated with
    // the matrix, i.e. block vectors compatible to matrix rows
    constexpr static int dbvBlockSize = nBlock;
    typedef Dune::FieldVector<Real,dbvBlockSize> DomainBlockVectorBlock;
    typedef Dune::BlockVector<DomainBlockVectorBlock> DomainBlockVector;

    // Type of vectors in the range of the linear map associated with
    // the matrix, i.e. block vectors compatible to matrix columns
    constexpr static int rbvBlockSize = mBlock;
    typedef Dune::FieldVector<Real,rbvBlockSize> RangeBlockVectorBlock;
    typedef Dune::BlockVector<RangeBlockVectorBlock> RangeBlockVector;

    // Types for vector index access
    typedef typename DomainBlockVector::size_type dbv_size_type;
    typedef typename RangeBlockVector::size_type rbv_size_type;
    typedef typename DomainBlockVectorBlock::size_type dbvb_size_type;
    typedef typename RangeBlockVectorBlock::size_type rbvb_size_type;

    // Get vector v from a block vector object which is compatible to
    // matrix rows
    static inline void
      domainBlockVectorToArray (const DomainBlockVector& dbv, Real* v)
    {
      for (dbv_size_type block = 0; block < dbv.N(); ++block)
        for (dbvb_size_type iBlock = 0; iBlock < dbvBlockSize; ++iBlock)
          v[block*dbvBlockSize + iBlock] = dbv[block][iBlock];
    }

    // Get vector v from a block vector object which is compatible to
    // matrix columns
    static inline void
      rangeBlockVectorToArray (const RangeBlockVector& rbv, Real* v)
    {
      for (rbv_size_type block = 0; block < rbv.N(); ++block)
        for (rbvb_size_type iBlock = 0; iBlock < rbvBlockSize; ++iBlock)
          v[block*rbvBlockSize + iBlock] = rbv[block][iBlock];
    }

  public:
    //! Get vector v as a block vector object which is compatible to
    //! matrix rows
    static inline void arrayToDomainBlockVector (const Real* v,
                                                 DomainBlockVector& dbv)
    {
      for (dbv_size_type block = 0; block < dbv.N(); ++block)
        for (dbvb_size_type iBlock = 0; iBlock < dbvBlockSize; ++iBlock)
          dbv[block][iBlock] = v[block*dbvBlockSize + iBlock];
    }

    template<typename Vector>
    static inline void arrayToVector(const Real* data, Vector& v)
    {
      std::copy(data,data + v.flatsize(),v.begin());
    }

    //! Get vector v as a block vector object which is compatible to
    //! matrix columns
    static inline void arrayToRangeBlockVector (const Real* v,
                                                RangeBlockVector& rbv)
    {
      for (rbv_size_type block = 0; block < rbv.N(); ++block)
        for (rbvb_size_type iBlock = 0; iBlock < rbvBlockSize; ++iBlock)
          rbv[block][iBlock] = v[block*rbvBlockSize + iBlock];
    }

  protected:
    // The DUNE-ISTL BCRSMatrix
    const BCRSMatrix& A_;
    const BCRSMatrix& B_;

    Dune::UMFPack<BCRSMatrix> A_solver;

    // Number of rows and columns in the matrix
    const int a_, n_;

    // Auxiliary block vector objects which are
    // compatible to matrix rows / columns
    mutable DomainBlockVector domainBlockVector;
    mutable RangeBlockVector rangeBlockVector;
  };


  /**
   * \brief A modified version of the corresponding ISTL class supporting
   * specific generalized eigenproblems as required for the GenEO coarse space.
   *
   * \note For a recent version of the ARPACK++ library working with recent
   *       compiler versions see "http://reuter.mit.edu/software/arpackpatch/"
   *       or the git repository "https://github.com/m-reuter/arpackpp.git".
   *
   *
   * \tparam BCRSMatrix  Type of a DUNE-ISTL BCRSMatrix whose eigenvalues
   *                     respectively singular values shall be considered;
   *                     is assumed to have blocklevel 2.
   * \tparam BlockVector Type of the associated vectors; compatible with the
   *                     rows of a BCRSMatrix object (if #rows >= #ncols) or
   *                     its columns (if #rows < #ncols).
   *
   * \author Sebastian Westerheide.
   */
  template <typename BCRSMatrix, typename BlockVectorWrapper>
  class ArPackPlusPlus_Algorithms
  {
  public:

    using BlockVector = Dune::PDELab::Backend::Native<BlockVectorWrapper>;
    typedef typename BlockVector::field_type Real;

  public:
    /**
     * \brief Construct from required parameters.
     *
     * \param[in] m               The DUNE-ISTL BCRSMatrix whose eigenvalues
     *                            resp. singular values shall be considered.
     * \param[in] nIterationsMax  Influences the maximum number of Arnoldi
     *                            update iterations allowed; depending on the
     *                            algorithm, c*nIterationsMax iterations may
     *                            be performed, where c is a natural number.
     * \param[in] verbosity_level Verbosity setting;
     *                            >= 1: algorithms print a preamble and
     *                                  the final result,
     *                            >= 2: algorithms print information about
     *                                  the problem solved using ARPACK++,
     *                            >= 3: the final result output includes
     *                                  the approximated eigenvector,
     *                            >= 4: sets the ARPACK(++) verbosity mode.
     */
    ArPackPlusPlus_Algorithms (const BCRSMatrix& m,
                               const unsigned int nIterationsMax = 100000,
                               const unsigned int verbosity_level = 0)
      : a_(m), nIterationsMax_(nIterationsMax),
        verbosity_level_(verbosity_level),
        nIterations_(0),
        title_("    ArPackPlusPlus_Algorithms: "),
        blank_(title_.length(),' ')
    {}


    inline void computeGenNonSymMinMagnitude (const BCRSMatrix& b_, const Real& epsilon,
                                              std::vector<BlockVectorWrapper>& x, std::vector<Real>& lambda, Real sigma) const
    {
      // print verbosity information
      if (verbosity_level_ > 0)
        std::cout << title_ << "Computing an approximation of the "
                  << "least dominant eigenvalue of a matrix which "
                  << "is assumed to be symmetric." << std::endl;



      // allocate memory for variables, set parameters
      const int nev = x.size();              // Number of eigenvalues to compute
      const int ncv = 0;                     // Number of Arnoldi vectors generated at each iteration (0 == auto)
      const Real tol = epsilon;              // Stopping tolerance (relative accuracy of Ritz values) (0 == machine precision)
      const int maxit = nIterationsMax_*nev; // Maximum number of Arnoldi update iterations allowed   (0 == 100*nev)
      auto ev = std::vector<Real>(nev);      // Computed generalized eigenvalues
      auto ev_imag = std::vector<Real>(nev); // Computed generalized eigenvalues
      const bool ivec = true;                // Flag deciding if eigenvectors shall be determined

      BCRSMatrix ashiftb(a_);
      ashiftb.axpy(-sigma,b_);

      // use type ArPackPlusPlus_BCRSMatrixWrapperGen to store matrix information
      // and to perform the product (A-sigma B)^-1 v (LU decomposition is not used)
      typedef ArPackPlusPlus_BCRSMatrixWrapperGen<BCRSMatrix> WrappedMatrix;
      WrappedMatrix A(ashiftb,b_);

      // get number of rows and columns in A
      const int nrows = A.nrows();
      const int ncols = A.ncols();

      // assert that A is square
      if (nrows != ncols)
        DUNE_THROW(Dune::ISTLError,"Matrix is not square ("
                   << nrows << "x" << ncols << ").");

      // define what we need: eigenvalues with smallest magnitude
      char which[] = "LM";
      //ARNonSymGenEig<Real,WrappedMatrix,WrappedMatrix>
      //  dprob(nrows, nev, &A, &WrappedMatrix::multMv, &A, &WrappedMatrix::multMvB, sigma, which, ncv, tol, maxit);

      //Non generalised problem
      ARNonSymStdEig<Real,WrappedMatrix>
        dprob(nrows, nev, &A, &WrappedMatrix::multMv, which, ncv, tol, maxit);


      // set ARPACK verbosity mode if requested
      if (verbosity_level_ > 3) dprob.Trace();

      // find eigenvalues and eigenvectors of A
      // The broken interface of ARPACK++ actually wants a reference to a pointer...
      auto ev_data = ev.data();
      auto ev_imag_data = ev_imag.data();
      dprob.Eigenvalues(ev_data,ev_imag_data,ivec);

      // Get sorting permutation for un-shifted eigenvalues
      std::vector<int> index(nev, 0);
      std::iota(index.begin(),index.end(),0);

      std::sort(index.begin(), index.end(),
          [&](const int& a, const int& b) {
            return (sigma+1./ev[a] < sigma+1./ev[b]);
          }
      );

      // Unshift eigenpairs
      for (int i = 0; i < nev; i++) {
        lambda[i] = sigma+1./ev[index[i]];
        Real* x_raw = dprob.RawEigenvector(index[i]);
        WrappedMatrix::arrayToVector(x_raw,x[i]);
      }

      // obtain number of Arnoldi update iterations actually taken
      nIterations_ = dprob.GetIter();

    }

    /**
     * \brief Return the number of iterations in last application of
     *        an algorithm.
     */
    inline unsigned int getIterationCount () const
    {
      if (nIterations_ == 0)
        DUNE_THROW(Dune::ISTLError,"No algorithm applied, yet.");

      return nIterations_;
    }

  protected:
    // parameters related to iterative eigenvalue algorithms
    const BCRSMatrix& a_;
    const unsigned int nIterationsMax_;

    // verbosity setting
    const unsigned int verbosity_level_;

    // memory for storing temporary variables (mutable as they shall
    // just be effectless auxiliary variables of the const apply*(...)
    // methods)
    mutable unsigned int nIterations_;

    // constants for printing verbosity information
    const std::string title_;
    const std::string blank_;
  };

}  // namespace Dune


#endif  // HAVE_ARPACKPP

#endif  // DUNE_PDELAB_BACKEND_ISTL_GENEO_ARPACK_GENEO_HH
