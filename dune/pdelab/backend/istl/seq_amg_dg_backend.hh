#ifndef DUNE_PDELAB_SEQ_AMG_DG_BACKEND_HH
#define DUNE_PDELAB_SEQ_AMG_DG_BACKEND_HH

#include <dune/common/power.hh>
#include <dune/common/parametertree.hh>

#include <dune/istl/matrixmatrix.hh>

#include <dune/grid/common/datahandleif.hh>

#include <dune/pdelab/backend/istl/vector.hh>
#include <dune/pdelab/backend/istl/bcrsmatrix.hh>
#include <dune/pdelab/backend/istl/ovlpistlsolverbackend.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>

namespace Dune {
  namespace PDELab {

    /** An ISTL preconditioner for DG based on AMG applied to CG subspace

        The template parameters are:
        DGMatrix    BCRSMatrix assembled with DG
        DGPrec      preconditioner to be used for DG
        CGPrec      preconditioner to be used on CG subspace
        P           BCRSMatrix for grid transfer
    */
    template<class DGMatrix, class DGPrec, class CGPrec, class P>
    class SeqDGAMGPrec : public Dune::Preconditioner<typename DGPrec::domain_type,typename DGPrec::range_type>
    {
      DGMatrix& dgmatrix;
      DGPrec& dgprec;
      CGPrec& cgprec;
      P& p;
      int n1,n2;
      bool firstapply;

    public:
      typedef typename DGPrec::domain_type X;
      typedef typename DGPrec::range_type Y;
      typedef typename CGPrec::domain_type CGX;
      typedef typename CGPrec::range_type CGY;

      // define the category
      enum {
        //! \brief The category the preconditioner is part of.
        category=Dune::SolverCategory::sequential
      };

      /*! \brief Constructor.

        Constructor gets all parameters to operate the prec.
        \param A The matrix to operate on.
        \param n The number of iterations to perform.
        \param w The relaxation factor.
      */
      SeqDGAMGPrec (DGMatrix& dgmatrix_, DGPrec& dgprec_, CGPrec& cgprec_, P& p_, int n1_, int n2_)
        : dgmatrix(dgmatrix_), dgprec(dgprec_), cgprec(cgprec_), p(p_), n1(n1_), n2(n2_),
          firstapply(true)
      {
      }

      /*!
        \brief Prepare the preconditioner.

        \copydoc Preconditioner::pre(X&,Y&)
      */
      virtual void pre (X& x, Y& b)
      {
        dgprec.pre(x,b);
        CGY cgd(p.M());
        cgd = 0.0;
        CGX cgv(p.M());
        cgv = 0.0;
        cgprec.pre(cgv,cgd);
      }

      /*!
        \brief Apply the precondioner.

        \copydoc Preconditioner::apply(X&,const Y&)
      */
      virtual void apply (X& x, const Y& b)
      {
        // need local copies to store defect and solution
        Y d(b);
        X v(x);

        // pre-smoothing on DG matrix
        for (int i=0; i<n1; i++)
          {
            v = 0.0;
            dgprec.apply(v,d);
            dgmatrix.mmv(v,d);
            x += v;
          }

        // restrict defect to CG subspace
        CGY cgd(p.M());
        p.mtv(d,cgd);
        CGX cgv(p.M());
        cgv = 0.0;

        // apply AMG
        cgprec.apply(cgv,cgd);

        // prolongate correction
        p.mv(cgv,v);
        dgmatrix.mmv(v,d);
        x += v;

        // pre-smoothing on DG matrix
        for (int i=0; i<n2; i++)
          {
            v = 0.0;
            dgprec.apply(v,d);
            dgmatrix.mmv(v,d);
            x += v;
          }
      }

      /*!
        \brief Clean up.

        \copydoc Preconditioner::post(X&)
      */
      virtual void post (X& x)
      {
        dgprec.post(x);
        CGX cgv(p.M());
        cgv = 0.0;
        cgprec.post(cgv);
      }
    };


    /** Sequential solver backend for using AMG for DG in PDELab

        The template parameters are:
        DGGO       GridOperator for DG discretization, allows access to matrix, vector and grid function space
        CGGFS      grid function space for CG subspace
        DGPrec     preconditioner for DG problem
        Solver     solver to be used on the complete problem

    */
    template<class DGGO, class CGGFS, class TransferLOP, template<class,class,class,int> class DGPrec, template<class> class Solver>
    class ISTLBackend_SEQ_AMG_4_DG : public Dune::PDELab::LinearResultStorage
    {
      // DG grid function space
      typedef typename DGGO::Traits::TrialGridFunctionSpace GFS;

      // vectors and matrices on DG level
      typedef typename DGGO::Traits::Jacobian M; // wrapped istl DG matrix
      typedef typename DGGO::Traits::Domain V;   // wrapped istl DG vector
      typedef typename M::BaseT Matrix;          // istl DG matrix
      typedef typename V::BaseT Vector;          // istl DG vector
      typedef typename Vector::field_type field_type;

      // vectors and matrices on CG level
      using CGV = Dune::PDELab::Backend::Vector<CGGFS,field_type>; // wrapped istl CG vector
      typedef typename CGV::BaseT CGVector;                               // istl CG vector

      // prolongation matrix
      typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
      typedef Dune::PDELab::EmptyTransformation CC;
      typedef TransferLOP CGTODGLOP; // local operator
      typedef Dune::PDELab::GridOperator<CGGFS,GFS,CGTODGLOP,MBE,field_type,field_type,field_type,CC,CC> PGO;
      typedef typename PGO::Jacobian PMatrix; // wrapped ISTL prolongation matrix
      typedef typename PMatrix::BaseT P;      // ISTL prolongation matrix

      // CG subspace matrix
      typedef typename Dune::TransposedMatMultMatResult<P,Matrix>::type PTADG;
      typedef typename Dune::MatMultMatResult<PTADG,P>::type ACG; // istl coarse space matrix
      typedef ACG CGMatrix; // another name

      DGGO& dggo;
      CGGFS& cggfs;
      unsigned maxiter;
      int verbose;
      bool usesuperlu;
      std::size_t low_order_space_entries_per_row;

      CGTODGLOP cgtodglop;  // local operator to assemble prolongation matrix
      PGO pgo;              // grid operator to assemble prolongation matrix
      PMatrix pmatrix;      // wrapped prolongation matrix

    public:
      ISTLBackend_SEQ_AMG_4_DG(DGGO& dggo_, CGGFS& cggfs_, unsigned maxiter_=5000, int verbose_=1, bool usesuperlu_=true)
        : dggo(dggo_)
        , cggfs(cggfs_)
        , maxiter(maxiter_)
        , verbose(verbose_)
        , usesuperlu(usesuperlu_)
        , low_order_space_entries_per_row(StaticPower<3,GFS::Traits::GridView::dimension>::power)
        , cgtodglop()
        , pgo(cggfs,dggo.trialGridFunctionSpace(),cgtodglop,MBE(low_order_space_entries_per_row))
        , pmatrix(pgo)
      {
#if !HAVE_SUPERLU
        if (usesuperlu == true)
          {
            std::cout << "WARNING: You are using AMG without SuperLU!"
                      << " Please consider installing SuperLU,"
                      << " or set the usesuperlu flag to false"
                      << " to suppress this warning." << std::endl;
          }
#endif
      }

      ISTLBackend_SEQ_AMG_4_DG(DGGO& dggo_, CGGFS& cggfs_, const ParameterTree& params)//unsigned maxiter_=5000, int verbose_=1, bool usesuperlu_=true)
        : dggo(dggo_)
        , cggfs(cggfs_)
        , maxiter(params.get<int>("max_iterations",5000))
        , verbose(params.get<int>("verbose",1))
        , usesuperlu(params.get<bool>("use_superlu",true))
        , low_order_space_entries_per_row(params.get<std::size_t>("low_order_space.entries_per_row",StaticPower<3,GFS::Traits::GridView::dimension>::power))
        , cgtodglop()
        , pgo(cggfs,dggo.trialGridFunctionSpace(),cgtodglop,MBE(low_order_space_entries_per_row))
        , pmatrix(pgo)
      {
#if !HAVE_SUPERLU
        if (usesuperlu == true)
          {
            std::cout << "WARNING: You are using AMG without SuperLU!"
                      << " Please consider installing SuperLU,"
                      << " or set the usesuperlu flag to false"
                      << " to suppress this warning." << std::endl;
          }
#endif


        // assemble prolongation matrix; this will not change from one apply to the next
        pmatrix = 0.0;
        if (verbose>0) std::cout << "allocated prolongation matrix of size " << pmatrix.N() << " x " << pmatrix.M() << std::endl;
        CGV cgx(cggfs,0.0);         // need vector to call jacobian
        pgo.jacobian(cgx,pmatrix);
      }

      /*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      typename V::ElementType norm (const V& v) const
      {
        return Backend::native(v).two_norm();
      }

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      void apply (M& A, V& z, V& r, typename V::ElementType reduction)
      {
        using Backend::native;
        // do triple matrix product ACG = P^T ADG P
        Dune::Timer watch;
        watch.reset();
        ACG acg;
        {
          PTADG ptadg;
          Dune::transposeMatMultMat(ptadg,native(pmatrix),native(A));
          Dune::matMultMat(acg,ptadg,native(pmatrix));
        }
        double triple_product_time = watch.elapsed();
        if (verbose>0) std::cout << "=== triple matrix product " << triple_product_time << " s" << std::endl;
        //Dune::printmatrix(std::cout,acg,"triple product matrix","row",10,2);

        // set up AMG solver for the CG subspace
        typedef ACG CGMatrix;
        typedef Dune::MatrixAdapter<CGMatrix,CGVector,CGVector> CGOperator;
        CGOperator cgop(acg);
        typedef Dune::Amg::Parameters Parameters; // AMG parameters (might be nice to change from outside)
        Parameters params(15,2000);
        params.setDefaultValuesIsotropic(CGGFS::Traits::GridViewType::Traits::Grid::dimension);
        params.setDebugLevel(verbose);
        params.setCoarsenTarget(1000);
        params.setMaxLevel(20);
        params.setProlongationDampingFactor(1.8);
        params.setNoPreSmoothSteps(2);
        params.setNoPostSmoothSteps(2);
        params.setGamma(1);
        params.setAdditive(false);
        typedef Dune::SeqSSOR<CGMatrix,CGVector,CGVector,1> Smoother;
        typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
        SmootherArgs smootherArgs;
        smootherArgs.iterations = 2;
        smootherArgs.relaxationFactor = 1.0;
        typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<CGMatrix,Dune::Amg::FirstDiagonal> > Criterion;
        Criterion criterion(params);
        typedef Dune::Amg::AMG<CGOperator,CGVector,Smoother> AMG;
        watch.reset();
        AMG amg(cgop,criterion,smootherArgs);
        double amg_setup_time = watch.elapsed();
        if (verbose>0) std::cout << "=== AMG setup " <<amg_setup_time << " s" << std::endl;

        // set up hybrid DG/CG preconditioner
        Dune::MatrixAdapter<Matrix,Vector,Vector> op(native(A));
        DGPrec<Matrix,Vector,Vector,1> dgprec(native(A),1,1);
        typedef SeqDGAMGPrec<Matrix,DGPrec<Matrix,Vector,Vector,1>,AMG,P> HybridPrec;
        HybridPrec hybridprec(native(A),dgprec,amg,native(pmatrix),2,2);

        // set up solver
        Solver<Vector> solver(op,hybridprec,reduction,maxiter,verbose);

        // solve
        Dune::InverseOperatorResult stat;
        watch.reset();
        solver.apply(native(z),native(r),stat);
        double amg_solve_time = watch.elapsed();
        if (verbose>0) std::cout << "=== Hybrid total solve time " << amg_solve_time+amg_setup_time+triple_product_time << " s" << std::endl;
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = amg_solve_time+amg_setup_time+triple_product_time;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }

    };
  }
}
#endif
