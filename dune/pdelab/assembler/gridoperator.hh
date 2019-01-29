#ifndef DUNE_PDELAB_ASSEMBLER_GRIDOPERATOR_HH
#define DUNE_PDELAB_ASSEMBLER_GRIDOPERATOR_HH

#include <tuple>

#include <dune/pdelab/common/typetraits.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperator/common/borderdofexchanger.hh>
#include <dune/pdelab/assembler/assembler.hh>
#include <dune/pdelab/assembler/residualengine.hh>
#include <dune/pdelab/assembler/patternengine.hh>
#include <dune/pdelab/assembler/jacobianengine.hh>
#include <dune/pdelab/assembler/applyjacobianengine.hh>


namespace Dune{
  namespace PDELab{

    /**
       \brief Standard grid operator implementation

       \tparam GFSU GridFunctionSpace for ansatz functions
       \tparam GFSV GridFunctionSpace for test functions
       \tparam MB The matrix backend to be used for representation of the jacobian
       \tparam DF The domain field type of the operator
       \tparam RF The range field type of the operator
       \tparam JF The jacobian field type
       \tparam CU   Constraints maps for the individual dofs (trial space)
       \tparam CV   Constraints maps for the individual dofs (test space)

    */
    template<typename GFSU, typename GFSV, typename LOP,
             typename MB, typename DF, typename RF, typename JF,
             typename CU=Dune::PDELab::EmptyTransformation,
             typename CV=Dune::PDELab::EmptyTransformation,
             bool instationary_ = false
             >
    class NewGridOperator
      : OnlyMovable
    {
    public:

      using EntitySet = typename GFSU::Traits::EntitySet;

      //! The global assembler type
      using Assembler = Dune::PDELab::Assembler<EntitySet>;

      //! The type of the domain (solution).
      using Domain = Dune::PDELab::Backend::Vector<GFSU,DF>;
      //! The type of the range (residual).
      using Range = Dune::PDELab::Backend::Vector<GFSV,RF>;
      //! The type of the jacobian.
      using Jacobian = Dune::PDELab::Backend::Matrix<MB,Domain,Range,JF>;

      //! The sparsity pattern container for the jacobian matrix
      using Pattern = typename MB::template Pattern<Jacobian,GFSV,GFSU>;

      // Fix this as soon as the default Partitions are constexpr
      typedef typename std::conditional<
        GFSU::Traits::EntitySet::Partitions::partitionIterator() == InteriorBorder_Partition,
        NonOverlappingBorderDOFExchanger<NewGridOperator>,
        OverlappingBorderDOFExchanger<NewGridOperator>
        >::type BorderDOFExchanger;

      //! The grid operator traits
      struct Traits
      {

        //! The trial grid function space.
        typedef GFSU TrialGridFunctionSpace;

        //! The test grid function space.
        typedef GFSV TestGridFunctionSpace;

        using EntitySet = typename GFSU::Traits::EntitySet;


        //! The type of the trial grid function space constraints.
        typedef CU TrialGridFunctionSpaceConstraints;

        //! The type of the test grid function space constraints.
        typedef CV TestGridFunctionSpaceConstraints;


        //! The matrix backend of the grid operator.
        typedef MB MatrixBackend;


        //! The field type of the domain (solution).
        typedef DF DomainField;

        //! The type of the domain (solution).
        using Domain = Dune::PDELab::Backend::Vector<GFSU,DF>;


        //! The field type of the range (residual).
        typedef RF RangeField;

        //! The type of the range (residual).
        using Range = Dune::PDELab::Backend::Vector<GFSV,RF>;


        //! The field type of the jacobian.
        typedef JF JacobianField;

        //! The type of the jacobian.
        using Jacobian = Dune::PDELab::Backend::Matrix<MB,Domain,Range,JF>;

        using LocalOperator = LOP;


        //! The global assembler of the grid operator.
        using Assembler = Dune::PDELab::Assembler<EntitySet>;

        using PatternEngine = Dune::PDELab::PatternEngine<
          GFSU,
          GFSV,
          LOP,
          MB,
          JF,
          CU,
          CV
          >;

        using ResidualEngine = Dune::PDELab::ResidualEngine<
          Domain,
          Range,
          LOP,
          CU,
          CV,
          instationary_
          >;

        using JacobianEngine = Dune::PDELab::JacobianEngine<
          Domain,
          Jacobian,
          LOP,
          CU,
          CV,
          instationary_
          >;

        using ApplyJacobianEngine = Dune::PDELab::ApplyJacobianEngine<
          Domain,
          Range,
          LOP,
          CU,
          CV,
          instationary_
          >;

      };

      using Real = typename Traits::DomainField;

      template <typename MFT>
      struct MatrixContainer{
        typedef typename Traits::Jacobian Type;
      };

      //! Constructor for non trivial constraints
      NewGridOperator(const GFSU& gfsu, CU& cu, const GFSV& gfsv, CV& cv, LOP& lop, const MB& mb = MB())
        : _assembler(gfsu.entitySet())
        , _gfsu(&gfsu)
        , _gfsv(&gfsv)
        , _lop(&lop)
        , _cu(cu)
        , _cv(cv)
        , _dof_exchanger(std::make_shared<BorderDOFExchanger>(*this))
        , _backend(mb)
      {}

      //! Constructor for empty constraints
      NewGridOperator(const GFSU& gfsu, const GFSV& gfsv, LOP& lop, const MB& mb = MB())
        : _assembler(gfsu.entitySet())
        , _gfsu(&gfsu)
        , _gfsv(&gfsv)
        , _lop(&lop)
        , _dof_exchanger(std::make_shared<BorderDOFExchanger>(*this))
        , _backend(mb)
      {}

      //! Get the trial grid function space
      const GFSU& trialGridFunctionSpace() const
      {
        return *_gfsu;
      }

      //! Get the test grid function space
      const GFSV& testGridFunctionSpace() const
      {
        return *_gfsv;
      }

      const CU& trialConstraints() const
      {
        return _cu;
      }

      CU& trialConstraints()
      {
        return _cu;
      }

      const CV& testConstraints() const
      {
        return _cv;
      }

      CV& testConstraints()
      {
        return _cv;
      }

      LOP& localOperator() const
      {
        return *_lop;
      }
/*
      const LOP& localOperator() const
      {
        return *_lop;
      }
*/
      //! Get dimension of space u
      typename GFSU::Traits::SizeType globalSizeU () const
      {
        return trialGridFunctionSpace().globalSize();
      }

      //! Get dimension of space v
      typename GFSV::Traits::SizeType globalSizeV () const
      {
        return testGridFunctionSpace().globalSize();
      }

      Assembler & assembler() { return _assembler; }

      const Assembler & assembler() const { return _assembler; }

      //! Interpolate the constrained dofs from given function
      template<typename F, typename X>
      void interpolate (const X& xold, F& f, X& x) const
      {
        // Interpolate f into grid function space and set corresponding coefficients
        Dune::PDELab::interpolate(f,trialGridFunctionSpace(),x);

        // Copy non-constrained dofs from old time step
        Dune::PDELab::copy_nonconstrained_dofs(trialConstraints(),xold,x);
      }

      auto makeMatrixPatternEngine() const
      {
        return std::make_shared<typename Traits::PatternEngine>(
          trialGridFunctionSpace(),
          testGridFunctionSpace(),
          localOperator(),
          matrixBackend(),
          trialConstraints(),
          testConstraints(),
          JF()
          );
      }

      auto matrixPatternEngine() const
      {
        if (not _matrix_pattern_engine)
          _matrix_pattern_engine = makeMatrixPatternEngine();
        return _matrix_pattern_engine;
      }

      //! Fill pattern of jacobian matrix
      auto matrixPattern() const
      {
        return assembler().assemble(*matrixPatternEngine());
      }

      auto makeResidualEngine() const
      {
        return std::make_shared<typename Traits::ResidualEngine>(
          localOperator(),
          trialConstraints(),
          testConstraints()
          );
      }

      auto residualEngine() const
      {
        if (not _residual_engine)
          _residual_engine = makeResidualEngine();
        return _residual_engine;
      }

      //! Assemble residual
      void residual(const Domain & x, Range & r) const
      {
        auto engine = residualEngine();
        engine->setArgument(x);
        engine->setResidual(r);
        assembler().assemble(*engine);
      }

      auto makeJacobianEngine() const
      {
        return std::make_shared<typename Traits::JacobianEngine>(
          localOperator(),
          trialConstraints(),
          testConstraints()
          );
      }

      auto jacobianEngine() const
      {
        if (not _jacobian_engine)
          _jacobian_engine = makeJacobianEngine();
        return _jacobian_engine;
      }

      //! Assembler jacobian
      void jacobian(const Domain & x, Jacobian & a) const
      {
        auto engine = jacobianEngine();
        engine->setLinearizationPoint(x);
        engine->setJacobian(a);
        assembler().assemble(*engine);
      }

      auto makeApplyJacobianEngine() const
      {
        return std::make_shared<typename Traits::ApplyJacobianEngine>(
          localOperator(),
          trialConstraints(),
          testConstraints()
          );
      }

      auto applyJacobianEngine() const
      {
        if (not _apply_jacobian_engine)
          _apply_jacobian_engine = makeApplyJacobianEngine();
        return _apply_jacobian_engine;
      }

      //! Apply jacobian matrix without explicitly assembling it
      void jacobian_apply(const Domain & z, Range & r) const
      {
        auto engine = applyJacobianEngine();
        engine->setArgument(z);
        engine->setResult(r);
        assembler().assemble(*engine);
      }

      //! Apply jacobian matrix without explicitly assembling it
      void nonlinear_jacobian_apply(const Domain & x, const Domain & z, Range & r) const
      {
        auto engine = applyJacobianEngine();
        engine->setArgument(x);
        engine->setLinearizationPoint(z);
        engine->setResult(r);
        assembler().assemble(*engine);
      }

      void make_consistent(Jacobian& a) const
      {
        _dof_exchanger->accumulateBorderEntries(*this,a);
      }

      void update()
      {
        // the DOF exchanger has matrix information, so we need to update it
        _dof_exchanger->update(*this);
      }

      //! Get the matrix backend for this grid operator.
      const typename Traits::MatrixBackend& matrixBackend() const
      {
        return _backend;
      }

      void setOneStepMethod(std::shared_ptr<OneStep::Method<typename Traits::RangeField>> osm)
      {
        residualEngine()->setOneStepMethod(osm);
        if (_jacobian_engine)
          _jacobian_engine->setOneStepMethod(osm);
        if (_apply_jacobian_engine)
          _apply_jacobian_engine->setOneStepMethod(osm);
      }

      int startStep(Real t0, Real dt)
      {
        int stages = _residual_engine()->startStep(t0,dt);
        if (_jacobian_engine)
          _jacobian_engine->startStep(t0,dt);
        if (_apply_jacobian_engine)
          _apply_jacobian_engine->startStep(t0,dt);
        return stages;
      }

      bool acceptStage(int stage, const typename Traits::Domain& solution)
      {
        bool changed = residualEngine()->acceptStage(stage,assembler(),solution);
        if (_jacobian_engine)
          changed |= _jacobian_engine->acceptStage(stage,assembler(),solution);
        if (_apply_jacobian_engine)
          changed |= _apply_jacobian_engine->acceptStage(stage,assembler(),solution);
        return changed;
      }


    private:
      Assembler _assembler;

      const GFSU* _gfsu;
      const GFSV* _gfsv;

      std::conditional_t<
        std::is_same<CU,EmptyTransformation>::value,
        CU,
        CU&
        > _cu;
      std::conditional_t<
        std::is_same<CV,EmptyTransformation>::value,
        CV,
        CV&
        > _cv;

      mutable LOP* _lop;

      std::shared_ptr<BorderDOFExchanger> _dof_exchanger;

      MB _backend;

      mutable std::shared_ptr<typename Traits::PatternEngine> _matrix_pattern_engine;
      mutable std::shared_ptr<typename Traits::ResidualEngine> _residual_engine;
      mutable std::shared_ptr<typename Traits::JacobianEngine> _jacobian_engine;
      mutable std::shared_ptr<typename Traits::ApplyJacobianEngine> _apply_jacobian_engine;

    };

  }
}
#endif // DUNE_PDELAB_ASSEMBLER_GRIDOPERATOR_HH
