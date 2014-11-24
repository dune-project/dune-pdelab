#ifndef DUNE_PDELAB_GRIDOPERATOR_TBB_HH
#define DUNE_PDELAB_GRIDOPERATOR_TBB_HH

#include <mutex>

#include <dune/common/tupleutility.hh>

#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperator/common/borderdofexchanger.hh>
#include <dune/pdelab/gridoperator/common/gridoperatorutilities.hh>
#include <dune/pdelab/gridoperator/tbb/assembler.hh>
#include <dune/pdelab/gridoperator/tbb/localassembler.hh>

namespace Dune{
  namespace PDELab{

    /**
       \brief Grid operator implementation using tbb parallel algorithms

       \tparam Partitioning Partitioning to use for parallel traversal.
       \tparam GFSU GridFunctionSpace for ansatz functions
       \tparam GFSV GridFunctionSpace for test functions
       \tparam MB The matrix backend to be used for representation of the jacobian
       \tparam DF The domain field type of the operator
       \tparam RF The range field type of the operator
       \tparam JF The jacobian field type
       \tparam CU   Constraints maps for the individual dofs (trial space)
       \tparam CV   Constraints maps for the individual dofs (test space)
       \tparam nonoverlapping_mode Switch for nonoverlapping grids

    */
    template<typename Partitioning, typename GFSU, typename GFSV, typename LOP,
             typename MB, typename DF, typename RF, typename JF,
             typename LM,
             typename CU=Dune::PDELab::EmptyTransformation,
             typename CV=Dune::PDELab::EmptyTransformation,
             bool nonoverlapping_mode = false>
    class TBBGridOperator
    {
    public:

      //! The global assembler type
      typedef TBBAssembler<Partitioning, GFSU,GFSV,CU,CV,nonoverlapping_mode>
        Assembler;

      //! The type of the domain (solution).
      typedef typename Dune::PDELab::BackendVectorSelector<GFSU,DF>::Type Domain;
      //! The type of the range (residual).
      typedef typename Dune::PDELab::BackendVectorSelector<GFSV,RF>::Type Range;
      //! The type of the jacobian.
      typedef typename Dune::PDELab::BackendMatrixSelector<MB,Domain,Range,JF>::Type Jacobian;

      //! The sparsity pattern container for the jacobian matrix
      typedef typename MB::template Pattern<Jacobian,GFSV,GFSU> Pattern;

      //! The local assembler type
      typedef TBBLocalAssembler<TBBGridOperator,LOP,nonoverlapping_mode>
      LocalAssembler;

      typedef typename conditional<
        nonoverlapping_mode,
        NonOverlappingBorderDOFExchanger<TBBGridOperator>,
        OverlappingBorderDOFExchanger<TBBGridOperator>
        >::type BorderDOFExchanger;

      //! The grid operator traits
      typedef Dune::PDELab::GridOperatorTraits
      <GFSU,GFSV,MB,DF,RF,JF,CU,CV,Assembler,LocalAssembler> Traits;

      template <typename MFT>
      struct MatrixContainer{
        typedef typename Traits::Jacobian Type;
      };

      typedef LM LockManager;

      //! Constructor for non trivial constraints
      TBBGridOperator(const GFSU & gfsu_, const CU & cu_, const GFSV & gfsv_,
                      const CV & cv_, LOP & lop_,
                      const shared_ptr<LockManager> &lock_manager_,
                      const MB& mb_ = MB())
        : global_assembler(gfsu_,gfsv_,cu_,cv_)
        , dof_exchanger(make_shared<BorderDOFExchanger>(*this))
        , local_assembler(lop_, cu_, cv_,dof_exchanger, lock_manager_)
        , backend(mb_)
      {}

      //! Constructor for empty constraints
      TBBGridOperator(const GFSU & gfsu_, const GFSV & gfsv_, LOP & lop_,
                      const shared_ptr<LockManager> &lock_manager_,
                      const MB& mb_ = MB())
        : global_assembler(gfsu_,gfsv_)
        , dof_exchanger(make_shared<BorderDOFExchanger>(*this))
        , local_assembler(lop_,dof_exchanger, lock_manager_)
        , backend(mb_)
      {}

      //! Get the trial grid function space
      const GFSU& trialGridFunctionSpace() const
      {
        return global_assembler.trialGridFunctionSpace();
      }

      //! Get the test grid function space
      const GFSV& testGridFunctionSpace() const
      {
        return global_assembler.testGridFunctionSpace();
      }

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

      Assembler & assembler() { return global_assembler; }

      const Assembler & assembler() const { return global_assembler; }

      LocalAssembler & localAssembler() const { return local_assembler; }


      //! Visitor which is called in the method setupGridOperators for
      //! each tuple element.
      template <typename GridOperatorTuple>
      struct SetupGridOperator {
        SetupGridOperator()
          : index(0), size(Dune::tuple_size<GridOperatorTuple>::value) {}

        template <typename T>
        void visit(T& elem) {
          elem.localAssembler().doPreProcessing = index == 0;
          elem.localAssembler().doPostProcessing = index == size-1;
          ++index;
        }

        int index;
        const int size;
      };

      //! Method to set up a number of grid operators which are used
      //! in a joint assembling. It is assumed that all operators are
      //! specializations of the same template type
      template<typename GridOperatorTuple>
      static void setupGridOperators(GridOperatorTuple tuple)
      {
        Dune::ForEachValue<GridOperatorTuple> forEach(tuple);
        SetupGridOperator<GridOperatorTuple> setup_visitor;
        forEach.apply(setup_visitor);
      }

      //! Interpolate the constrained dofs from given function
      template<typename F, typename X>
      void interpolate (const X& xold, F& f, X& x) const
      {
        // Interpolate f into grid function space and set corresponding coefficients
        Dune::PDELab::interpolate(f,global_assembler.trialGridFunctionSpace(),x);

        // Copy non-constrained dofs from old time step
        Dune::PDELab::copy_nonconstrained_dofs(local_assembler.trialConstraints(),xold,x);
      }

      //! Fill pattern of jacobian matrix
      void fill_pattern(Pattern & p) const {
        typedef typename LocalAssembler::LocalPatternAssemblerEngine PatternEngine;
        global_assembler.assemble
          ([&,p](LocalAssembler &la) mutable -> PatternEngine&
                { return la.localPatternAssemblerEngine(p); },
           local_assembler);
      }

      //! Assemble residual
      void residual(const Domain & x, Range & r) const {
        typedef typename LocalAssembler::LocalResidualAssemblerEngine ResidualEngine;
        global_assembler.assemble
          ([&,r,x](LocalAssembler &la) mutable -> ResidualEngine&
                { return la.localResidualAssemblerEngine(r,x); },
           local_assembler);
      }

      //! Assembler jacobian
      void jacobian(const Domain & x, Jacobian & a) const {
        typedef typename LocalAssembler::LocalJacobianAssemblerEngine JacobianEngine;
        global_assembler.assemble
          ([&,a,x](LocalAssembler &la) mutable -> JacobianEngine&
                { return la.localJacobianAssemblerEngine(a,x); },
           local_assembler);
      }

      //! Apply jacobian matrix without explicitly assembling it
      void jacobian_apply(const Domain & x, Range & r) const {
        typedef typename LocalAssembler::LocalJacobianApplyAssemblerEngine JacobianApplyEngine;
        global_assembler.assemble
          ([&,r,x](LocalAssembler &la) -> JacobianApplyEngine&
                { return la.localJacobianApplyAssemblerEngine(r,x); },
           local_assembler);
      }

      void make_consistent(Jacobian& a) const {
        dof_exchanger->accumulateBorderEntries(*this,a);
      }

      void update()
      {
        // the DOF exchanger has matrix information, so we need to update it
        dof_exchanger->update(*this);
      }

      //! Get the matrix backend for this grid operator.
      const typename Traits::MatrixBackend& matrixBackend() const
      {
        return backend;
      }

    private:
      Assembler global_assembler;
      shared_ptr<BorderDOFExchanger> dof_exchanger;

      mutable LocalAssembler local_assembler;
      MB backend;
    };

    //! \brief Grid operator implementation using tbb parallel algorithms,
    //!        with coloring
    /**
     * \tparam Coloring Coloring to use for parallel traversal.
     * \tparam Partitioning Partitioning to use for parallel traversal.
     * \tparam GFSU GridFunctionSpace for ansatz functions
     * \tparam GFSV GridFunctionSpace for test functions
     * \tparam MB The matrix backend to be used for representation of the jacobian
     * \tparam DF The domain field type of the operator
     * \tparam RF The range field type of the operator
     * \tparam JF The jacobian field type
     * \tparam CU   Constraints maps for the individual dofs (trial space)
     * \tparam CV   Constraints maps for the individual dofs (test space)
     * \tparam nonoverlapping_mode Switch for nonoverlapping grids
     */
    template<typename Coloring, typename Partitioning, typename GFSU,
             typename GFSV, typename LOP, typename MB, typename DF,
             typename RF, typename JF,
             typename CU=Dune::PDELab::EmptyTransformation,
             typename CV=Dune::PDELab::EmptyTransformation,
             bool nonoverlapping_mode = false>
    class ColoredTBBGridOperator
    {
    public:

      //! The global assembler type
      typedef ColoredTBBAssembler<Coloring, Partitioning, GFSU,GFSV,CU,CV,
                                  nonoverlapping_mode> Assembler;

      //! The type of the domain (solution).
      typedef typename Dune::PDELab::BackendVectorSelector<GFSU,DF>::Type Domain;
      //! The type of the range (residual).
      typedef typename Dune::PDELab::BackendVectorSelector<GFSV,RF>::Type Range;
      //! The type of the jacobian.
      typedef typename Dune::PDELab::BackendMatrixSelector<MB,Domain,Range,JF>::Type Jacobian;

      //! The sparsity pattern container for the jacobian matrix
      typedef typename MB::template Pattern<Jacobian,GFSV,GFSU> Pattern;

      //! The local assembler type
      typedef ColoredTBBLocalAssembler<ColoredTBBGridOperator,LOP,
                                       nonoverlapping_mode> LocalAssembler;

      typedef typename conditional<
        nonoverlapping_mode,
        NonOverlappingBorderDOFExchanger<ColoredTBBGridOperator>,
        OverlappingBorderDOFExchanger<ColoredTBBGridOperator>
        >::type BorderDOFExchanger;

      //! The grid operator traits
      typedef Dune::PDELab::GridOperatorTraits
      <GFSU,GFSV,MB,DF,RF,JF,CU,CV,Assembler,LocalAssembler> Traits;

      template <typename MFT>
      struct MatrixContainer{
        typedef typename Traits::Jacobian Type;
      };

      //! Constructor for non trivial constraints
      ColoredTBBGridOperator(const GFSU & gfsu_, const CU & cu_,
                             const GFSV & gfsv_, const CV & cv_, LOP & lop_,
                             const MB& mb_ = MB())
        : global_assembler(gfsu_,gfsv_,cu_,cv_)
        , dof_exchanger(make_shared<BorderDOFExchanger>(*this))
        , local_assembler(lop_, cu_, cv_,dof_exchanger)
        , backend(mb_)
      {}

      //! Constructor for empty constraints
      ColoredTBBGridOperator(const GFSU & gfsu_, const GFSV & gfsv_,
                             LOP & lop_, const MB& mb_ = MB())
        : global_assembler(gfsu_,gfsv_)
        , dof_exchanger(make_shared<BorderDOFExchanger>(*this))
        , local_assembler(lop_,dof_exchanger)
        , backend(mb_)
      {}

      //! Get the trial grid function space
      const GFSU& trialGridFunctionSpace() const
      {
        return global_assembler.trialGridFunctionSpace();
      }

      //! Get the test grid function space
      const GFSV& testGridFunctionSpace() const
      {
        return global_assembler.testGridFunctionSpace();
      }

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

      Assembler & assembler() { return global_assembler; }

      const Assembler & assembler() const { return global_assembler; }

      LocalAssembler & localAssembler() const { return local_assembler; }


      //! Visitor which is called in the method setupGridOperators for
      //! each tuple element.
      template <typename GridOperatorTuple>
      struct SetupGridOperator {
        SetupGridOperator()
          : index(0), size(Dune::tuple_size<GridOperatorTuple>::value) {}

        template <typename T>
        void visit(T& elem) {
          elem.localAssembler().doPreProcessing = index == 0;
          elem.localAssembler().doPostProcessing = index == size-1;
          ++index;
        }

        int index;
        const int size;
      };

      //! Method to set up a number of grid operators which are used
      //! in a joint assembling. It is assumed that all operators are
      //! specializations of the same template type
      template<typename GridOperatorTuple>
      static void setupGridOperators(GridOperatorTuple tuple)
      {
        Dune::ForEachValue<GridOperatorTuple> forEach(tuple);
        SetupGridOperator<GridOperatorTuple> setup_visitor;
        forEach.apply(setup_visitor);
      }

      //! Interpolate the constrained dofs from given function
      template<typename F, typename X>
      void interpolate (const X& xold, F& f, X& x) const
      {
        // Interpolate f into grid function space and set corresponding coefficients
        Dune::PDELab::interpolate(f,global_assembler.trialGridFunctionSpace(),x);

        // Copy non-constrained dofs from old time step
        Dune::PDELab::copy_nonconstrained_dofs(local_assembler.trialConstraints(),xold,x);
      }

      //! Fill pattern of jacobian matrix
      void fill_pattern(Pattern & p) const {
        typedef typename LocalAssembler::LocalPatternAssemblerEngine PatternEngine;
        global_assembler.assemble
          ([&,p](LocalAssembler &la) mutable -> PatternEngine&
                { return la.localPatternAssemblerEngine(p); },
           local_assembler);
      }

      //! Assemble residual
      void residual(const Domain & x, Range & r) const {
        typedef typename LocalAssembler::LocalResidualAssemblerEngine ResidualEngine;
        global_assembler.assemble
          ([&,r,x](LocalAssembler &la) mutable -> ResidualEngine&
                { return la.localResidualAssemblerEngine(r,x); },
           local_assembler);
      }

      //! Assembler jacobian
      void jacobian(const Domain & x, Jacobian & a) const {
        typedef typename LocalAssembler::LocalJacobianAssemblerEngine JacobianEngine;
        global_assembler.assemble
          ([&,a,x](LocalAssembler &la) mutable -> JacobianEngine&
                { return la.localJacobianAssemblerEngine(a,x); },
           local_assembler);
      }

      //! Apply jacobian matrix without explicitly assembling it
      void jacobian_apply(const Domain & x, Range & r) const {
        typedef typename LocalAssembler::LocalJacobianApplyAssemblerEngine JacobianApplyEngine;
        global_assembler.assemble
          ([&,r,x](LocalAssembler &la) -> JacobianApplyEngine&
                { return la.localJacobianApplyAssemblerEngine(r,x); },
           local_assembler);
      }

      void make_consistent(Jacobian& a) const {
        dof_exchanger->accumulateBorderEntries(*this,a);
      }

      void update()
      {
        // the DOF exchanger has matrix information, so we need to update it
        dof_exchanger->update(*this);
      }

      //! Get the matrix backend for this grid operator.
      const typename Traits::MatrixBackend& matrixBackend() const
      {
        return backend;
      }

    private:
      Assembler global_assembler;
      shared_ptr<BorderDOFExchanger> dof_exchanger;

      mutable LocalAssembler local_assembler;
      MB backend;

    };

    /**
       \brief Grid operator implementation using tbb parallel algorithms

       \tparam Partitioning Partitioning to use for parallel traversal.
       \tparam GFSU GridFunctionSpace for ansatz functions
       \tparam GFSV GridFunctionSpace for test functions
       \tparam MB The matrix backend to be used for representation of the jacobian
       \tparam DF The domain field type of the operator
       \tparam RF The range field type of the operator
       \tparam JF The jacobian field type
       \tparam CU   Constraints maps for the individual dofs (trial space)
       \tparam CV   Constraints maps for the individual dofs (test space)
       \tparam nonoverlapping_mode Switch for nonoverlapping grids

    */
    template<typename Partitioning, typename GFSU, typename GFSV, typename LOP,
             typename MB, typename DF, typename RF, typename JF,
             typename CU=Dune::PDELab::EmptyTransformation,
             typename CV=Dune::PDELab::EmptyTransformation,
             bool nonoverlapping_mode = false>
    class BatchedTBBGridOperator
    {
    public:

      typedef std::mutex Mutex;

      //! The global assembler type
      typedef TBBAssembler<Partitioning, GFSU,GFSV,CU,CV,nonoverlapping_mode>
        Assembler;

      //! The type of the domain (solution).
      typedef typename Dune::PDELab::BackendVectorSelector<GFSU,DF>::Type Domain;
      //! The type of the range (residual).
      typedef typename Dune::PDELab::BackendVectorSelector<GFSV,RF>::Type Range;
      //! The type of the jacobian.
      typedef typename Dune::PDELab::BackendMatrixSelector<MB,Domain,Range,JF>::Type Jacobian;

      //! The sparsity pattern container for the jacobian matrix
      typedef typename MB::template Pattern<Jacobian,GFSV,GFSU> Pattern;

      //! The local assembler type
      typedef BatchedTBBLocalAssembler<BatchedTBBGridOperator,LOP,
                                       nonoverlapping_mode> LocalAssembler;

      typedef typename conditional<
        nonoverlapping_mode,
        NonOverlappingBorderDOFExchanger<BatchedTBBGridOperator>,
        OverlappingBorderDOFExchanger<BatchedTBBGridOperator>
        >::type BorderDOFExchanger;

      //! The grid operator traits
      typedef Dune::PDELab::GridOperatorTraits
      <GFSU,GFSV,MB,DF,RF,JF,CU,CV,Assembler,LocalAssembler> Traits;

      template <typename MFT>
      struct MatrixContainer{
        typedef typename Traits::Jacobian Type;
      };

      //! Constructor for non trivial constraints
      BatchedTBBGridOperator(const GFSU & gfsu_, const CU & cu_,
                             const GFSV & gfsv_, const CV & cv_, LOP & lop_,
                             const MB& mb_ = MB())
        : global_assembler(gfsu_,gfsv_,cu_,cv_)
        , dof_exchanger(make_shared<BorderDOFExchanger>(*this))
        , local_assembler(lop_, cu_, cv_,dof_exchanger)
        , backend(mb_)
      {}

      //! Constructor for empty constraints
      BatchedTBBGridOperator(const GFSU & gfsu_, const GFSV & gfsv_,
                             LOP & lop_, const MB& mb_ = MB())
        : global_assembler(gfsu_,gfsv_)
        , dof_exchanger(make_shared<BorderDOFExchanger>(*this))
        , local_assembler(lop_,dof_exchanger)
        , backend(mb_)
      {}

      //! Get the trial grid function space
      const GFSU& trialGridFunctionSpace() const
      {
        return global_assembler.trialGridFunctionSpace();
      }

      //! Get the test grid function space
      const GFSV& testGridFunctionSpace() const
      {
        return global_assembler.testGridFunctionSpace();
      }

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

      Assembler & assembler() { return global_assembler; }

      const Assembler & assembler() const { return global_assembler; }

      LocalAssembler & localAssembler() const { return local_assembler; }


      //! Visitor which is called in the method setupGridOperators for
      //! each tuple element.
      template <typename GridOperatorTuple>
      struct SetupGridOperator {
        SetupGridOperator()
          : index(0), size(Dune::tuple_size<GridOperatorTuple>::value) {}

        template <typename T>
        void visit(T& elem) {
          elem.localAssembler().doPreProcessing = index == 0;
          elem.localAssembler().doPostProcessing = index == size-1;
          ++index;
        }

        int index;
        const int size;
      };

      //! Method to set up a number of grid operators which are used
      //! in a joint assembling. It is assumed that all operators are
      //! specializations of the same template type
      template<typename GridOperatorTuple>
      static void setupGridOperators(GridOperatorTuple tuple)
      {
        Dune::ForEachValue<GridOperatorTuple> forEach(tuple);
        SetupGridOperator<GridOperatorTuple> setup_visitor;
        forEach.apply(setup_visitor);
      }

      //! Interpolate the constrained dofs from given function
      template<typename F, typename X>
      void interpolate (const X& xold, F& f, X& x) const
      {
        // Interpolate f into grid function space and set corresponding coefficients
        Dune::PDELab::interpolate(f,global_assembler.trialGridFunctionSpace(),x);

        // Copy non-constrained dofs from old time step
        Dune::PDELab::copy_nonconstrained_dofs(local_assembler.trialConstraints(),xold,x);
      }

      //! Fill pattern of jacobian matrix
      void fill_pattern(Pattern & p) const {
        typedef typename LocalAssembler::LocalPatternAssemblerEngine PatternEngine;
        global_assembler.assemble
          ([&,p](LocalAssembler &la) mutable -> PatternEngine&
                { return la.localPatternAssemblerEngine(p); },
           local_assembler);
      }

      //! Assemble residual
      void residual(const Domain & x, Range & r) const {
        typedef typename LocalAssembler::LocalResidualAssemblerEngine ResidualEngine;
        global_assembler.assemble
          ([&,r,x](LocalAssembler &la) mutable -> ResidualEngine&
                {
                  ResidualEngine & engine =
                    la.localResidualAssemblerEngine(r,x);
                  engine.setAutocommit(1000);
                  return engine;
                },
           local_assembler);
      }

      //! Assembler jacobian
      void jacobian(const Domain & x, Jacobian & a) const {
        typedef typename LocalAssembler::LocalJacobianAssemblerEngine JacobianEngine;
        global_assembler.assemble
          ([&,a,x](LocalAssembler &la) mutable -> JacobianEngine&
                {
                  JacobianEngine & engine =
                    la.localJacobianAssemblerEngine(a,x);
                  engine.setAutocommit(1000);
                  return engine;
                },
           local_assembler);
      }

      //! Apply jacobian matrix without explicitly assembling it
      void jacobian_apply(const Domain & x, Range & r) const {
        typedef typename LocalAssembler::LocalJacobianApplyAssemblerEngine JacobianApplyEngine;
        global_assembler.assemble
          ([&,r,x](LocalAssembler &la) -> JacobianApplyEngine&
                { return la.localJacobianApplyAssemblerEngine(r,x); },
           local_assembler);
      }

      void make_consistent(Jacobian& a) const {
        dof_exchanger->accumulateBorderEntries(*this,a);
      }

      void update()
      {
        // the DOF exchanger has matrix information, so we need to update it
        dof_exchanger->update(*this);
      }

      //! Get the matrix backend for this grid operator.
      const typename Traits::MatrixBackend& matrixBackend() const
      {
        return backend;
      }

    private:
      Assembler global_assembler;
      shared_ptr<BorderDOFExchanger> dof_exchanger;

      mutable LocalAssembler local_assembler;
      MB backend;

    };

  }
}
#endif //  DUNE_PDELAB_GRIDOPERATOR_TBB_HH
