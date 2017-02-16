#ifndef DUNE_PDELAB_GRIDOPERATOR_GRIDOPERATOR_HH
#define DUNE_PDELAB_GRIDOPERATOR_GRIDOPERATOR_HH

// TODO: Remove deprecated nonoverlapping parameter. This will break
// boilerplate code and needs to be done with care.

#include <tuple>

#include <dune/common/hybridutilities.hh>

#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperator/common/borderdofexchanger.hh>
#include <dune/pdelab/gridoperator/common/gridoperatorutilities.hh>
#include <dune/pdelab/gridoperator/default/assembler.hh>
#include <dune/pdelab/gridoperator/default/localassembler.hh>

namespace Dune{
  namespace PDELab{

    namespace impl {

      template<int i>
      struct warn_on_deprecated_nonoverlapping_mode_parameter
      {

        warn_on_deprecated_nonoverlapping_mode_parameter() DUNE_DEPRECATED_MSG(
"the nonoverlapping_mode parameter on the GridOperator has been deprecated and will be removed after PDELab 2.4."
"The correct mode is now automatically deduced from the EntitySet of the function space.")
        {}

      };

      template<>
      struct warn_on_deprecated_nonoverlapping_mode_parameter<-1>
      {};

    }

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
             int nonoverlapping_mode = -1>
    class GridOperator
      : public impl::warn_on_deprecated_nonoverlapping_mode_parameter<nonoverlapping_mode>
    {
    public:

      static_assert(nonoverlapping_mode == -1 ||
                    nonoverlapping_mode == 0 ||
                    nonoverlapping_mode == 1,
                    "invalid value for nonoverlapping_mode! This parameter is also deprecated in PDELab 2.4, so please remove it from your typedefs!");

      //! The global assembler type
      typedef DefaultAssembler<GFSU,GFSV,CU,CV,static_cast<bool>(nonoverlapping_mode)> Assembler;

      //! The type of the domain (solution).
      using Domain = Dune::PDELab::Backend::Vector<GFSU,DF>;
      //! The type of the range (residual).
      using Range = Dune::PDELab::Backend::Vector<GFSV,RF>;
      //! The type of the jacobian.
      using Jacobian = Dune::PDELab::Backend::Matrix<MB,Domain,Range,JF>;

      //! The sparsity pattern container for the jacobian matrix
      typedef typename MB::template Pattern<Jacobian,GFSV,GFSU> Pattern;

      //! The local assembler type
      typedef DefaultLocalAssembler<
        GridOperator,
        LOP,
        GFSU::Traits::EntitySet::Partitions::partitionIterator() == InteriorBorder_Partition
        > LocalAssembler;

      // Fix this as soon as the default Partitions are constexpr
      typedef typename std::conditional<
        GFSU::Traits::EntitySet::Partitions::partitionIterator() == InteriorBorder_Partition,
        NonOverlappingBorderDOFExchanger<GridOperator>,
        OverlappingBorderDOFExchanger<GridOperator>
        >::type BorderDOFExchanger;

      //! The grid operator traits
      typedef Dune::PDELab::GridOperatorTraits
      <GFSU,GFSV,MB,DF,RF,JF,CU,CV,Assembler,LocalAssembler> Traits;

      template <typename MFT>
      struct MatrixContainer{
        typedef typename Traits::Jacobian Type;
      };

      //! Constructor for non trivial constraints
      GridOperator(const GFSU & gfsu_, const CU & cu_, const GFSV & gfsv_, const CV & cv_, LOP & lop_, const MB& mb_ = MB())
        : global_assembler(gfsu_,gfsv_,cu_,cv_)
        , dof_exchanger(std::make_shared<BorderDOFExchanger>(*this))
        , local_assembler(lop_, cu_, cv_,dof_exchanger)
        , backend(mb_)
      {}

      //! Constructor for empty constraints
      GridOperator(const GFSU & gfsu_, const GFSV & gfsv_, LOP & lop_, const MB& mb_ = MB())
        : global_assembler(gfsu_,gfsv_)
        , dof_exchanger(std::make_shared<BorderDOFExchanger>(*this))
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
          : index(0), size(std::tuple_size<GridOperatorTuple>::value) {}

        template <typename T>
        void visit(T& elem) {
          elem.localAssembler().preProcessing(index == 0);
          elem.localAssembler().postProcessing(index == size-1);
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
        SetupGridOperator<GridOperatorTuple> setup_visitor;
        Hybrid::forEach(tuple,
                        [&](auto &el) { setup_visitor.visit(el); });
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
        PatternEngine & pattern_engine = local_assembler.localPatternAssemblerEngine(p);
        global_assembler.assemble(pattern_engine);
      }

      //! Assemble residual
      void residual(const Domain & x, Range & r) const {
        typedef typename LocalAssembler::LocalResidualAssemblerEngine ResidualEngine;
        ResidualEngine & residual_engine = local_assembler.localResidualAssemblerEngine(r,x);
        global_assembler.assemble(residual_engine);
      }

      //! Assembler jacobian
      void jacobian(const Domain & x, Jacobian & a) const {
        typedef typename LocalAssembler::LocalJacobianAssemblerEngine JacobianEngine;
        JacobianEngine & jacobian_engine = local_assembler.localJacobianAssemblerEngine(a,x);
        global_assembler.assemble(jacobian_engine);
      }

      //! Apply jacobian matrix without explicitly assembling it
      void jacobian_apply(const Domain & z, Range & r) const {
        typedef typename LocalAssembler::LocalJacobianApplyAssemblerEngine JacobianApplyEngine;
        JacobianApplyEngine & jacobian_apply_engine = local_assembler.localJacobianApplyAssemblerEngine(r,z);
        global_assembler.assemble(jacobian_apply_engine);
      }

      //! Apply jacobian matrix without explicitly assembling it
      void nonlinear_jacobian_apply(const Domain & x, const Domain & z, Range & r) const {
        global_assembler.assemble(local_assembler.localNonlinearJacobianApplyAssemblerEngine(r,x,z));
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
#endif // DUNE_PDELAB_GRIDOPERATOR_GRIDOPERATOR_HH
