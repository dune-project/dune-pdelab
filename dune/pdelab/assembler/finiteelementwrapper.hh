// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_FINITEELEMENTWRAPPER_HH
#define DUNE_PDELAB_ASSEMBLER_FINITEELEMENTWRAPPER_HH

#include <cassert>
#include <type_traits>
#include <optional>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/pdelab/assembler/utility.hh>

namespace Dune {
  namespace PDELab {

    struct OnlyMovable
    {
      OnlyMovable() = default;
      OnlyMovable(const OnlyMovable&) = delete;
      OnlyMovable(OnlyMovable&&) = default;
      OnlyMovable& operator=(const OnlyMovable&) = delete;
      OnlyMovable& operator=(OnlyMovable&&) = delete;
    };

    template<typename FE, typename Cell>
    class FiniteElementWrapper;


    template<typename Basis_,typename Cell>
    class BasisWrapper
      : public OnlyMovable
    {

      using Switch = BasisInterfaceSwitch<Basis_>;

      template<typename, typename>
      friend class FiniteElementWrapper;

    public:

      using Basis = Basis_;
      using Native = Basis;
      using DomainField = typename Switch::DomainField;
      static const int dimDomainLocal = Switch::dimDomainLocal;
      using DomainLocal = typename Switch::DomainLocal;
      using RangeField = typename Switch::RangeField;
      static const int dimRange = Switch::dimRange;
      using Range = typename Switch::Range;
      using size_type = std::size_t;

      const Basis& native() const
      {
        return *_basis;
      }

      size_type size() const
      {
        return _basis->size();
      }

      size_type order() const
      {
        return _basis->order();
      }

      template<typename In, typename Out>
      void evaluateFunction(const In& in, Out& out) const
      {
        _basis->evaluateFunction(in,out);
      }

      template<typename QP>
      const std::vector<Range>& operator()(const QP& qp)
      {
        _values.resize(size());
        _basis->evaluateFunction(Cell::quadratureCoordinate(qp),_values);
        return _values;
      }

    private:

      void setBasis(const Native& basis)
      {
        _basis = &basis;
      }

      const Basis* _basis;
      std::vector<Range> _values;

    };


    struct set_finite_elements
      : public TypeTree::TreePairVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename LFS, typename FE, typename TreePath>
      void leaf(const LFS& lfs, FE& fe, TreePath treePath) const
      {
        fe.setFiniteElement(lfs.finiteElement());
      }

    };

    template<typename FE, typename Cell>
    class FiniteElementWrapper
      : public TypeTree::LeafNode
      , public OnlyMovable
    {

      using Switch = FiniteElementInterfaceSwitch<FE>;
      friend struct set_finite_elements;

    public:

      using Basis         = BasisWrapper<typename Switch::Basis,Cell>;
      using Interpolation = typename Switch::Interpolation;
      using Coefficients  = typename Switch::Coefficients;

      using FiniteElement = FE;
      using Native = FiniteElement;

      const FiniteElement& native() const
      {
        return *_fe;
      }

      Basis& basis()
      {
        return _basis;
      }

      const Basis& basis() const
      {
        return _basis;
      }

      const Interpolation& interpolation() const
      {
        return Switch::interpolation(*_fe);
      }

      const Coefficients& coefficients() const
      {
        return Switch::coefficients(*_fe);
      }

    private:

      void setFiniteElement(const Native& fe)
      {
        _fe = &fe;
        _basis.setBasis(Switch::basis(fe));
      }

      const FiniteElement* _fe;
      Basis _basis;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ASSEMBLER_FINITEELEMENTWRAPPER_HH
