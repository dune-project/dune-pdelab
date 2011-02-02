// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TYPETREE_COMPOSITENODE_MACROS_HH
#define DUNE_PDELAB_COMMON_TYPETREE_COMPOSITENODE_MACROS_HH

#if HAVE_VARIADIC_TEMPLATES

#include <dune/pdelab/common/typetree/variadiccompositenode.hh>

#define DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN typename... Children

#define DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN_FOR_SPECIALIZATION typename... Children

#define DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES Children...

#define DUNE_TYPETREE_COMPOSITENODE_BASETYPE Dune::PDELab::TypeTree::VariadicCompositeNode<Children...>

#define DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES children...

#define DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES_THROUGH_FUNCTION(f) f(children)...

#define DUNE_TYPETREE_COMPOSITENODE_CONSTRUCTOR_SIGNATURE Children&... children

#define DUNE_TYPETREE_COMPOSITENODE_STORAGE_CONSTRUCTOR_SIGNATURE shared_ptr<Children>... children

#else

#include <dune/pdelab/common/typetree/compositenode.hh>

#define DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN \
  typename C0,\
  typename C1 = Dune::PDELab::TypeTree::EmptyNode,       \
  typename C2 = Dune::PDELab::TypeTree::EmptyNode,       \
  typename C3 = Dune::PDELab::TypeTree::EmptyNode,       \
  typename C4 = Dune::PDELab::TypeTree::EmptyNode,       \
  typename C5 = Dune::PDELab::TypeTree::EmptyNode,       \
  typename C6 = Dune::PDELab::TypeTree::EmptyNode,       \
  typename C7 = Dune::PDELab::TypeTree::EmptyNode,       \
  typename C8 = Dune::PDELab::TypeTree::EmptyNode,       \
  typename C9 = Dune::PDELab::TypeTree::EmptyNode

#define DUNE_TYPETREE_COMPOSITENODE_TEMPLATE_CHILDREN_FOR_SPECIALIZATION \
  typename C0,\
  typename C1,\
  typename C2,\
  typename C3,\
  typename C4,\
  typename C5,\
  typename C6,\
  typename C7,\
  typename C8,\
  typename C9

#define DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES C0,C1,C2,C3,C4,C5,C6,C7,C8,C9

#define DUNE_TYPETREE_COMPOSITENODE_BASETYPE Dune::PDELab::TypeTree::CompositeNode<C0,C1,C2,C3,C4,C5,C6,C7,C8,C9>

#define DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES c0,c1,c2,c3,c4,c5,c6,c7,c8,c9

#define DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES_THROUGH_FUNCTION(f) f(c0),f(c1),f(c2),f(c3),f(c4),f(c5),f(c6),f(c7),f(c8),f(c9)

#define DUNE_TYPETREE_COMPOSITENODE_CONSTRUCTOR_SIGNATURE \
  C0& c0, \
  typename Dune::PDELab::TypeTree::OptionalChild<C1>::type c1 = Dune::PDELab::TypeTree::OptionalChild<C1>::default_value(), \
  typename Dune::PDELab::TypeTree::OptionalChild<C2>::type c2 = Dune::PDELab::TypeTree::OptionalChild<C2>::default_value(), \
  typename Dune::PDELab::TypeTree::OptionalChild<C3>::type c3 = Dune::PDELab::TypeTree::OptionalChild<C3>::default_value(), \
  typename Dune::PDELab::TypeTree::OptionalChild<C4>::type c4 = Dune::PDELab::TypeTree::OptionalChild<C4>::default_value(), \
  typename Dune::PDELab::TypeTree::OptionalChild<C5>::type c5 = Dune::PDELab::TypeTree::OptionalChild<C5>::default_value(), \
  typename Dune::PDELab::TypeTree::OptionalChild<C6>::type c6 = Dune::PDELab::TypeTree::OptionalChild<C6>::default_value(), \
  typename Dune::PDELab::TypeTree::OptionalChild<C7>::type c7 = Dune::PDELab::TypeTree::OptionalChild<C7>::default_value(), \
  typename Dune::PDELab::TypeTree::OptionalChild<C8>::type c8 = Dune::PDELab::TypeTree::OptionalChild<C8>::default_value(), \
  typename Dune::PDELab::TypeTree::OptionalChild<C9>::type c9 = Dune::PDELab::TypeTree::OptionalChild<C9>::default_value()

#define DUNE_TYPETREE_COMPOSITENODE_STORAGE_CONSTRUCTOR_SIGNATURE \
  shared_ptr<C0> c0, \
  shared_ptr<C1> c1, \
  shared_ptr<C2> c2, \
  shared_ptr<C3> c3, \
  shared_ptr<C4> c4, \
  shared_ptr<C5> c5, \
  shared_ptr<C6> c6, \
  shared_ptr<C7> c7, \
  shared_ptr<C8> c8, \
  shared_ptr<C9> c9

#endif

#endif // DUNE_PDELAB_COMMON_TYPETREE_COMPOSITENODE_MACROS_HH
