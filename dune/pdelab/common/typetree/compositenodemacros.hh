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

#define DUNE_TYPETREE_COMPOSITENODE_STORAGE_CONSTRUCTOR_SIGNATURE Dune::shared_ptr<Children>... children

#define DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES_MEMBER(member) children.member...

#define DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES_NESTED_TYPE(nested_type) typename Children::nested_type...

#define DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES_NESTED_STATIC_MEMBER(member) Children::member...

#define DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES_THROUGH_META_FUNCTION(meta_function) typename meta_function::template apply<Children>::type...

template<typename... T> struct extract_first_child;

template<typename T0, typename... T>
struct extract_first_child<T0,T...>
{
  typedef T0 type;
};

#define DUNE_TYPETREE_COMPOSITENODE_FIRST_CHILD extract_first_child<Children...>::type

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
  Dune::shared_ptr<C0> c0,                                        \
  Dune::shared_ptr<C1> c1, \
  Dune::shared_ptr<C2> c2, \
  Dune::shared_ptr<C3> c3, \
  Dune::shared_ptr<C4> c4, \
  Dune::shared_ptr<C5> c5, \
  Dune::shared_ptr<C6> c6, \
  Dune::shared_ptr<C7> c7, \
  Dune::shared_ptr<C8> c8, \
  Dune::shared_ptr<C9> c9

#define DUNE_TYPETREE_COMPOSITENODE_CHILDVARIABLES_MEMBER(member) \
  c0.member, \
  c1.member, \
  c2.member, \
  c3.member, \
  c4.member, \
  c5.member, \
  c6.member, \
  c7.member, \
  c8.member, \
  c9.member

#define DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES_NESTED_TYPE(nested_type) \
  typename C0::nested_type, \
  typename C1::nested_type, \
  typename C2::nested_type, \
  typename C3::nested_type, \
  typename C4::nested_type, \
  typename C5::nested_type, \
  typename C6::nested_type, \
  typename C7::nested_type, \
  typename C8::nested_type, \
  typename C9::nested_type

#define DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES_NESTED_STATIC_MEMBER(member) \
  C0::member, \
  C1::member, \
  C2::member, \
  C3::member, \
  C4::member, \
  C5::member, \
  C6::member, \
  C7::member, \
  C8::member, \
  C9::member

#define DUNE_TYPETREE_COMPOSITENODE_CHILDTYPES_THROUGH_META_FUNCTION(meta_function) \
  typename meta_function::template apply<C0>::type, \
  typename meta_function::template apply<C1>::type, \
  typename meta_function::template apply<C2>::type, \
  typename meta_function::template apply<C3>::type, \
  typename meta_function::template apply<C4>::type, \
  typename meta_function::template apply<C5>::type, \
  typename meta_function::template apply<C6>::type, \
  typename meta_function::template apply<C7>::type, \
  typename meta_function::template apply<C8>::type, \
  typename meta_function::template apply<C9>::type

#define DUNE_TYPETREE_COMPOSITENODE_FIRST_CHILD C0

#endif

#endif // DUNE_PDELAB_COMMON_TYPETREE_COMPOSITENODE_MACROS_HH
