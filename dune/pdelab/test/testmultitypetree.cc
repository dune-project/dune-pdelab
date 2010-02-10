// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include <iostream>

#include <dune/common/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include"../common/countingptr.hh"
#include"../common/multitypetree.hh"
#include"../common/cpstoragepolicy.hh"

class A : public Dune::PDELab::Countable
{
public:
  void hello () const
  {
	std::cout << "object of class A" << std::endl;
  }
  int number () const
  {
	return 17;
  }
};

class B : public Dune::PDELab::Countable, public Dune::PDELab::LeafNode
{
public:
  void hello () const
  {
	std::cout << "object of class B" << std::endl;
  }
  int number () const
  {
	return 33;
  }
};


int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

	// make objects
	A a;
	B b;

	// power node
	Dune::PDELab::PowerNode<A,2,Dune::PDELab::CountingPointerStoragePolicy> tree(a,a);
 	std::cout << "made power node with 2 children" << std::endl;
  	if (tree.getChild<0>().number()!=17)
	  return 1;
   	if (tree.getChild<1>().number()!=17)
	  return 1;

    // power node with const members
    Dune::PDELab::PowerNode<const A,2,
      Dune::PDELab::CountingPointerStoragePolicy> ctree(a,a);
    std::cout << "made power node with 2 const children" << std::endl;
    if (ctree.getChild<0>().number()!=17)
      return 1;
    if (ctree.getChild<1>().number()!=17)
      return 1;

	// composite node 
	Dune::PDELab::CompositeNode<Dune::PDELab::CountingPointerStoragePolicy,
	  A,B,A,B,A,B,A,B,A> composite(a,b,a,b,a,b,a,b,a);
   	std::cout << "made composite node with 9 children" << std::endl;
 	if (composite.getChild<0>().number()!=17)
	  return 3;
   	if (composite.getChild<1>().number()!=33)
	  return 4;

	Dune::PDELab::CompositeNode<Dune::PDELab::CountingPointerStoragePolicy,
	  A,B> composite2(a,b);
  	std::cout << "made composite node with 2 children" << std::endl;
	Dune::PDELab::CompositeNode<Dune::PDELab::CountingPointerStoragePolicy,
	  A,B,A> composite3(a,b,a);
  	std::cout << "made composite node with 3 children" << std::endl;
	Dune::PDELab::CompositeNode<Dune::PDELab::CountingPointerStoragePolicy,
	  A,B,A,B> composite4(a,b,a,b);
  	std::cout << "made composite node with 4 children" << std::endl;
	Dune::PDELab::CompositeNode<Dune::PDELab::CountingPointerStoragePolicy,
	  A,B,A,B,A> composite5(a,b,a,b,a);
  	std::cout << "made composite node with 5 children" << std::endl;
	Dune::PDELab::CompositeNode<Dune::PDELab::CountingPointerStoragePolicy,
	  A,B,A,B,A,B> composite6(a,b,a,b,a,b);
  	std::cout << "made composite node with 6 children" << std::endl;
	Dune::PDELab::CompositeNode<Dune::PDELab::CountingPointerStoragePolicy,
	  A,B,A,B,A,B,A> composite7(a,b,a,b,a,b,a);
  	std::cout << "made composite node with 7 children" << std::endl;
	Dune::PDELab::CompositeNode<Dune::PDELab::CountingPointerStoragePolicy,
	  A,B,A,B,A,B,A,B> composite8(a,b,a,b,a,b,a,b);
  	std::cout << "made composite node with 8 children" << std::endl;
	Dune::PDELab::CompositeNode<Dune::PDELab::CountingPointerStoragePolicy,
	  A,B,A,B,A,B,A,B,A> composite9(a,b,a,b,a,b,a,b,a);
  	std::cout << "made composite node with 9 children" << std::endl;
    if (a.get_reference_counter()!=33)
	  return 5;

    // composite node (const members)
    Dune::PDELab::CompositeNode<Dune::PDELab::CountingPointerStoragePolicy,
      const A,const B,const A,const B,const A,const B,const A,const B,const A>
      const_composite(a,b,a,b,a,b,a,b,a);
    std::cout << "made composite node with 9 const children" << std::endl;
    if (const_composite.getChild<0>().number()!=17)
      return 3;
    if (const_composite.getChild<1>().number()!=33)
      return 4;

    Dune::PDELab::CompositeNode<Dune::PDELab::CountingPointerStoragePolicy,
      const A,const B> const_composite2(a,b);
    std::cout << "made composite node with 2 const children" << std::endl;
    Dune::PDELab::CompositeNode<Dune::PDELab::CountingPointerStoragePolicy,
      const A,const B,const A> const_composite3(a,b,a);
    std::cout << "made composite node with 3 const children" << std::endl;
    Dune::PDELab::CompositeNode<Dune::PDELab::CountingPointerStoragePolicy,
      const A,const B,const A,const B> const_composite4(a,b,a,b);
    std::cout << "made composite node with 4 const children" << std::endl;
    Dune::PDELab::CompositeNode<Dune::PDELab::CountingPointerStoragePolicy,
      const A,const B,const A,const B,const A> const_composite5(a,b,a,b,a);
    std::cout << "made composite node with 5 const children" << std::endl;
    Dune::PDELab::CompositeNode<Dune::PDELab::CountingPointerStoragePolicy,
      const A,const B,const A,const B,const A,const B>
      const_composite6(a,b,a,b,a,b);
    std::cout << "made composite node with 6 const children" << std::endl;
    Dune::PDELab::CompositeNode<Dune::PDELab::CountingPointerStoragePolicy,
      const A,const B,const A,const B,const A,const B,const A>
      const_composite7(a,b,a,b,a,b,a);
    std::cout << "made composite node with 7 const children" << std::endl;
    Dune::PDELab::CompositeNode<Dune::PDELab::CountingPointerStoragePolicy,
      const A,const B,const A,const B,const A,const B,const A,const B>
      const_composite8(a,b,a,b,a,b,a,b);
    std::cout << "made composite node with 8 const children" << std::endl;
    Dune::PDELab::CompositeNode<Dune::PDELab::CountingPointerStoragePolicy,
      const A,const B,const A,const B,const A,const B,const A,const B,const A>
      const_composite9(a,b,a,b,a,b,a,b,a);
    std::cout << "made composite node with 9 const children" << std::endl;
    if (a.get_reference_counter()!=62)
      return 6;

	// test passed
	return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
	return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
	return 1;
  }
} 
