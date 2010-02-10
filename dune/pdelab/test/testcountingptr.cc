// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include <iostream>

#include <dune/common/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include"../common/countingptr.hh"

class A : public Dune::PDELab::Countable
{
public:
  void hello () const
  {
	std::cout << "object of class A" << std::endl;
  }
  void hello ()
  {
    std::cout << "object of const class A" << std::endl;
  }
};

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

	// make object
	A a;
	std::cout << "reference count is " << a.get_reference_counter() << std::endl;
	if (a.get_reference_counter()!=0)
	  return 1;

    // mutable stuff

	// make first pointer
	Dune::PDELab::CP<A> cp1(&a);
	std::cout << "reference count is " << a.get_reference_counter() << std::endl;
	if (a.get_reference_counter()!=1)
	  return 2;

	// make a second pointer
	Dune::PDELab::CP<A> cp2;
	cp2 = &a;
	std::cout << "reference count is " << a.get_reference_counter() << std::endl;
	if (a.get_reference_counter()!=2)
	  return 3;

	// reset first pointer
	cp1 = 0;
	std::cout << "reference count is " << a.get_reference_counter() << std::endl;
	if (a.get_reference_counter()!=1)
	  return 4;

    // assign mutable pointer from mutable pointer
    cp1 = cp2;
    std::cout << "reference count is " << a.get_reference_counter()
              << std::endl;
    if (a.get_reference_counter()!=2)
      return 4;

    // const stuff

    // make pointer to const from mutable object
    Dune::PDELab::CP<const A> ccp1(&a);
    std::cout << "reference count is " << a.get_reference_counter()
              << std::endl;
    if (a.get_reference_counter()!=3)
      return 5;

    const A ca;

    // make pointer to const from const object
    Dune::PDELab::CP<const A> ccp2(&ca);
    std::cout << "reference count is " << ca.get_reference_counter()
              << std::endl;
    if (ca.get_reference_counter()!=1)
      return 6;

    // assign a const pointer from a const pointer
    ccp1 = ccp2;
    std::cout << "reference count is " << a.get_reference_counter()
              << std::endl;
    if (a.get_reference_counter()!=2)
      return 7;
    std::cout << "reference count is " << ca.get_reference_counter()
              << std::endl;
    if (ca.get_reference_counter()!=2)
      return 8;

    // assign a const pointer from a mutable pointer
    ccp1 = cp2;
    std::cout << "reference count is " << a.get_reference_counter()
              << std::endl;
    if (a.get_reference_counter()!=3)
      return 9;
    std::cout << "reference count is " << ca.get_reference_counter()
              << std::endl;
    if (ca.get_reference_counter()!=1)
      return 10;

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
