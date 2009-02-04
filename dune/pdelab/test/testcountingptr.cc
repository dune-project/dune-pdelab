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
};

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

	// make object
	A a;
	std::cout << "reference count is " << a.get_reference_counter() << std::endl;
	if (a.get_reference_counter()!=0)
	  return 1;

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
