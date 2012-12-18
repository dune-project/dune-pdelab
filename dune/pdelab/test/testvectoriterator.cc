// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#include "config.h"

#include <iostream>
#include <algorithm>

#include <dune/pdelab/backend/istlvectorbackend.hh>


int main(int argc, char** argv)
{

  {
    typedef Dune::BlockVector<Dune::FieldVector<int,1> > V;
    typedef Dune::PDELab::istl::vector_iterator<V> Iterator;

    V v(20);
    v = 1;

    std::cout << v.one_norm() << std::endl;

    for (Iterator it = Iterator(v,false),
           end = Iterator(v,true);
         it != end;
         ++it)
      (*it) *= 3;

    std::cout << v.one_norm() << std::endl;
  }

  {
    typedef Dune::BlockVector<Dune::FieldVector<int,3> > V;
    typedef Dune::PDELab::istl::vector_iterator<V> Iterator;

    V v(20);

    v = 1;

    std::cout << v.one_norm() << std::endl;

    for (Iterator it = Iterator(v,false),
           end = Iterator(v,true);
         it != end;
         ++it)
      (*it) *= 3;

    std::cout << v.one_norm() << std::endl;
  }

  {
    typedef Dune::BlockVector<Dune::BlockVector<Dune::FieldVector<int,1> > > V;
    typedef Dune::PDELab::istl::vector_iterator<V> Iterator;

    V v(20);
    int i = 1;
    for (V::iterator it = v.begin(),
           end = v.end();
         it != end;
         ++it, ++i)
      it->resize(i);

    v = 1;

    std::cout << v.one_norm() << std::endl;

    for (Iterator it = Iterator(v,false),
           end = Iterator(v,true);
         it != end;
         ++it)
      (*it) *= 3;

    std::cout << v.one_norm() << std::endl;
  }

  {
    typedef Dune::BlockVector<Dune::BlockVector<Dune::FieldVector<int,3> > > V;
    typedef Dune::PDELab::istl::vector_iterator<V> Iterator;

    V v(20);
    int i = 1;
    for (V::iterator it = v.begin(),
           end = v.end();
         it != end;
         ++it, ++i)
      it->resize(i);

    v = 1;

    std::cout << v.one_norm() << std::endl;

    for (Iterator it = Iterator(v,false),
           end = Iterator(v,true);
         it != end;
         ++it)
      (*it) *= 3;

    std::cout << v.one_norm() << std::endl;
  }

  return 0;
}
