// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#include "config.h"

#include <iostream>
#include <algorithm>
#include <numeric>

#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>
#include <dune/istl/bvector.hh>

#include <dune/pdelab/backend/istl.hh>

int main(int argc, char** argv)
{

  {
    typedef Dune::BlockVector<Dune::FieldVector<int,1> > V;
    typedef Dune::PDELab::ISTL::vector_iterator<V> Iterator;

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
    typedef Dune::PDELab::ISTL::vector_iterator<V> Iterator;

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
    typedef Dune::PDELab::ISTL::vector_iterator<V> Iterator;

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
    typedef Dune::PDELab::ISTL::vector_iterator<V> Iterator;

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
    typedef Dune::BlockVector<Dune::DynamicVector<int> > V;
    typedef Dune::PDELab::ISTL::vector_iterator<V> Iterator;

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
    typedef Dune::BlockVector<Dune::BlockVector<Dune::DynamicVector<int> > > V;
    typedef Dune::PDELab::ISTL::vector_iterator<V> Iterator;

    V v(20);
    int i = 1;
    for (V::iterator it = v.begin(),
           end = v.end();
         it != end;
         ++it, ++i)
      {
        it->resize(i);
        int j = 1;
        for (V::block_type::iterator it2 = it->begin(),
           end2 = it->end();
         it2 != end2;
         ++it2, ++j)
          it2->resize(j);
      }

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
    typedef Dune::BlockVector<Dune::BlockVector<Dune::DynamicVector<int> > > V;
    typedef Dune::PDELab::ISTL::vector_iterator<V> Iterator;

    V v(20);
    int i = 19;
    for (V::iterator it = v.begin(),
           end = v.end();
         it != end;
         ++it, --i)
      {
        it->resize(i+1);
        int j = i;
        for (V::block_type::iterator it2 = it->begin(),
           end2 = it->end();
         it2 != end2;
         ++it2, --j)
          it2->resize(j);
      }

    v = 1;

    std::cout << v.one_norm() << std::endl;

    for (Iterator it = Iterator(v,false),
           end = Iterator(v,true);
         it != end;
         ++it)
      (*it) *= 3;

    std::cout << v.one_norm() << std::endl;
  }

  // const test

  {
    typedef Dune::BlockVector<Dune::FieldVector<int,1> > V;
    typedef Dune::PDELab::ISTL::vector_iterator<const V> Iterator;

    V v(20);
    v = 1;

    std::cout << v.one_norm() << std::endl;

    int sum = std::accumulate(Iterator(v,false),Iterator(v,true),0);

    std::cout << sum << std::endl;
  }

  // non-const -> const assignment test

  {
    typedef Dune::BlockVector<Dune::FieldVector<int,1> > V;
    typedef Dune::PDELab::ISTL::vector_iterator<V> Iterator;
    typedef Dune::PDELab::ISTL::vector_iterator<const V> ConstIterator;

    V v(20);
    v = 1;

    std::cout << v.one_norm() << std::endl;

    int sum = 0;

    for (Iterator it = Iterator(v,false),
           end = Iterator(v,true);
         it != end;
         ++it)
      {
        ConstIterator cit(end);
        cit = it;
        if (cit != it)
          {
            std::cerr << "Equality postcondition failed!" << std::endl;
            return 1;
          }
        sum += *cit;
      }

    std::cout << sum << std::endl;
  }


  return 0;
}
