// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_TEST_FMT_HH
#define DUNE_PDELAB_TEST_FMT_HH

#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/fvector.hh>

template<typename T, int n>
std::string fmt(const Dune::FieldVector<T, n>& v) {
  std::ostringstream s;
  s << "(";
  if(n > 0)
    s << v[0];
  for(unsigned i = 1; i < n; ++i)
    s << ", " << v[i];
  s << ")";
  return s.str();
}

template<typename U, typename V>
std::string fmt(const std::pair<U, V>& p) {
  std::ostringstream s;
  s << "("<< p.first << ", " << p.second << ")";
  return s.str();
}

template<typename T>
std::string fmt(const std::vector<T> &v) {
  std::ostringstream s;
  s << "(";
  if(v.size() > 0)
    s << v[0];
  for(unsigned i = 1; i < v.size(); ++i)
    s << ", " << v[i];
  s << ")";
  return s.str();
}

#endif // DUNE_PDELAB_TEST_FMT_HH
