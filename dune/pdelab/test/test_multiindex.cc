#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/common/multiindex.hh>

#include <dune/pdelab/concepts/multiindex.hh>

#include <dune/common/test/testsuite.hh>

#include <unordered_map>

int main()
{
  using namespace Dune::Assembler;
  using namespace Dune::Indices;

  Dune::TestSuite suite("Check MultiIndex()");

  constexpr auto mi = multiIndex(_1,3,_2,5);
  static_assert(Concept::MultiIndex<std::remove_cvref_t<decltype(mi)>>);

  suite.check(mi[_0] == 1);
  suite.check(mi[3] == 5);

  suite.check(back(mi) == 5);
  static_assert(back(push_back(mi, _3)) == 3);
  suite.check(back(push_back(mi, 3)) == 3);
  static_assert(back(pop_back(mi)) == 2);
  suite.check(back(pop_back(mi)) == 2);
  suite.check(back(pop_back(pop_back(mi))) == 3);

  static_assert(front(mi) == 1);
  static_assert(front(push_front(mi,_0)) == 0);
  suite.check(front(push_front(mi,0)) == 0);
  suite.check(front(pop_front(mi)) == 3);
  static_assert(front(pop_front(pop_front(mi))) == 2);
  suite.check(front(pop_front(pop_front(mi))) == 2);

  using MI = std::decay_t<decltype(mi)>;
  std::unordered_map<MI, int, Dune::Assembler::Hash<MI>> map;

  map[mi] = 2;
  map[multiIndex(_1,4,_2,6)] = 3;

  suite.check(map[mi] == 2);
  suite.check(map[multiIndex(_1,4,_2,6)] == 3);

}
