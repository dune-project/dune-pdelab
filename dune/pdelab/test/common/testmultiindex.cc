#include <cassert>
#include <iostream>
#include <sstream>
#include <dune/typetree/hybridmultiindex.hh>
#include <dune/pdelab/common/multiindex.hh>

#include <dune/common/test/testsuite.hh>

int main()
{
  Dune::TestSuite test("MultiIndexTestSuite");

  using MI = Dune::PDELab::MultiIndex<int, 4>;

  // Default construction
  {
    MI mi;
    test.check(mi.empty(), "Default constructed MultiIndex should be empty");
    test.check(mi.size() == 0, "Default constructed MultiIndex size should be 0");
  }

  // Push back and access
  {
    MI mi;
    mi.push_back(1);
    mi.push_back(2);
    mi.push_back(3);
    test.check(mi.size() == 3, "Push back 3 elements");
    test.check(mi[0] == 1, "First element after push_back");
    test.check(mi[1] == 2, "Second element after push_back");
    test.check(mi[2] == 3, "Third element after push_back");
  }

  // Push front
  {
    MI mi;
    mi.push_back(1);
    mi.push_back(2);
    mi.push_back(3);
    mi.push_front(0);
    test.check(mi.size() == 4, "Push front increases size");
    test.check(mi[0] == 0, "First element after push_front");
    test.check(mi[1] == 1, "Second element after push_front");
  }

  // Pop back
  {
    MI mi;
    mi.push_back(0);
    mi.push_back(1);
    mi.push_back(2);
    mi.push_back(3);
    MI mi2 = pop_back(mi);
    test.check(mi.size() == 4, "Original size after pop_back");
    test.check(mi2.size() == 3, "Size after pop_back");
    test.check(mi2[0] == 0, "First element after pop_back");
    test.check(mi2[2] == 2, "Last element after pop_back");
  }

  // Pop front
  {
    MI mi;
    mi.push_back(0);
    mi.push_back(1);
    mi.push_back(2);
    MI mi3 = pop_front(mi);
    test.check(mi.size() == 3, "Original size after pop_front");
    test.check(mi3.size() == 2, "Size after pop_front");
    test.check(mi3[0] == 1, "First element after pop_front");
  }

  // Accumulate back
  {
    MI mi;
    mi.push_back(2);
    mi.push_back(5);
    test.check(accumulate_back(mi, 5).back() == 10, "Accumulate back");
    test.check(mi.back() == 5, "Accumulate back");
  }

  // Accumulate front
  {
    MI mi;
    mi.push_back(10);
    mi.push_back(1);
    test.check(accumulate_front(mi, 1).front() == 11, "Accumulate front");;
    test.check(mi.front() == 10, "Accumulate front");
  }

  // Join
  {
    MI mi3;
    mi3.push_back(11);
    mi3.push_back(7);
    MI mi4;
    mi4.push_back(8);
    mi4.push_back(9);
    MI joined = join(mi3, mi4);
    test.check(joined.size() == 4, "Join size");
    test.check(joined[0] == 11, "Join first element");
    test.check(joined[2] == 8, "Join third element");
    test.check(joined[3] == 9, "Join fourth element");
  }

  // Reverse
  {
    MI joined;
    joined.push_back(11);
    joined.push_back(7);
    joined.push_back(8);
    joined.push_back(9);
    MI rev = reverse(joined);
    test.check(rev[0] == 9, "Reverse first element");
    test.check(rev[3] == 11, "Reverse last element");
  }

  // View
  {
    MI joined;
    joined.push_back(11);
    joined.push_back(7);
    joined.push_back(8);
    joined.push_back(9);
    auto v = joined.view();
    test.check(v.size() == joined.size(), "View size");
    test.check(v[1] == joined[1], "View element");
  }

  // Copy constructor from ReservedVector
  {
    Dune::ReservedVector<int, 4> rv;
    rv.push_back(4);
    rv.push_back(5);
    MI mi5(rv);
    test.check(mi5.size() == 2, "Copy from ReservedVector size");
    test.check(mi5[0] == 4, "Copy from ReservedVector first element");
  }

  // Equality and inequality
  {
    MI mi5;
    mi5.push_back(4);
    mi5.push_back(5);
    MI mi6;
    mi6.push_back(4);
    mi6.push_back(5);
    test.check(mi5 == mi6, "Equality operator");
    mi6.push_back(6);
    test.check(mi5 != mi6, "Inequality operator");
  }

  // Output operator
  {
    MI mi5;
    mi5.push_back(4);
    mi5.push_back(5);
    std::ostringstream oss;
    oss << mi5;
    test.check(oss.str().find("4") != std::string::npos, "Output operator contains 4");
  }

  // Hash
  {
    MI mi5;
    mi5.push_back(4);
    mi5.push_back(5);
    MI mi6;
    mi6.push_back(4);
    mi6.push_back(5);
    mi6.push_back(6);
    std::size_t h1 = hash_value(mi5);
    std::size_t h2 = hash_value(mi6);
    test.check(h1 != h2, "Hash values differ");
  }

  // Construction from HybridTreePath
  {
    Dune::TypeTree::HybridTreePath<std::size_t, std::size_t, std::size_t> tp(7, 8, 9);
    Dune::PDELab::MultiIndex<int, 3> mi7(tp);
    test.check(mi7.size() == 3, "HybridTreePath size");
    test.check(mi7[0] == 7, "HybridTreePath first element");
    test.check(mi7[2] == 9, "HybridTreePath last element");
  }

  return test.exit();
}
