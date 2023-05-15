#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/common/property.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/overloadset.hh>

#include <gtest/gtest.h>

#include <list>

TEST(TestProperty, References) {
  Dune::PDELab::Property val;
  int& i = unwrap_property_ref<int>(val = 1);
  EXPECT_EQ(i, 1);
  Dune::PDELab::Property ref;
  int& i2 = unwrap_property_ref<int>(ref = std::ref(i));
  EXPECT_EQ(i2, 1);
  EXPECT_EQ(&i2, &i);
  int&& m2 = unwrap_property_ref<int>(std::move(ref));
  EXPECT_EQ(m2, 1);
  Dune::PDELab::Property cref;
  const int& i3 = unwrap_property_ref<const int>(cref = std::cref(i2));
  EXPECT_EQ(&m2, &i3);
  auto i4 = std::make_shared<int>(2);
  Dune::PDELab::Property ptr;
  int& i5 = unwrap_property_ref<int>(ptr = i4);
  EXPECT_EQ(&*i4, &i5);
  EXPECT_EQ(i5, 2);
  const int& i6 = unwrap_property_ref<const int>(ptr);
  EXPECT_EQ(&*i4, &i6);
  int&& i7 = unwrap_property_ref<int>(std::move(ptr));
  EXPECT_EQ(i7, 2);
}

TEST(TestProperty, Vector) {
  using namespace Dune::PDELab;
  Property ppt;
  std::vector<Property>& vector = ppt.as_vector();
  vector.resize(7);

  vector[0] = 1;
  vector[2] = true;
  vector[3] = 'c';
  vector[4] = "foo1";
  vector[5].documentation = "Heterogeneous and nested vector";
  std::vector<Property>& vector5 = vector[5].as_vector();
  vector5.resize(2);
  vector5[2] = "foo2";
  std::vector<Property>& vector51 = vector5[1].as_vector();
  vector51.resize(1);
  vector51[0] = "foo3";
  vector[6] = {"foo4", "foo5"};

  EXPECT_EQ(property_cast<char const *>(vector[6].as_vector()[1]), "foo5");

  std::vector<Property>& vector50 = vector5[0].as_vector();
  vector50.resize(10);
  for (std::size_t i = 0; i != 10; ++i)
    vector50[i] = i;

  std::vector<std::size_t> vec(vector50.size());
  for (std::size_t i = 0; i != vec.size(); ++i)
    vec[i] = property_cast<std::size_t>(vector50[i]);

  for (std::size_t i = 0; i != vec.size(); ++i)
    EXPECT_EQ(vec[i], i);

  auto view = vector5[0].vector_view<std::size_t>();
  for (auto& entry : view) entry += 2;

  auto v = vector5[0].vector_to<std::vector<std::size_t>>(5);
  // auto v = view | std::ranges::to<std::vector<std::size_t>>(5); // c++23
  EXPECT_EQ(v.size(), 15);
  for (std::size_t i = 0; i != 10; ++i) EXPECT_EQ(v[i+5], i+2);

  auto l = vector5[0].vector_to<std::list<std::size_t>>();
  // auto l = view | std::ranges::to<std::list<std::size_t>>(); // c++23
  EXPECT_EQ(l.size(), 10);
  for (std::size_t count = 2; std::size_t entry : l) EXPECT_EQ(entry, count++);

  auto a = vector5[0].vector_to<std::array<std::size_t,5>>();
  // auto a = view | std::ranges::to<std::array<std::size_t,5>>(); // c++23
  EXPECT_EQ(a.size(), 5);
  for (std::size_t i = 0; i != a.size(); ++i) EXPECT_EQ(a[i], i+2);

  auto f = vector5[0].vector_to<Dune::FieldVector<std::size_t,3>>();
  // auto f = view | std::ranges::to<Dune::FieldVector<std::size_t,10>>(); // c++23
  EXPECT_EQ(f.size(), 3);
  for (std::size_t i = 0; i != f.size(); ++i) EXPECT_EQ(f[i], i+2);

  EXPECT_EQ(property_cast<std::size_t>(ppt.as_vector()[5].as_vector()[0].as_vector()[9]), 11);
  std::cout << ppt << std::endl;
}

TEST(TestProperty, SetterGetter) {
  using Dune::PDELab::Property;
  Property ppt;
  bool const_getter_used = false;
  bool mut_getter_used = false;
  bool setter_used = false;

  ppt.documentation = "Integer with setter, getter and documentation";
  ppt.setter = [&](Property& val){ setter_used = true; };
  ppt.getter = [&](auto variant){
    std::visit(Dune::overload(
      [&](const Property& val) { const_getter_used = true; },
      [&](Property& val)       { mut_getter_used = true; val = 3; }
    ), variant);
  };
  EXPECT_FALSE(setter_used or mut_getter_used or const_getter_used);

  ppt = 2;
  EXPECT_TRUE(setter_used);
  EXPECT_FALSE(mut_getter_used or const_getter_used);

  [[maybe_unused]] int& val = property_cast<int&>(ppt);
  EXPECT_TRUE(setter_used and mut_getter_used);
  EXPECT_FALSE(const_getter_used);
  EXPECT_EQ(val, 3);

  ppt = 2;
  [[maybe_unused]] const int& cval = property_cast<const int&>(std::as_const(ppt));
  EXPECT_TRUE(setter_used and mut_getter_used and const_getter_used);
  EXPECT_EQ(cval, 2);

  std::cout << ppt << std::endl;
}
