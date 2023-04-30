#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/common/property_tree.hh>

#include <dune/common/exceptions.hh>

#include <gtest/gtest.h>

TEST(TestProperty, References) {
  Dune::PDELab::Property val = 1;
  int& i = unwrap_property_ref<int>(val);
  EXPECT_EQ(i, 1);
  Dune::PDELab::Property ref = std::ref(i);
  int& i2 = unwrap_property_ref<int>(ref);
  EXPECT_EQ(i2, 1);
  EXPECT_EQ(&i2, &i);
  int&& m2 = unwrap_property_ref<int>(std::move(ref));
  EXPECT_EQ(m2, 1);
  Dune::PDELab::Property cref = std::cref(i2);
  const int& i3 = unwrap_property_ref<const int>(cref);
  EXPECT_EQ(&m2, &i3);
  auto i4 = std::make_shared<int>(2);
  Dune::PDELab::Property ptr = i4;
  int& i5 = unwrap_property_ref<int>(ptr);
  EXPECT_EQ(&*i4, &i5);
  EXPECT_EQ(i5, 2);
  const int& i6 = unwrap_property_ref<const int>(ptr);
  EXPECT_EQ(&*i4, &i6);
  int&& i7 = unwrap_property_ref<int>(std::move(ptr));
  EXPECT_EQ(i7, 2);
}


TEST(TestPropertyTree, SubPropertyTree) {
  Dune::PDELab::PropertyTree foo;
  foo["value1"] = 1;
  foo["value2"] = 2;

  EXPECT_EQ(foo.template get<int>("value1"), 1);
  EXPECT_EQ(foo.template get<int>("value2"), 2);

  foo["value2"] = std::monostate{}; // override type
  foo["value2"] = 3; // override type again
  EXPECT_EQ(foo.template get<int>("value2"), 3);

  foo.erase("value2");
  EXPECT_FALSE(foo.hasKey("value2"));

  foo["bla.bla.bla"] = "2";
  EXPECT_EQ(&foo.sub("bla").sub("bla").get("bla"), &foo["bla.bla.bla"]);
  EXPECT_TRUE(foo.hasSub("bla"));
  EXPECT_TRUE(foo.hasSub("bla.bla"));
  EXPECT_FALSE(foo.hasSub("bla.bla.bla"));
  std::cout << foo << std::endl;
}

TEST(TestPropertyTree, References) {
  using namespace Dune::PDELab;
  std::shared_ptr<PropertyTree> foo_ptr = std::make_shared<PropertyTree>();
  PropertyTree& foo = *foo_ptr;
  foo["value"] = 1;
  int i = 2;
  foo["value_ptr"] = &i;
  EXPECT_EQ(foo.template get<int*>("value_ptr"), &i);

  auto shared = std::make_shared<int>(4);
  foo["shared"] = shared;
  EXPECT_EQ(foo.get("shared", 2), 4);
  EXPECT_EQ(foo.template get<int>("shared"), 4);

  foo["weak"] = std::weak_ptr(shared);
  EXPECT_EQ(foo.get("weak", 2), 4);
  EXPECT_EQ(foo.template get<int>("weak"), 4);

  foo.erase("shared");
  shared = nullptr;
  EXPECT_THROW(foo.get("weak", 2), Dune::PDELab::BadPropertyReference);

  foo["shared"] = std::make_shared<const int>(4);
  EXPECT_THROW(foo.get("shared", 2), Dune::PDELab::BadPropertyCast);

  EXPECT_THROW(foo.sub("value"), Dune::RangeError);
  EXPECT_THROW(foo.sub("value_ptr"), Dune::RangeError);

  foo["sub_ref"] = std::weak_ptr(foo_ptr);
  EXPECT_EQ(&foo.sub("sub_ref.sub_ref.sub_ref.sub_ref"), &foo);
  EXPECT_EQ(foo.template get<int>("sub_ref.sub_ref.sub_ref.sub_ref.value"), 1);
  foo["copy"] = foo;
  foo["copy.value"] = 2;
  EXPECT_NE(&foo["value"], &foo["copy.value"]);
  EXPECT_EQ(&foo["value"], &foo["copy.sub_ref.value"]);
  std::cout << foo << std::endl;
}

TEST(TestPropertyTree, BaseClass) {
  using namespace Dune::PDELab;
  PropertyTree boo;
  struct Foo : public PropertyTree {};
  auto foo = std::make_shared<Foo>();
  foo->get("value") = 1;
  foo->get("self") = std::weak_ptr(foo);
  boo["foo1"] = std::weak_ptr(foo);
  boo["foo2"] = foo;
  foo->get("value") = 2;
  EXPECT_EQ(boo.sub("foo1.self.self").template get<int>("value"), 2);
  EXPECT_EQ(boo.sub("foo2.self.self").template get<int>("value"), 2);
  std::cout << boo << std::endl;
}

TEST(TestPropertyTree, Cycle) {
  using namespace Dune::PDELab;
  auto foo = std::make_shared<PropertyTree>();
  {
    auto boo = std::make_shared<PropertyTree>();
    foo->get("value") = 1;
    foo->get("boo") = std::weak_ptr(boo);
    boo->get("foo") = std::weak_ptr(foo);
    foo->get("value") = 2;
    EXPECT_EQ(boo->sub("foo.boo.foo").template get<int>("value"), 2);
    std::cout << *boo << std::endl;
  }
  EXPECT_THROW(foo->sub("boo.foo"), Dune::InvalidStateException);
}

TEST(TestPropertyTree, SetterGetter) {
  using namespace Dune::PDELab;
  PropertyTree dict;
  bool getter_used = false;
  bool setter_used = false;

  dict["value"].documentation = "Integer with setter, getter and documentation";
  dict["value"].setter = [&](const Property& val){ setter_used = true; };
  dict["value"].getter = [&](const Property& val){ getter_used = true; };
  EXPECT_FALSE(setter_used);
  EXPECT_FALSE(getter_used);

  dict["value"] = 2;
  EXPECT_TRUE(setter_used);
  EXPECT_FALSE(getter_used);

  [[maybe_unused]] int two = dict.get("value", 2);
  EXPECT_TRUE(setter_used);
  EXPECT_TRUE(getter_used);

  std::cout << dict << std::endl;
}


TEST(TestPropertyTree, Array) {
  using namespace Dune::PDELab;
  PropertyTree dict;

  dict["value"] = 1;
  dict["array"][0] = 1;
  dict["array"][2] = true;
  dict["array"][3] = 'c';
  dict["array"][4] = "foo1";
  dict["array"][5][2] = "foo2";
  dict["array"][5][2].documentation = "Heterogeneous array";
  dict["array"][5][1][0] = "foo3";

  for (std::size_t i = 0; i != 10; ++i)
    dict["array"][5][0][i] = i;

  std::vector<std::size_t> vec(dict["array"][5][0].size());
  for (std::size_t i = 0; i != vec.size(); ++i)
    vec[i] = property_cast<std::size_t>(dict["array"][5][0][i]);

  for (std::size_t i = 0; i != vec.size(); ++i)
      EXPECT_EQ(vec[i], i);

  std::cout << dict << std::endl;
}
