#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/common/property_tree.hh>

#include <dune/common/exceptions.hh>

#include <gtest/gtest.h>

#include <list>

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
  boo["foo1"] = std::weak_ptr<const Foo>(foo);
  boo["foo2"] = foo;
  foo->get("value") = 2;
  EXPECT_THROW(boo.sub("foo1"), Dune::RangeError);
  EXPECT_NO_THROW(std::as_const(boo).sub("foo1"));
  EXPECT_THROW(std::as_const(boo).sub("foo3"), Dune::RangeError);
  EXPECT_EQ(std::as_const(boo).sub("foo1.self.self").template get<int>("value"), 2);
  EXPECT_EQ(boo.sub("foo2.self.self").template get<int>("value"), 2);
  std::cout << boo << std::endl;
}

TEST(TestPropertyTree, Cycle) {
  using namespace Dune::PDELab;
  struct Foo : PropertyTree {};
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
