#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/common/dictionary.hh>

#include <dune/common/exceptions.hh>

#include <gtest/gtest.h>


TEST(TestDictionary, SubDictionary) {
  using namespace Dune::PDELab;
  Dictionary foo;
  foo["value1"] = 1;
  foo["value2"] = 2;

  EXPECT_EQ(foo.template get<int>("value1"), 1);
  EXPECT_EQ(foo.template get<int>("value2"), 2);

  foo["value2"] = std::monostate{}; // override type
  foo["value2"] = 3;
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

TEST(TestDictionary, References) {
  using namespace Dune::PDELab;
  std::shared_ptr<Dictionary> foo_ptr = std::make_shared<Dictionary>();
  Dictionary& foo = *foo_ptr;
  foo["value"] = 1;
  int i = 2;
  foo["value_ptr"] = &i;
  EXPECT_EQ(foo.template get<int*>("value_ptr"), &i);

  EXPECT_THROW(foo.sub("value"), Dune::RangeError);
  EXPECT_THROW(foo.sub("value_ptr"), Dune::RangeError);

  foo["sub_ref"] = Dictionary::make_reference(foo_ptr);
  EXPECT_EQ(&foo.sub("sub_ref.sub_ref.sub_ref.sub_ref"), &foo);
  EXPECT_EQ(foo.template get<int>("sub_ref.sub_ref.sub_ref.sub_ref.value"), 1);
  foo["copy"] = foo;
  foo["copy.value"] = 2;
  EXPECT_NE(&foo["value"], &foo["copy.value"]);
  EXPECT_EQ(&foo["value"], &foo["copy.sub_ref.value"]);
  std::cout << foo << std::endl;
}

TEST(TestDictionary, BaseClass) {
  using namespace Dune::PDELab;
  Dictionary boo;
  struct Foo : public Dictionary {};
  auto foo = std::make_shared<Foo>();
  foo->get("value") = 1;
  foo->get("foo1") = Dictionary::make_reference(foo);
  foo->get("foo2") = Dictionary::make_shared(foo);
  boo["foo1"] = Dictionary::make_reference(foo);
  boo["foo2"] = Dictionary::make_shared(foo);
  foo->get("value") = 2;
  EXPECT_EQ(boo.sub("foo1").template get<int>("value"), 2);
  EXPECT_EQ(boo.sub("foo2").template get<int>("value"), 2);
  std::cout << boo << std::endl;
}



TEST(TestDictionary, Cyclic) {
  using namespace Dune::PDELab;
  auto foo = std::make_shared<Dictionary>();
  {
    auto boo = std::make_shared<Dictionary>();
    foo->get("value") = 1;
    foo->get("boo") = Dictionary::make_reference(boo);
    boo->get("foo") = Dictionary::make_reference(foo);
    foo->get("value") = 2;
    EXPECT_EQ(boo->sub("foo.boo.foo").template get<int>("value"), 2);
    std::cout << *boo << std::endl;
  }
  EXPECT_THROW(foo->sub("boo.foo"), Dune::InvalidStateException);
}
