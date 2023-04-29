#ifndef DUNE_PDELAB_COMMON_PROPERTY_HH
#define DUNE_PDELAB_COMMON_PROPERTY_HH

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>

#include <unordered_map>
#include <any>
#include <memory>
#include <string>
#include <string_view>
#include <vector>
#include <typeindex>
#include <functional>
#include <sstream>
#include <set>

namespace Dune::PDELab::inline Experimental {

class PropertyTree;

class Property {

public:
  ~Property() {
    *this = nullptr;
  }

  Property& operator=(std::nullptr_t) {
    if (_clean_up) {_clean_up(); _clean_up = nullptr;};
    object.reset();
    return *this;
  }

  template<class T>
  Property& operator=(T&& t) {
    *this = nullptr;
    object = std::forward<T>(t);
    if (setter) setter(object);
    return *this;
  }

  template<class T>
  requires (std::is_base_of_v<PropertyTree, T> and not std::same_as<PropertyTree,T>)
  Property& operator=(std::weak_ptr<T> ptree_ptr) {
    auto observe_ptr = ptree_ptr.lock();
    *this = nullptr;
    object = std::move(ptree_ptr);
    if (observe_ptr)
      register_ptree(*observe_ptr, true);
    else
      DUNE_THROW(InvalidStateException, "Weak poiner does not point anywhere");
    if (setter) setter(object);
    return *this;
  }

  template<class T>
  requires (std::is_base_of_v<PropertyTree, T> and not std::same_as<PropertyTree,T>)
  Property& operator=(std::shared_ptr<T> ptree_ptr) {
    PropertyTree& ptree = *ptree_ptr;
    *this = nullptr;
    object = std::move(ptree_ptr);
    register_ptree(ptree, true);
    if (setter) setter(object);
    return *this;
  }

  template<class T>
  requires (std::is_base_of_v<PropertyTree, std::remove_cv_t<T>> and not std::same_as<PropertyTree,std::remove_cv_t<T>>)
  Property& operator=(T&& ptree_ptr) {
    *this = nullptr;
    object = std::forward<T>(ptree_ptr);
    register_ptree(std::any_cast<PropertyTree&>(object), false);
    if (setter) setter(object);
    return *this;
  }

  template <class T>
  friend T any_cast(Property& val) {
    if (val.getter) val.getter(val.object);
    try {
      return std::any_cast<T>(val.object);
    } catch (std::bad_any_cast& ex) {
      val.report_bad_cast(className<T>());
      throw ex;
    }
  }

  template <class T>
  friend T any_cast(const Property& val) {
    if (val.getter) val.getter(val.object);
    try {
      return std::any_cast<T>(val.object);
    } catch (std::bad_any_cast& ex) {
      val.report_bad_cast(className<T>());
      throw ex;
    }
  }

  template <class T>
  friend T any_cast(Property&& val) {
    if (val.getter) val.getter(val.object);
    try {
      return std::any_cast<T>(std::move(val.object));
    } catch (std::bad_any_cast& ex) {
      val.report_bad_cast(className<T>());
      throw ex;
    }
  }

  template <class T, class U>
  friend T any_cast(Property& val, U&& default_val) {
    if (not val.object.has_value())
      val = std::forward<U>(default_val);
    return any_cast<T>(val);
  }

  Property& operator[](std::size_t i);
  const Property& operator[](std::size_t i) const;

  std::size_t size() const { return as_array().size(); }

  bool has_property_tree() const;
  PropertyTree& as_property_tree();
  const PropertyTree& as_property_tree() const;

  bool has_array() const;
  std::vector<Property>& as_array();
  const std::vector<Property>& as_array() const;

  friend std::ostream& operator<<(std::ostream& out, const Property& ptree);

  std::string documentation;
  std::function<void(const std::any&)> getter;
  std::function<void(std::any&)> setter;

  std::set<PropertyTree const *> report(std::ostream& out, std::string, std::string indent) const;

private:

  void report_bad_cast(const std::string& type) const;

  void register_ptree(PropertyTree& ptree, bool reference);
  void static register_format(const std::type_index& type, std::function<std::string(const Property&)> f);

  template<class T>
  void static register_format(std::function<std::string(const Property&)> f) {
    register_format(typeid(T), f);
  }

  template<class T>
  void static register_format() {
    register_format(typeid(T), [](const Property& val) {
      std::stringstream ss;
      ss << any_cast<const T&>(val);
      return ss.str();
    });

    register_format(typeid(T*), [](const Property& val){
      std::stringstream ss;
      auto ptr = any_cast<T*>(val);
      ss << "[" << ptr << "]";
      if (ptr) ss << " " << *ptr;
      return ss.str();
    });

    if constexpr (not std::is_const_v<T>)
      register_format<const T>();
  }

  std::any object;
  std::function<void()> _clean_up;
};

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_PROPERTY_HH
