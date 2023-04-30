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
#include <typeinfo>
#include <functional>
#include <sstream>
#include <set>

namespace Dune::PDELab::inline Experimental {

class PropertyTree;


struct BadPropertyCast : Dune::Exception {};
struct BadPropertyReference : Dune::InvalidStateException {};

class Property {

public:

  Property() = default;
  Property(std::string_view);
  Property(const Property&);
  Property(Property&&);
  ~Property();

  Property& operator=(const Property&);
  Property& operator=(Property&&);

  operator std::any&() = delete;
  operator const std::any&() const = delete;

  std::size_t size() const { return as_array().size(); }
  const std::type_info& type() const { return _object.type(); };
  std::string_view name() const { return _name; }

  bool has_tree() const;
  PropertyTree& as_tree();
  const PropertyTree& as_tree() const;

  bool has_array() const;
  std::vector<Property>& as_array();
  const std::vector<Property>& as_array() const;

  Property& operator[](std::size_t i);
  const Property& operator[](std::size_t i) const;

  friend std::ostream& operator<<(std::ostream& out, const Property& ptree);
  std::set<PropertyTree const *> report(std::ostream& out, std::string indent) const;

  template<class T>
  Property(T&& arg) {
    *this = std::forward<T>(arg);
  }

  Property& operator=(std::nullptr_t) {
    if (auto clean_up = std::move(_clean_up)) clean_up();
    _object.reset();
    return *this;
  }

  template<class T>
  Property& operator=(T&& t) {
    *this = nullptr;
    _object = std::make_any<T>(std::forward<T>(t));
    if (setter) setter(*this);
    return *this;
  }

  template<class T>
  requires (std::is_base_of_v<PropertyTree, T> and not std::same_as<PropertyTree,T>)
  Property& operator=(std::weak_ptr<T> ptree_ptr) {
    *this = nullptr;
    _object = ptree_ptr;
    register_ptree(ptree_ptr);
    if (setter) setter(*this);
    return *this;
  }

  template<class T>
  requires (std::is_base_of_v<PropertyTree, T> and not std::same_as<PropertyTree,T>)
  Property& operator=(std::shared_ptr<T> ptree_ptr) {
    std::weak_ptr<PropertyTree> ptree_weak = ptree_ptr;
    *this = nullptr;
    _object = std::move(ptree_ptr);
    register_ptree(ptree_weak);
    if (setter) setter(*this);
    return *this;
  }

  template<class T, class U>
  requires std::derived_from<std::remove_reference_t<U>, Property>
  friend T property_cast(U&& val) {
    if (val.getter) val.getter(val);
    try {
      return std::any_cast<T>(std::forward<U>(val)._object);
    } catch (std::bad_any_cast& ex) {
      throw val.bad_cast_exception(className<T>());
    }
  }

  template<class T>
  friend std::remove_cvref_t<T>& property_cast(Property& val, T&& default_val) {
    if (not val._object.has_value())
      val = std::forward<T>(default_val);
    return property_cast<std::remove_cvref_t<T>&>(val);
  }

  template<class T>
  requires (not std::is_reference_v<T>)
  friend const T& unwrap_property_ref(const Property& val) {
    using U = std::remove_const_t<T>;
    if (val.type() == typeid(std::shared_ptr<U>))
      return *property_cast<const std::shared_ptr<U>&>(val);
    if (val.type() == typeid(std::shared_ptr<const U>))
      return *property_cast<const std::shared_ptr<const U>&>(val);
    else if (val.type() == typeid(std::weak_ptr<U>)) {
      if (auto observe_ptr = property_cast<const std::weak_ptr<U>&>(val).lock())
        return *observe_ptr;
      else
        DUNE_THROW(BadPropertyReference, "Weak pointer to \'" + className<U>() + "\' not exist anymore");
    } else if (val.type() == typeid(std::weak_ptr<const U>)) {
      if (auto observe_ptr = property_cast<const std::weak_ptr<const U>&>(val).lock())
        return *observe_ptr;
      else
        DUNE_THROW(BadPropertyReference, "Weak pointer to \'" + className<U>() + "\' not exist anymore");
    } else if (val.type() ==  typeid(std::reference_wrapper<U>))
      return property_cast<std::reference_wrapper<U>>(val).get();
    else if (val.type() ==  typeid(std::reference_wrapper<const U>))
      return property_cast<std::reference_wrapper<const U>>(val).get();
    else
      return property_cast<const T&>(val);
  }

  template<class T>
  requires (not std::is_reference_v<T>)
  friend T&& unwrap_property_ref(Property&& val) {
    if (val.type() == typeid(std::shared_ptr<T>))
      return std::move(*property_cast<const std::shared_ptr<T>&>(val));
    else if (val.type() == typeid(std::weak_ptr<T>)) {
      if (auto observe_ptr = property_cast<const std::weak_ptr<T>&>(val).lock())
        return std::move(*observe_ptr);
      else
        DUNE_THROW(BadPropertyReference, "Weak pointer to \'" + className<T>() + "\' not exist anymore");
    } else if (val.type() ==  typeid(std::reference_wrapper<T>))
      return std::move(property_cast<std::reference_wrapper<T>>(val).get());
    else
      return property_cast<T&&>(std::move(val));
  }

  template<class T>
  requires (not std::is_reference_v<T> && not std::is_const_v<T>)
  friend T& unwrap_property_ref(Property& val) {
    if (val.type() == typeid(std::shared_ptr<T>))
      return *property_cast<const std::shared_ptr<T>&>(val);
    else if (val.type() == typeid(std::weak_ptr<T>)) {
      if (auto observe_ptr = property_cast<const std::weak_ptr<T>&>(val).lock())
        return *observe_ptr;
      else
        DUNE_THROW(BadPropertyReference, "Weak pointer to \'" + className<T>() + "\' not exist anymore");
    } else if (val.type() ==  typeid(std::reference_wrapper<T>))
      return property_cast<std::reference_wrapper<T>>(val).get();
    else
      return property_cast<T&>(val);
  }


  template<class T>
  friend std::remove_cvref_t<T>& unwrap_property_ref(Property& val, T&& default_val) {
    if (val._object.has_value())
      return unwrap_property_ref<std::remove_cvref_t<T>>(val);
    else
      return property_cast<std::remove_cvref_t<T>&>(val = std::forward<T>(default_val));
  }

private:

  BadPropertyCast bad_cast_exception(const std::string& type) const;

  void register_ptree(std::weak_ptr<PropertyTree> ptree);
  void static register_format(const std::type_info& type, std::function<std::string(const Property&)> f);

  template<class T>
  void static register_format(std::function<std::string(const Property&)> f) {
    register_format(typeid(T), f);
  }

  template<class T>
  void static register_format() {
    register_format(typeid(T), [](const Property& val) {
      std::stringstream ss;
      ss << property_cast<const T&>(val);
      return ss.str();
    });

    register_format(typeid(T*), [](const Property& val){
      std::stringstream ss;
      auto ptr = property_cast<T*>(val);
      ss << "[" << ptr << "]";
      if (ptr) ss << " " << *ptr;
      return ss.str();
    });

    if constexpr (not std::is_const_v<T>)
      register_format<const T>();
  }

public:
  std::string documentation;
  std::function<void(const Property&)> getter;
  std::function<void(Property&)> setter;
private:
  std::any _object;
  std::string _name;
  std::function<void()> _clean_up;
};

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_PROPERTY_HH
