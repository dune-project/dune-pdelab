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
#include <ranges>
#include <algorithm>
#include <variant>

namespace Dune::PDELab::inline Experimental {


// Forward declaration of a parameter tree
class PropertyTree;

//! Exception thrown when a Property is badly casted
struct BadPropertyCast : Dune::Exception {};

//! Exception thrown when a property tries to acceses an invalidated reference
struct BadPropertyReference : Dune::InvalidStateException {};


namespace Impl {
  using WeakPropertyTree = std::variant<std::weak_ptr<PropertyTree>, std::weak_ptr<const PropertyTree>>;

  //! helper concept to detect if type if a managed property tree
  template<class T>
  concept PropertyTreePtr =
    (std::same_as<T, std::shared_ptr<typename std::remove_reference_t<T>::element_type>&> ||
     std::same_as<T, std::weak_ptr<typename std::remove_reference_t<T>::element_type>&>)  &&
     std::derived_from<typename std::remove_reference_t<T>::element_type, PropertyTree>;
}



/**
 * @brief Property holding a value, a name, a setter, a getter and a documentation.
 * @details The value of the property is dynamically stored in std::any and its contents
 * must be accessed via helper functions like property_cast and unwrap_property_ref.
 * Whenever a the property value is assigned, the setter function (if any) is invoked.
 * Similarly, whenever a property is casted to its underlying contents, the getter function
 * (if any) is invoked.
 *
 * The assignement operator sets the value of the property with a type of the (deduced) template argument.
 * On the other hand, the function call property_cast retreives the content of value
 * with the type of the (deduced) template argument. Both the assignement operator
 * and the property_cast calls must have exactly the same template argument for the same property,
 * otherwise, a BadPropertyCast will be throwed. Most of the cast supported by std::any_cast are
 * supported by property_cast.
 *
 * To ease the access to (maybe) reference types (i.e. std::[shared|weak]_ptr, or std::reference_wrapper)
 * the unwrap_property_ref function unwraps these objects or an into a reference for homogeneous access.
 *
 * Finally, since Property and PropertyTree have an special relation between them,
 * Property may be casted into a PropertyTree via the function as_tree. This cast is
 * allowed if the stored value of the property is PropertyTree or std::[shared|weak]_ptr<[const] T>
 * whenever T is a base class of PropertyTree (i.e. polymorphism may be used via smart pointers).
 *
 * Properties are not thread-safe per se, but they are guaranteed to be
 */
class Property {

public:
  //! The type of the setter functor
  using Setter = std::function<void(std::reference_wrapper<Property>)>;
  //! The type of the getter functor
  using Getter = std::function<void(std::variant<std::reference_wrapper<Property>, std::reference_wrapper<const Property>>)>;

  /**
   * @brief Construct an empty Property object
   *
   * @param name_           The name of the property
   * @param documentation_  The documentation of the Property
   * @param setter_         The setter functor of the Property
   * @param getter_         The getter functor of the Property
   */
  Property(std::string name_ = {}, std::string documentation_ = {}, Setter setter_ = {}, Getter getter_ = {});

  //! Copy constructs a property from another property
  Property(const Property&);

  //! Move constructs a property out of another property
  Property(Property&&);

  //! @brief Destroy the Property object
  ~Property();

  Property& operator=(const Property&);
  Property& operator=(Property&);
  Property& operator=(Property&&);
  Property& operator=(std::nullptr_t);

  template<class T>
  Property& operator=(std::initializer_list<T> list);

  template<class T>
  Property& operator=(T&& t);

  template<class T>
  requires std::derived_from<std::remove_const_t<T>, PropertyTree>
  Property& operator=(std::weak_ptr<T> ptree_ptr);

  template<class T>
  requires std::derived_from<std::remove_const_t<T>, PropertyTree>
  Property& operator=(std::shared_ptr<T> ptree_ptr);

  operator std::any&() = delete;
  operator const std::any&() const = delete;

  //! Checks whether the object contains a value.
  bool has_value() const { return _object.has_value(); };
  //! Returns the typeid of the contained value.
  const std::type_info& type() const { return _object.type(); };

  //! Checks whether the object contains a PropertyTree or `std::[shared|weak]_ptr<[const] T>` whenever T is a base class of PropertyTree.
  bool has_tree() const;

  //! Returns a ParameterTree. If property has no value, a new ParameterTree is created.
  PropertyTree& as_tree();
  //! Returns a ParameterTree.
  const PropertyTree& as_tree() const;

  //! Checks whether the object contains a std::vector<Property>
  bool has_vector() const;
  //! Returns a std::vector<Property>. If property has no value, a new std::vector<Property> is created.
  std::vector<Property>& as_vector() &;
  //! Returns a std::vector<Property>. If property has no value, a new std::vector<Property> is created.
  std::vector<Property>&& as_vector() &&;
  //! Returns a std::vector<Property>
  const std::vector<Property>& as_vector() const &;

  friend std::ostream& operator<<(std::ostream& out, const Property& ptree);
  std::set<PropertyTree const *> report(std::ostream& out, std::string indent) const;

  /**
   * @brief Performs type-safe access to the contained property value.
   * @warning The cast when T is `std::[shared|weak]_ptr<[const] V>&` where
   * V is a base class of PropertyTree is deleted (note that other cv-ref
   * qualifiers are allowed). This case would allow users to directily manipulate
   * the contents of the smart pointer, thus, allowing to unintentionally invalidate
   * the internal property registry.
   *
   * @tparam T    Type to cast the property value
   * @tparam U    Type of the property
   * @param ppt   Target property object
   * @return      Property value with type T
   * @throws      Throws BadPropertyCast if the typeid of the requested T does not match that of the contents of ppt.
   */
  template<class T, class U>
  requires std::derived_from<std::remove_reference_t<U>, Property>
  friend T property_cast(U&& ppt);

  // Deleted property cast:  T is `std::[shared|weak]_ptr<[const] V>&` where V is a base class of PropertyTree.
  template<class T, class U>
  requires (std::derived_from<std::remove_reference_t<U>, Property> && Impl::PropertyTreePtr<T>)
  friend T property_cast(U&& val) = delete;

  /**
   * @brief Performs type-safe access to the contained property value with a default value.
   * @details
   *
   * @tparam T
   * @param val
   * @param default_val
   * @return std::remove_cvref_t<T>&
   */
  template<class T>
  friend std::remove_cvref_t<T>& property_cast(Property& val, T&& default_val);

  // Deleted property cast: T is `std::[shared|weak]_ptr<[const] V>&` where V is a base class of PropertyTree.
  template<Impl::PropertyTreePtr T>
  friend std::remove_cvref_t<T>& property_cast(Property& val, T&& default_val) = delete;

  template<class T = Property>
  std::ranges::view auto vector_view() const &;

  template<class T>
  std::ranges::view auto vector_view() &;

  template<class T>
  std::ranges::view auto vector_view() &&;

  // simple fallback for simple containers (not sfinae friendly)
  template<class Container, class... Args>
  requires (std::ranges::sized_range<Container> &&
           (not std::is_reference_v<Container>) &&
           (not std::ranges::view<Container>))
  Container vector_to(Args&&... args) const;

  BadPropertyCast bad_cast_exception(const std::type_info& info, std::optional<std::string> other_info = std::nullopt) const;

private:
  // TODO: this kind of formatting is not compatible with std::format(...);
  void static register_format(const std::type_info& type, std::function<std::string(const Property&)> f);

  template<class T>
  void static register_format(std::function<std::string(const Property&)> f);

  template<class T>
  void static register_format();

  void register_derived_ptree(std::optional<Impl::WeakPropertyTree> ptree);

  template<class T>
  void maybe_register_derived_ptree(std::weak_ptr<T> ptree_ptr) {
    if constexpr (not std::same_as<PropertyTree,std::remove_const_t<T>>) {
      if constexpr (std::is_const_v<T>)
        register_derived_ptree(std::weak_ptr<const PropertyTree>(ptree_ptr));
      else
        register_derived_ptree(std::weak_ptr<PropertyTree>(ptree_ptr));
    }
  }

public:
  std::string name;
  std::string documentation;
  Setter setter;
  Getter getter;

private:
  std::any _object;
  std::function<void()> _clean_up;
};


template<class T>
Property& Property::operator=(std::initializer_list<T> list) {
  std::vector<Property> vec(list.size());
  auto list_it = list.begin();
  for(std::size_t i = 0; i != list.size(); ++i, ++list_it)
    vec[i] = *list_it;
  return *this = std::move(vec);
}

template<class T>
Property& Property::operator=(T&& t) {
  *this = nullptr;
  _object = std::make_any<T>(std::forward<T>(t));
  if (setter) setter(*this);
  return *this;
}

template<class T>
requires std::derived_from<std::remove_const_t<T>, PropertyTree>
Property& Property::operator=(std::weak_ptr<T> ptree_ptr) {
  *this = nullptr;
  _object = ptree_ptr;
  maybe_register_derived_ptree(ptree_ptr);
  if (setter) setter(*this);
  return *this;
}

template<class T>
requires std::derived_from<std::remove_const_t<T>, PropertyTree>
Property& Property::operator=(std::shared_ptr<T> ptree_ptr) {
  *this = nullptr;
  _object = ptree_ptr;
  maybe_register_derived_ptree(std::weak_ptr(ptree_ptr));
  if (setter) setter(*this);
  return *this;
}

template<class T, class U>
requires std::derived_from<std::remove_reference_t<U>, Property>
T property_cast(U&& val) {
  try {
    if (val.getter) val.getter(std::ref(val));
    return std::any_cast<T>(std::forward<U>(val)._object);
  } catch (std::bad_any_cast& ex) {
    throw val.bad_cast_exception(typeid(T));
  }
}

template<class T>
std::remove_cvref_t<T>& property_cast(Property& val, T&& default_val) {
  if (not val.has_value())
    val = std::forward<T>(default_val);
  return property_cast<std::remove_cvref_t<T>&>(val);
}

template<class T>
requires std::same_as<T, std::add_const_t<std::remove_reference_t<T>>>
T& unwrap_property_ref(const Property& val) {
  using U = std::remove_cvref_t<T>;
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
  else if (val.type() ==  typeid(T))
    return property_cast<T&>(val);
  throw val.bad_cast_exception(typeid(T&));
}

template<class T>
requires std::same_as<T, std::remove_cvref_t<T>>
T&& unwrap_property_ref(Property&& val) {
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
requires std::same_as<T, std::remove_cvref_t<T>>
T& unwrap_property_ref(Property& val) {
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
std::remove_cvref_t<T>& unwrap_property_ref(Property& val, T&& default_val) {
  if (not val.has_value())
    val = std::forward<T>(default_val);
  return unwrap_property_ref<std::remove_cvref_t<T>>(val);
}

template<class T>
std::remove_cvref_t<T>& unwrap_property_ref(Property& val, std::shared_ptr<T> default_val) {
  if (not val.has_value())
    val = std::move(default_val);
  return unwrap_property_ref<std::remove_cvref_t<T>>(val);
}

template<class T>
std::remove_cvref_t<T>& unwrap_property_ref(Property& val, std::weak_ptr<T> default_val) {
  if (not val.has_value())
    val = std::move(default_val);
  return unwrap_property_ref<std::remove_cvref_t<T>>(val);
}

template<class T>
std::remove_cvref_t<T>& unwrap_property_ref(Property& val, std::reference_wrapper<T> default_val) {
  if (not val.has_value())
    val = std::move(default_val);
  return unwrap_property_ref<std::remove_cvref_t<T>>(val);
}



template<class T>
std::ranges::view auto Property::vector_view() const & {
  std::size_t bound;
  std::function<T(std::size_t)> transform;
  if (has_vector()) {
    const auto& vec = as_vector();
    bound = vec.size();
    transform = [&](std::size_t i){ return unwrap_property_ref<const T>(vec[i]); };
  } else {
    bound = 1;
    transform = [&](std::size_t i){ return unwrap_property_ref<const T>(*this); };
  }
  return std::views::iota(std::size_t{0}, bound) | std::views::transform(transform);
}

template<class T>
std::ranges::view auto Property::vector_view() & {
  std::size_t bound;
  std::function<T&(std::size_t)> transform;
  if (has_vector()) {
    auto& vec = as_vector();
    bound = vec.size();
    transform = [&](std::size_t i) -> T& { return unwrap_property_ref<T>(vec[i]); };
  } else {
    bound = 1;
    transform = [&](std::size_t i) -> T& { return unwrap_property_ref<T>(*this); };
  }
  return std::views::iota(std::size_t{0}, bound) | std::views::transform(transform);
}

template<class T>
std::ranges::view auto Property::vector_view() && {
  std::size_t bound;
  std::function<T&&(std::size_t)> transform;
  if (has_vector()) {
    auto&& vec = as_vector();
    bound = vec.size();
    transform = [&](std::size_t i) -> T&& { return unwrap_property_ref<T>(std::move(vec[i])); };
  } else {
    bound = 1;
    transform = [&](std::size_t i) -> T&& { return unwrap_property_ref<T>(std::move(*this)); };
  }
  return std::views::iota(std::size_t{0}, bound) | std::views::transform(transform);
}

template<class Container, class... Args>
requires (std::ranges::sized_range<Container> &&
          (not std::is_reference_v<Container>) &&
          (not std::ranges::view<Container>))
Container Property::vector_to(Args&&... args) const {
  // warning: this is a very naive implementation which happens to work with
  // basic stuff (std::vector, std::array, Dune::FieldVector)
  // but is probably wrong for many kind of ranges!
  Container c(std::forward<Args>(args)...);

  using Value = std::ranges::range_value_t<Container>;
  auto view = vector_view<Value>();
  if constexpr (requires { c.reserve(view.size()); })
    c.reserve(view.size());

  if constexpr( requires { c.push_back(std::declval<Value>()); } )
    for(auto&& entry : view) c.push_back(entry);
  else if constexpr( requires { c.insert(c.end(), std::declval<Value>()); } )
    for(auto&& entry : view) c.push_back(entry);
  else {
    auto cit = c.begin();
    auto vit = view.begin();
    while (cit != c.end() and vit != view.end())
      *(cit++) = *(vit++);
  }

  return c;
}

template<class T>
void Property::register_format(std::function<std::string(const Property&)> f) {
  register_format(typeid(T), f);
}

template<class T>
void Property::register_format() {
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

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_PROPERTY_HH
