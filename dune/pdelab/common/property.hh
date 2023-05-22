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

  //! helper variant to hold constant and mutable property tree weak pointers
  using WeakPropertyTree = std::variant<std::weak_ptr<PropertyTree>, std::weak_ptr<const PropertyTree>>;

  //! helper concept to detect if type is a managed property tree
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
 * Properties are not thread-safe per se, but they are guaranteed to be re-entrant
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

  //! Destroy the Property object
  ~Property();

  //! Copy property out of another property
  Property& operator=(const Property&);

  //! Copy property out of another property
  Property& operator=(Property&);

  //! Move a property out of another property
  Property& operator=(Property&&);

  //! Resets property value to an empty value
  Property& operator=(std::nullptr_t);

  //! Assigns a vector of properties (see as_vector())
  template<class T>
  Property& operator=(std::initializer_list<T> list);

  //! Assigns an arbitrary object of type `std::decay_t<T>` as the property value.
  template<class T>
  Property& operator=(T&& t);

  //! Assigns a ParameterTree smart pointer as the property value.
  template<class T>
  requires std::derived_from<std::remove_const_t<T>, PropertyTree>
  Property& operator=(std::weak_ptr<T> ptree_ptr);

  //! Assigns a ParameterTree smart pointer as the property value.
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
  bool is_tree() const;

  //! Returns a ParameterTree. If property has no value, a new ParameterTree is created.
  PropertyTree& as_tree();
  //! Returns a ParameterTree.
  const PropertyTree& as_tree() const;

  //! Checks whether the object contains a std::vector<Property>
  bool is_vector() const;
  //! Returns a std::vector<Property>. If property has no value, a new std::vector<Property> is created.
  std::vector<Property>& as_vector() &;
  //! Returns a std::vector<Property>. If property has no value, a new std::vector<Property> is created.
  std::vector<Property>&& as_vector() &&;
  //! Returns a std::vector<Property>
  const std::vector<Property>& as_vector() const &;

  //! Report property contents into an output stream
  friend std::ostream& operator<<(std::ostream& out, const Property& ptree);

  //! Report property contents into an output stream
  std::set<PropertyTree const *> report(std::ostream& out, std::string indent) const;

  // Allow property_cast to see underyling object
  template<class T, class U>
  requires std::derived_from<std::remove_reference_t<U>, Property>
  friend T property_cast(U&& ppt);

  // Delete property cast:  T is `std::[shared|weak]_ptr<[const] V>&` where V is a base class of PropertyTree.
  template<class T, class U>
  requires (std::derived_from<std::remove_reference_t<U>, Property> && Impl::PropertyTreePtr<T>)
  friend T property_cast(U&& val) = delete;

  // Allow property_cast to see underyling object
  template<class T>
  friend std::remove_cvref_t<T>& property_cast(Property& val, T&& default_val);

  // Delete property cast: T is `std::[shared|weak]_ptr<[const] V>&` where V is a base class of PropertyTree.
  template<Impl::PropertyTreePtr T>
  friend std::remove_cvref_t<T>& property_cast(Property& val, T&& default_val) = delete;

  //! View property contents as a non-mutable range of properties
  template<class T = Property>
  std::ranges::view auto vector_view() const &;

  //! View property contents as a mutable range of properties
  template<class T>
  std::ranges::view auto vector_view() &;

  //! View property contents as a movable range of properties
  template<class T>
  std::ranges::view auto vector_view() &&;

  /**
   * @brief Constructs a range object with the contents of the property.
   * @warning This is a simple fallback for simple containers (not sfinae friendly)
   *
   * @tparam Container  Container to store the properties
   * @tparam Args       Arguments to construct the container
   */
  template<class Container, class... Args>
  requires (std::ranges::sized_range<Container> &&
           (not std::is_reference_v<Container>) &&
           (not std::ranges::view<Container>))
  Container vector_to(Args&&... args) const;

  /**
   * @brief Creates a bad cast exception containing a message with diagnostics of the problem
   *
   * @param info          run-time type information of the attempted type to cast
   * @param other_info    other additional information to include in the diagnistics
   * @return BadPropertyCast  The exception
   */
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

  //! Name of the property
  std::string name;

  //! Documentation explaining how to use the property and (if necessary) the allowed types
  std::string documentation;

  //! Functor called on every value assignment of the property value
  Setter setter;

  //! Functor called on every value cast of the property value
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

/**
 * @brief Performs type-safe access to the contained property value
 * @warning When T is `std::[shared|weak]_ptr<[const] V>&` where
 * V is a base class of PropertyTree is deleted (note that other cv-ref
 * qualifiers are allowed). This case would allow users to directily manipulate
 * the contents of the smart pointer, thus, allowing to unintentionally invalidate
 * the internal property registry.
 *
 * @tparam T    Type to cast the property value
 * @tparam U    A qualified type of Property
 * @param ppt   Target property object
 * @return      Property value with type T
 * @throws      Throws BadPropertyCast if the typeid of the requested T does not match that of the contents of ppt.
 */
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

/**
 * @brief Performs type-safe access/assign to the contained property value with a default value
 * @details If the property value has not been set, assigns value to default_val.
 * @warning When T is `std::[shared|weak]_ptr<[const] V>&` where
 * V is a base class of PropertyTree is deleted (note that other cv-ref
 * qualifiers are allowed). This case would allow users to directily manipulate
 * the contents of the smart pointer, thus, allowing to unintentionally invalidate
 * the internal property registry.
 *
 * @tparam T            Type of the default value. Access will be attempted with `std::decay_t<T>`.
 * @tparam U            A mutable reference of a property
 * @param ppt           Target property object
 * @param default_val   A default value to assign
 * @throws              Throws BadPropertyCast if the typeid of the requested T does not match that of the contents of ppt.
 * @return              Property value with type std::remove_cvref_t<T>&
 */
template<class T>
std::remove_cvref_t<T>& property_cast(Property& val, T&& default_val) {
  if (not val.has_value())
    val = std::forward<T>(default_val);
  return property_cast<std::remove_cvref_t<T>&>(val);
}

/**
 * @brief Performs type-safe access to a non-mutable reference of the contained property value.
 * @details This will try to access the contents of the property with different
 * reference-like types and unwrap it if necessary to access the reference of the
 * value itself.
 *
 * Reference-like types are:
 *  - std::[shared|weak]_ptr<[const] T>
 *  - std::reference_wrapper<[const] T>
 *
 * If contained type is not a reference-like object, a normal cast will be attempted.
 * @tparam T    The type of the object to unwrap
 * @param ppt   Property to access
 * @return      Reference of type T& to the property value of ppt
 */
template<class T>
requires std::same_as<T, std::add_const_t<std::remove_reference_t<T>>>
T& unwrap_property_ref(const Property& ppt) {
  using U = std::remove_cvref_t<T>;
  if (ppt.type() == typeid(std::shared_ptr<U>))
    return *property_cast<const std::shared_ptr<U>&>(ppt);
  if (ppt.type() == typeid(std::shared_ptr<const U>))
    return *property_cast<const std::shared_ptr<const U>&>(ppt);
  else if (ppt.type() == typeid(std::weak_ptr<U>)) {
    if (auto observe_ptr = property_cast<const std::weak_ptr<U>&>(ppt).lock())
      return *observe_ptr;
    else
      DUNE_THROW(BadPropertyReference, "Weak pointer to \'" + className<U>() + "\' not exist anymore");
  } else if (ppt.type() == typeid(std::weak_ptr<const U>)) {
    if (auto observe_ptr = property_cast<const std::weak_ptr<const U>&>(ppt).lock())
      return *observe_ptr;
    else
      DUNE_THROW(BadPropertyReference, "Weak pointer to \'" + className<U>() + "\' not exist anymore");
  } else if (ppt.type() ==  typeid(std::reference_wrapper<U>))
    return property_cast<std::reference_wrapper<U>>(ppt).get();
  else if (ppt.type() ==  typeid(std::reference_wrapper<const U>))
    return property_cast<std::reference_wrapper<const U>>(ppt).get();
  else if (ppt.type() ==  typeid(T))
    return property_cast<T&>(ppt);
  throw ppt.bad_cast_exception(typeid(T&));
}

/**
 * @brief Performs type-safe access to a movable reference of the contained property value.
 * @details This will try to access the contents of the property with different
 * reference-like types and unwrap it if necessary to access the reference of the
 * value itself.
 *
 * Reference-like types are:
 *  - std::[shared|weak]_ptr<T>
 *  - std::reference_wrapper<T>
 *
 * If contained type is not a reference-like object, a normal cast will be attempted.
 * @tparam T    The type of the object to unwrap
 * @param ppt   Property to access
 * @return      Reference of type T& to the property value of ppt
 */
template<class T>
requires std::same_as<T, std::remove_cvref_t<T>>
T&& unwrap_property_ref(Property&& ppt) {
  if (ppt.type() == typeid(std::shared_ptr<T>))
    return std::move(*property_cast<const std::shared_ptr<T>&>(ppt));
  else if (ppt.type() == typeid(std::weak_ptr<T>)) {
    if (auto observe_ptr = property_cast<const std::weak_ptr<T>&>(ppt).lock())
      return std::move(*observe_ptr);
    else
      DUNE_THROW(BadPropertyReference, "Weak pointer to \'" + className<T>() + "\' not exist anymore");
  } else if (ppt.type() ==  typeid(std::reference_wrapper<T>))
    return std::move(property_cast<std::reference_wrapper<T>>(ppt).get());
  else
    return property_cast<T&&>(std::move(ppt));
}

/**
 * @brief Performs type-safe access to a non-mutable reference of the contained property value.
 * @details This will try to access the contents of the property with different
 * reference-like types and unwrap it if necessary to access the reference of the
 * value itself.
 *
 * Reference-like types are:
 *  - std::[shared|weak]_ptr<T>
 *  - std::reference_wrapper<T>
 *
 * If contained type is not a reference-like object, a normal cast will be attempted.
 * @tparam T    The type of the object to unwrap
 * @param ppt   Property to access
 * @return      Reference of type T& to the property value of ppt
 */
template<class T>
requires std::same_as<T, std::remove_cvref_t<T>>
T& unwrap_property_ref(Property& ppt) {
  if (ppt.type() == typeid(std::shared_ptr<T>))
    return *property_cast<const std::shared_ptr<T>&>(ppt);
  else if (ppt.type() == typeid(std::weak_ptr<T>)) {
    if (auto observe_ptr = property_cast<const std::weak_ptr<T>&>(ppt).lock())
      return *observe_ptr;
    else
      DUNE_THROW(BadPropertyReference, "Weak pointer to \'" + className<T>() + "\' not exist anymore");
  } else if (ppt.type() ==  typeid(std::reference_wrapper<T>))
    return property_cast<std::reference_wrapper<T>>(ppt).get();
  else
    return property_cast<T&>(ppt);
}


/**
 * @brief Performs type-safe access/assign to the contained property value with a default value
 * @details Similar to property_cast with default argument with the semantics of
 * reference-like access of unwrap_property_ref without default argument.
 * @tparam T            Type of the default value. Access will be attempted with `std::decay_t<T>`.
 * @tparam U            A mutable reference of a property
 * @param ppt           Target property object
 * @param default_val   A default value to assign
 * @throws              Throws BadPropertyCast if the typeid of the requested T does not match that of the contents of ppt.
 * @return              Property value with type std::remove_cvref_t<T>&
 */
template<class T>
std::remove_cvref_t<T>& unwrap_property_ref(Property& ppt, T&& default_val) {
  if (not ppt.has_value())
    ppt = std::forward<T>(default_val);
  return unwrap_property_ref<std::remove_cvref_t<T>>(ppt);
}

/**
 * @brief Performs type-safe access/assign to the contained property value with a default value
 * @details Similar to property_cast with default argument with the semantics of
 * reference-like access of unwrap_property_ref without default argument.
 * @tparam T            Type of the default value. Access will be attempted with `std::decay_t<T>`.
 * @tparam U            A mutable reference of a property
 * @param ppt           Target property object
 * @param default_val   A default value to assign
 * @throws              Throws BadPropertyCast if the typeid of the requested T does not match that of the contents of ppt.
 * @return              Property value with type std::remove_cvref_t<T>&
 */
template<class T>
std::remove_cvref_t<T>& unwrap_property_ref(Property& ppt, std::shared_ptr<T> default_val) {
  if (not ppt.has_value())
    ppt = std::move(default_val);
  return unwrap_property_ref<std::remove_cvref_t<T>>(ppt);
}

/**
 * @brief Performs type-safe access/assign to the contained property value with a default value
 * @details Similar to property_cast with default argument with the semantics of
 * reference-like access of unwrap_property_ref without default argument.
 * @tparam T            Type of the default value. Access will be attempted with `std::decay_t<T>`.
 * @tparam U            A mutable reference of a property
 * @param ppt           Target property object
 * @param default_val   A default value to assign
 * @throws              Throws BadPropertyCast if the typeid of the requested T does not match that of the contents of ppt.
 * @return              Property value with type std::remove_cvref_t<T>&
 */
template<class T>
std::remove_cvref_t<T>& unwrap_property_ref(Property& ppt, std::weak_ptr<T> default_val) {
  if (not ppt.has_value())
    ppt = std::move(default_val);
  return unwrap_property_ref<std::remove_cvref_t<T>>(ppt);
}

/**
 * @brief Performs type-safe access/assign to the contained property value with a default value
 * @details Similar to property_cast with default argument with the semantics of
 * reference-like access of unwrap_property_ref without default argument.
 * @tparam T            Type of the default value. Access will be attempted with `std::decay_t<T>`.
 * @tparam U            A mutable reference of a property
 * @param ppt           Target property object
 * @param default_val   A default value to assign
 * @throws              Throws BadPropertyCast if the typeid of the requested T does not match that of the contents of ppt.
 * @return              Property value with type std::remove_cvref_t<T>&
 */
template<class T>
std::remove_cvref_t<T>& unwrap_property_ref(Property& ppt, std::reference_wrapper<T> default_val) {
  if (not ppt.has_value())
    ppt = std::move(default_val);
  return unwrap_property_ref<std::remove_cvref_t<T>>(ppt);
}


template<class T>
std::ranges::view auto Property::vector_view() const & {
  std::size_t bound;
  std::function<T(std::size_t)> transform;
  if (is_vector()) {
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
  if (is_vector()) {
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
  if (is_vector()) {
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
