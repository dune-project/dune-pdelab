#ifndef DUNE_PDELAB_COMMON_PROPERTY_TREE_HH
#define DUNE_PDELAB_COMMON_PROPERTY_TREE_HH

#include <dune/pdelab/common/property.hh>

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

/**
 * @brief A tree of properties
 * @details This object is equivalent to a map of std::string to Property.
 * A Property may contain a PropertyTree as its value and its contents may be
 * accessed by concatenating the keys with a dot (i.e., `.`).
 *
 * @note Since Property may hold shared pointers to PropertyTrees and still be
 * interpreted as ParameterTree (see Property::as_tree()) this object may
 * technically also represent a generic directed graph of Property.
 * @warning In case that this object contains a directed graph of Property, the
 * caller needs to ensure that no memory leaks are present within the contained
 * data (e.g., two std::shared_ptr<PropertyTree> poiting to each other). To prevent
 * this issue, use std::weak_ptr<PropertyTree> instead.
 */
class PropertyTree {
public:

  /**
   * @brief Access/Assign a mutable PropertyTree with a given key
   * @details Equivalent to `get(key).as_tree()`. Unassigned
   * keys will be assigned and interpreted as ParameterTree.
   *
   * @param key               The keyword to find the property tree
   * @return PropertyTree&    The property tree object at key
   */
  PropertyTree& sub(std::string_view key);

  //! Access a non-mutable PropertyTree with a given key

  /**
   * @brief Access a non-mutable PropertyTree with a given key
   * @details Equivalent to `get(key).as_tree()`. Unassigned keys will throw an
   * exception.
   *
   * @param key               The keyword to find the property tree
   * @return PropertyTree&    The property tree object at key
   */
  const PropertyTree& sub(std::string_view key) const;

  /**
   * @brief Find if a key is present in the map.
   * @note If the key contains a prefix of the form `[prefix.]key`, the result
   * is equivalent to write `sub(prefix).hasKey(key)`.
   *
   * @param key   The keyword to be found
   * @return true   If key is contained in the object
   * @return false  If key is not contained in the object
   */
  bool hasKey(std::string_view key) const;

  /**
   * @brief Find if a key is present in the map and is a PropertyTree.
   * @note If the key contains a prefix of the form `[prefix.]key`, the result
   * is equivalent to write `sub(prefix).hasSub(key)`.
   *
   * @param key   The keyword to be found
   * @return true   If key is contained in the object and is a PropertyTree
   * @return false  If key is not contained in the object and is a PropertyTree
   */
  bool hasSub(std::string_view key) const;

  /**
   * @brief Obtains Property in the map assigned to key.
   * @note If the key contains a prefix of the form `[prefix.]key`, the result
   * is equivalent to write `sub(prefix).get(key)`. Unassigned keys will be
   * default assigned with a property with name of the key.
   *
   * @param key           The keyword to identify the property
   * @return Property&    The property stored by the key.
   */
  Property& get(std::string_view key);

  /**
   * @brief Obtains Property in the map assigned to key.
   * @note If the key contains a prefix of the form `[prefix.]key`, the result
   * is equivalent to write `sub(prefix).get(key)`. Unassigned keys will throw an
   * exception.
   *
   * @param key           The keyword to identify the property
   * @return Property&    The property stored by the key.
   */
  const Property& get(std::string_view key) const;

  //! Equivalent to get(key)
  Property& operator[] (std::string_view key);

  //! Equivalent to get(key)
  const Property& operator[] (std::string_view key) const;

  /**
   * @brief Obtains a mutable reference to the value contained by the proerty at key
   * @details Equivalent to unwrap_property_ref<T>(get(key))
   *
   * @tparam T      Type to interpret property
   * @param key     Keyword to find property
   * @return T&     (Unwrapped) reference to object contained by the property
   */
  template <class T>
  T& get(std::string_view key) & {
    return unwrap_property_ref<T>(get(key));
  }

  /**
   * @brief Obtains a rvalue reference to the value contained by the proerty at key
   * @details Equivalent to unwrap_property_ref<T>(std::move(get(key)))
   *
   * @tparam T      Type to interpret property
   * @param key     Keyword to find property
   * @return T&     (Unwrapped) reference to object contained by the property
   */
  template <class T>
  T&& get(std::string_view key) && {
    return unwrap_property_ref<T>(std::move(get(key)));
  }

  /**
   * @brief Obtains a non-mutable reference to the value contained by the proerty at key
   * @details Equivalent to unwrap_property_ref<const T>(get(key))
   *
   * @tparam T      Type to interpret property
   * @param key     Keyword to find property
   * @return T&     (Unwrapped) reference to object contained by the property
   */
  template <class T>
  const T& get(std::string_view key) const {
    return unwrap_property_ref<const T>(get(key));
  }

  /**
   * @brief Obtains a mutable reference to the value contained by the proerty at key
   * @details If not assigned, the default value will be used to assign the property value.
   * Equivalent to unwrap_property_ref(get(key), std::forward<T>(default_value)).
   *
   * @tparam T              Type to of the default value
   * @param key             Keyword to find property
   * @param default_value   Default value to assign in case of that key is not assigned yet.
   * @return T&             (Unwrapped) reference to object contained by the property of type std::remove_cvref_t<T>&
   */
  template <class T>
  std::remove_cvref_t<T>& get(std::string_view key, T&& default_value) {
    return unwrap_property_ref(get(key), std::forward<T>(default_value));
  }

  //! List of keys contaied in the container
  std::vector<std::string_view> keys() const;

  //! Removes a propery from the container
  void erase(std::string_view key);

  //! Reports the contents of the container into an output stream
  friend std::ostream& operator<<(std::ostream& out, const PropertyTree& ptree);

  //! Reports the contents of the container into an output stream
  std::set<PropertyTree const *> report(std::ostream& out, std::string indent = "") const;

private:
  std::unordered_map<std::string, Property> _param;
};

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_PROPERTY_TREE_HH
