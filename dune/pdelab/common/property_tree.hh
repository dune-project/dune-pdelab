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
 * @brief Structure with dynamically assigned data members
 * @details This structure holds a ptree of data members and its assigned value.
 * Additionally, it allows for recursive definition of ptreeionaries.
 */
class PropertyTree {
public:

  PropertyTree& sub(std::string_view key);
  const PropertyTree& sub(std::string_view key) const;

  bool hasKey(std::string_view key) const;
  bool hasSub(std::string_view key) const;

  Property& get(std::string_view key);
  const Property& get(std::string_view key) const;

  Property& operator[] (std::string_view key);
  const Property& operator[] (std::string_view key) const;

  template <class T>
  T& get(std::string_view key) & {
    return any_cast<T&>(get(key));
  }

  template <class T>
  T&& get(std::string_view key) && {
    return any_cast<T&&>(get(key));
  }

  template <class T>
  T& get(std::string_view key, T&& default_value) {
    return any_cast<T&>(get(key), std::forward<T>(default_value));
  }

  template <class T>
  const T& get(std::string_view key) const {
    const Property& value = get(key);
    try {
      return any_cast<const T&>(value);
    } catch (std::bad_any_cast& exception) {
      try {
        report_bad_key(key);
      } catch (...) {}
      throw exception;
    }
  }

  std::vector<std::string_view> keys() const;
  std::vector<std::string_view> subs() const;

  void erase(std::string_view key);

  friend std::ostream& operator<<(std::ostream& out, const PropertyTree& ptree);

  std::set<PropertyTree const *> report(std::ostream& out, std::string indent = "") const;

private:
  void report_bad_key(std::string_view key) const;

  std::unordered_map<std::string, Property> _param;
};

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_PROPERTY_TREE_HH
