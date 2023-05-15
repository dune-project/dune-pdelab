#ifndef DUNE_PDELAB_COMMON_DICTIONARY_HH
#define DUNE_PDELAB_COMMON_DICTIONARY_HH

#include <dune/common/classname.hh>

#include <unordered_map>
#include <any>
#include <memory>
#include <string>
#include <string_view>
#include <vector>
#include <typeindex>
#include <functional>
#include <iostream>
#include <sstream>
#include <set>
#include <mutex>

namespace Dune::PDELab::inline Experimental {

/**
 * @brief Structure with dynamically assigned data members
 * @details This structure holds a dictionary of data members and its assigned value.
 * Additionally, it allows for recursive definition of dictionaries.
 */
class Dictionary {
public:

  Dictionary& sub(std::string_view key);
  const Dictionary& sub(std::string_view key) const;

  bool hasKey(std::string_view key) const;
  bool hasSub(std::string_view key) const;

  std::any& get(std::string_view key);
  const std::any& get(std::string_view key) const;

  std::any& operator[] (std::string_view key);
  const std::any& operator[] (std::string_view key) const;

  template <class T>
  T& get(std::string_view key) {
    return std::any_cast<T&>(get(key));
  }

  template <class T>
  const T& get(std::string_view key) const {
    auto& val = get(key);
    try {
      return std::any_cast<const T&>(val);
    } catch (...) {
      std::cerr << "Bad cast of key '" << key << "' from 'std::any := "
                << Impl::demangle(val.type().name()) << "' to '" << className<const T&>() << "'\n";
      throw;
    }
  }

  std::vector<std::string_view> keys() const;
  std::vector<std::string_view> subs() const;

  void erase(std::string_view key);

  friend std::ostream& operator<<(std::ostream& out, const Dictionary& dict);

  template<class T, class Proj>
  requires std::is_invocable_r_v<Dictionary&, Proj, T&>
  static std::shared_ptr<Dictionary> make_shared(const std::shared_ptr<T>& ref, Proj proj) {
    return {ref, &std::invoke(proj, *ref)};
  }

  template<class T>
  requires std::is_base_of_v<Dictionary, T>
  static std::shared_ptr<Dictionary> make_shared(const std::shared_ptr<T>& ref) {
    return make_shared(ref, std::identity{});
  }


  template<class T, class Proj>
  requires std::is_invocable_r_v<Dictionary&, Proj, T&>
  static std::weak_ptr<Dictionary> make_reference(const std::shared_ptr<T>& ref, Proj proj) {
    return make_shared(ref, proj);
  }

  template<class T>
  requires std::is_base_of_v<Dictionary, T>
  static std::weak_ptr<Dictionary> make_reference(const std::shared_ptr<T>& ref) {
    return make_shared(ref);
  }

  template<class T>
  void static register_format(std::function<std::string(const std::any&)> f) {
    _format[typeid(T)] = f;
  }

  template<class T>
  void static register_format() {
    _format[typeid(T)] = [](const std::any& val) {
      std::stringstream ss;
      ss << std::any_cast<const T&>(val);
      return ss.str();
    };

    _format[typeid(T*)] = [](const std::any& val){
      std::stringstream ss;
      auto ptr = std::any_cast<T*>(val);
      ss << "[" << ptr << "]";
      if (ptr) ss << " " << *ptr;
      return ss.str();
    };

    if constexpr (not std::is_const_v<T>)
      register_format<const T>();
  }

private:
  std::set<Dictionary const *> report(std::ostream& out, std::string indent = "") const;

  inline static std::once_flag _default_format_flag;
  inline static std::unordered_map<std::type_index, std::function<std::string(const std::any&)>> _format;
  std::unordered_map<std::string, std::any> _param;
};

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_DICTIONARY_HH
