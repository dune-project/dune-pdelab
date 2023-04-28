#ifndef DUNE_PDELAB_COMMON_DICTIONARY_HH
#define DUNE_PDELAB_COMMON_DICTIONARY_HH

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
#include <iostream>


namespace Dune::PDELab::inline Experimental {

/**
 * @brief Structure with dynamically assigned data members
 * @details This structure holds a dictionary of data members and its assigned value.
 * Additionally, it allows for recursive definition of dictionaries.
 */
class Dictionary {
public:

  class Value {

  public:
    ~Value() {
      if (_clean_up) _clean_up();
    }

    template<class T>
    requires (std::is_base_of_v<Dictionary, T> and not std::same_as<Dictionary,T>)
    std::any& operator=(std::weak_ptr<T> dict_ptr) {
      auto observe_ptr = dict_ptr.lock();
      object = std::move(dict_ptr);
      if (observe_ptr)
        register_dictionary(*observe_ptr, true);
      else
        DUNE_THROW(InvalidStateException, "Weak poiner does not point anywhere");
      if (setter) setter(object);
      return object;
    }

    template<class T>
    requires (std::is_base_of_v<Dictionary, T> and not std::same_as<Dictionary,T>)
    std::any& operator=(std::shared_ptr<T> dict_ptr) {
      Dictionary& dict = *dict_ptr;
      object = std::move(dict_ptr);
      register_dictionary(dict, true);
      if (setter) setter(object);
      return object;
    }

    template<class T>
    requires (std::is_base_of_v<Dictionary, std::remove_cv_t<T>> and not std::same_as<Dictionary,std::remove_cv_t<T>>)
    std::any& operator=(T&& dict_ptr) {
      object = std::forward<T>(dict_ptr);
      register_dictionary(std::any_cast<Dictionary&>(object), false);
      if (setter) setter(object);
      return object;
    }

    template<class T>
    std::any& operator=(T&& t) {
      if (_clean_up) {_clean_up(); _clean_up = nullptr;};
      object = std::forward<T>(t);
      if (setter) setter(object);
      return object;
    }


    template <class T>
    friend T any_cast(Value& val) {
      if (val.getter) val.getter(val.object);
      return std::any_cast<T>(val.object);
    }

    template <class T>
    friend T any_cast(const Value& val) {
      if (val.getter) val.getter(val.object);
      return std::any_cast<T>(val.object);
    }

    template <class T>
    friend T any_cast(Value&& val) {
      if (val.getter) val.getter(val.object);
      return std::any_cast<T>(std::move(val.object));
    }

    template <class T, class U>
    friend T any_cast(Value& val, U&& default_val) {
      if (not val.object.has_value())
        val = std::forward<U>(default_val);
      return any_cast<T>(val);
    }

    std::string documentation;
    std::function<void(const std::any&)> getter;
    std::function<void(std::any&)> setter;

  private:
    friend class Dictionary;

    void register_dictionary(Dictionary& dict, bool reference) {
      if (_clean_up) {_clean_up(); _clean_up = nullptr;};
      Dictionary::register_dictionary(*this, dict, reference);
      _clean_up = [this]{ Dictionary::unregister_dictionary(*this); };
    }

    std::any object;
    std::function<void()> _clean_up;
  };

  Dictionary& sub(std::string_view key);
  const Dictionary& sub(std::string_view key) const;

  bool hasKey(std::string_view key) const;
  bool hasSub(std::string_view key) const;

  Value& get(std::string_view key);
  const Value& get(std::string_view key) const;

  Value& operator[] (std::string_view key);
  const Value& operator[] (std::string_view key) const;

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
    const Value& value = get(key);
    try {
      return any_cast<const T&>(value);
    } catch (std::bad_any_cast& exception) {
      // printing diagnostic message may throw again, so we guard it just in
      // case of something worse than a cast is going on
      try {
        std::cerr << "=> Bad cast of key '" << key << "' from 'std::any := "
                  << Impl::demangle(value.object.type().name()) << "' to '" << className<const T&>() << "' <=\n";
        if (not value.documentation.empty())
          std::cerr << "\n************* Documentation of '" << key << "' *************\n"
                    << value.documentation << "\n"
                    << std::string(key.size()+ 47, '*') + "\n\n";
      } catch (...) {}
      throw exception;
    }
  }

  std::vector<std::string_view> keys() const;
  std::vector<std::string_view> subs() const;

  void erase(std::string_view key);

  friend std::ostream& operator<<(std::ostream& out, const Dictionary& dict);

  void static register_format(const std::type_index& type, std::function<std::string(const std::any&)> f);

  template<class T>
  void static register_format(std::function<std::string(const std::any&)> f) {
    register_format(typeid(T), f);
  }

  template<class T>
  void static register_format() {
    register_format(typeid(T), [](const std::any& val) {
      std::stringstream ss;
      ss << std::any_cast<const T&>(val);
      return ss.str();
    });

    register_format(typeid(T*), [](const std::any& val){
      std::stringstream ss;
      auto ptr = std::any_cast<T*>(val);
      ss << "[" << ptr << "]";
      if (ptr) ss << " " << *ptr;
      return ss.str();
    });

    if constexpr (not std::is_const_v<T>)
      register_format<const T>();
  }

private:
  static void register_dictionary(const Value& node, Dictionary& dict, bool reference);
  static void unregister_dictionary(const Value& node);

  std::set<Dictionary const *> report(std::ostream& out, std::string indent = "") const;

  std::unordered_map<std::string, Value> _param;
};

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_DICTIONARY_HH
