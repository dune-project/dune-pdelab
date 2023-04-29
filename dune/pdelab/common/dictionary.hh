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
#include <optional>


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
    Value& operator=(std::weak_ptr<T> dict_ptr) {
      auto observe_ptr = dict_ptr.lock();
      object = std::move(dict_ptr);
      if (observe_ptr)
        register_dictionary(*observe_ptr, true);
      else
        DUNE_THROW(InvalidStateException, "Weak poiner does not point anywhere");
      if (setter) setter(object);
      return *this;
    }

    template<class T>
    requires (std::is_base_of_v<Dictionary, T> and not std::same_as<Dictionary,T>)
    Value& operator=(std::shared_ptr<T> dict_ptr) {
      Dictionary& dict = *dict_ptr;
      object = std::move(dict_ptr);
      register_dictionary(dict, true);
      if (setter) setter(object);
      return *this;
    }

    template<class T>
    requires (std::is_base_of_v<Dictionary, std::remove_cv_t<T>> and not std::same_as<Dictionary,std::remove_cv_t<T>>)
    Value& operator=(T&& dict_ptr) {
      object = std::forward<T>(dict_ptr);
      register_dictionary(std::any_cast<Dictionary&>(object), false);
      if (setter) setter(object);
      return *this;
    }

    template<class T>
    Value& operator=(T&& t) {
      if (_clean_up) {_clean_up(); _clean_up = nullptr;};
      object = std::forward<T>(t);
      if (setter) setter(object);
      return *this;
    }

    template <class T>
    friend T any_cast(Value& val) {
      if (val.getter) val.getter(val.object);
      try {
        return std::any_cast<T>(val.object);
      } catch (std::bad_any_cast& ex) {
        val.debug(std::cerr, className<T>());
        throw ex;
      }
    }

    template <class T>
    friend T any_cast(const Value& val) {
      if (val.getter) val.getter(val.object);
      try {
        return std::any_cast<T>(val.object);
      } catch (std::bad_any_cast& ex) {
        val.debug(std::cerr, className<T>());
        throw ex;
      }
    }

    template <class T>
    friend T any_cast(Value&& val) {
      if (val.getter) val.getter(val.object);
      try {
        return std::any_cast<T>(std::move(val.object));
      } catch (std::bad_any_cast& ex) {
        val.debug(std::cerr, className<T>());
        throw ex;
      }
    }

    template <class T, class U>
    friend T any_cast(Value& val, U&& default_val) {
      if (not val.object.has_value())
        val = std::forward<U>(default_val);
      return any_cast<T>(val);
    }

    Value& operator[](std::size_t i);
    const Value& operator[](std::size_t i) const;

    std::size_t size() const { return as_array().size(); }

    friend std::ostream& operator<<(std::ostream& out, const Value& dict);

    std::string documentation;
    std::function<void(const std::any&)> getter;
    std::function<void(std::any&)> setter;

  private:
    friend class Dictionary;

    void debug(std::ostream& out, const std::string& type) const;

    bool has_dictionary() const;
    Dictionary& as_dictionary();
    const Dictionary& as_dictionary() const;

    std::vector<Value>& as_array();
    const std::vector<Value>& as_array() const;

    void register_dictionary(Dictionary& dict, bool reference) {
      if (_clean_up) {_clean_up(); _clean_up = nullptr;};
      Dictionary::register_dictionary(*this, dict, reference);
      _clean_up = [this]{ Dictionary::unregister_dictionary(*this); };
    }


    void static register_format(const std::type_index& type, std::function<std::string(const Value&)> f);

    template<class T>
    void static register_format(std::function<std::string(const Value&)> f) {
      register_format(typeid(T), f);
    }

    template<class T>
    void static register_format() {
      register_format(typeid(T), [](const Value& val) {
        std::stringstream ss;
        ss << any_cast<const T&>(val);
        return ss.str();
      });

      register_format(typeid(T*), [](const Value& val){
        std::stringstream ss;
        auto ptr = any_cast<T*>(val);
        ss << "[" << ptr << "]";
        if (ptr) ss << " " << *ptr;
        return ss.str();
      });

      if constexpr (not std::is_const_v<T>)
        register_format<const T>();
    }


    std::set<Dictionary const *> report(std::ostream& out, std::string, std::string indent) const;

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
      try {
        std::cerr << "=> Bad cast of key '" << key << "' <=\n";
      } catch (...) {}
      throw exception;
    }
  }

  std::vector<std::string_view> keys() const;
  std::vector<std::string_view> subs() const;

  void erase(std::string_view key);

  friend std::ostream& operator<<(std::ostream& out, const Dictionary& dict);

private:
  static void register_dictionary(const Value& node, Dictionary& dict, bool reference);
  static void unregister_dictionary(const Value& node);

  std::set<Dictionary const *> report(std::ostream& out, std::string indent = "") const;

  std::unordered_map<std::string, Value> _param;
};

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_COMMON_DICTIONARY_HH
