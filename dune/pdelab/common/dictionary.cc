#include <dune/pdelab/common/dictionary.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/classname.hh>

#include <unordered_set>
#include <mutex>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <sstream>

namespace Dune::PDELab::inline Experimental {


inline static std::once_flag _dictionary_default_format_flag;
inline static std::unordered_map<std::type_index, std::function<std::string(const typename Dictionary::Value&)>> _dictionary_format;
inline static std::unordered_map<typename Dictionary::Value const *, std::pair<Dictionary*,bool>> _dictionary_registry;


bool Dictionary::Value::has_dictionary() const {
  const std::type_index& type = object.type();
  return   (type == typeid(Dictionary))
        or (type == typeid(std::shared_ptr<Dictionary>))
        or (type == typeid(std::weak_ptr<Dictionary>))
        or (_dictionary_registry.contains(this));
}

const Dictionary& Dictionary::Value::as_dictionary() const {
  Dictionary const * dict = nullptr;
  if (object.type() == typeid(Dictionary)) {
    dict = &any_cast<const Dictionary&>(*this);
  } else if (object.type() == typeid(std::shared_ptr<Dictionary>)) {
    dict = any_cast<const std::shared_ptr<Dictionary>&>(*this).get();
  } else if (object.type() == typeid(std::weak_ptr<Dictionary>)) {
    if (auto observe_ptr = any_cast<const std::weak_ptr<Dictionary>&>(*this).lock())
      dict = observe_ptr.get();
    else
      DUNE_THROW(InvalidStateException, "Reference to sub-dictionary does not exist anymore");
  } else if (auto it = _dictionary_registry.find(this); it != _dictionary_registry.end()) {
    dict = it->second.first;
  } else {
    DUNE_THROW(RangeError, "Value (" << Impl::demangle(object.type().name()) << ") is not a dictionary");
  }
  assert(dict);
  return *dict;
}

Dictionary& Dictionary::Value::as_dictionary() {
  Dictionary* dict = nullptr;
  if (not object.has_value())
    *this = std::make_any<Dictionary>();

  if (object.type() == typeid(Dictionary)) {
    dict = &any_cast<Dictionary&>(*this);
  } else if (object.type() == typeid(std::shared_ptr<Dictionary>)) {
    dict = any_cast<std::shared_ptr<Dictionary>&>(*this).get();
  } else if (object.type() == typeid(std::weak_ptr<Dictionary>)) {
    if (auto observe_ptr = any_cast<std::weak_ptr<Dictionary>&>(*this).lock())
      dict = observe_ptr.get();
    else
      DUNE_THROW(InvalidStateException, "Reference to sub-dictionary does not exist anymore");
  } else if (auto it = _dictionary_registry.find(this); it != _dictionary_registry.end()) {
    dict = it->second.first;
  } else {
    DUNE_THROW(RangeError, "Value (" << Impl::demangle(object.type().name()) << ") is not a dictionary");
  }
  assert(dict);
  return *dict;
}

const std::vector<Dictionary::Value>& Dictionary::Value::as_array() const {
  return any_cast<const std::vector<Dictionary::Value>&>(*this);
}

std::vector<Dictionary::Value>& Dictionary::Value::as_array() {
  return any_cast<std::vector<Dictionary::Value>&>(*this);
}

const Dictionary::Value& Dictionary::Value::operator[](std::size_t i) const {
  return as_array()[i];
}

Dictionary::Value& Dictionary::Value::operator[](std::size_t i) {
  if (not object.has_value())
    *this = std::vector<Dictionary::Value>(i+1);
  auto& vec = as_array();
  if (vec.size() <= i) vec.resize(i+1);
  return vec[i];
}

std::set<Dictionary const *> Dictionary::Value::report(std::ostream& out, std::string name, std::string indent) const {

  // register default formatters
  std::call_once(_dictionary_default_format_flag, []{
    auto trivial = std::tuple<std::string,
      double, signed char, unsigned char, short int,
      unsigned short int, int, unsigned int, long int,
      unsigned long int, long long int, unsigned long long int>{};
    std::apply([]<class... T>(T...){
      (register_format<T>(), ...);
    }, trivial);
    register_format<void*>();

    register_format<bool>([](const Dictionary::Value& val){
      std::stringstream ss;
      ss << std::boolalpha << any_cast<bool>(val);
      return ss.str();
    });

    register_format<char>([](const Dictionary::Value& val){
      std::stringstream ss;
      ss << "\'" << any_cast<char>(val) << "\'";
      return ss.str();
    });

    register_format<char const*>([](const Dictionary::Value& val){
      std::stringstream ss;
      ss << "\"" << any_cast<char const*>(val) << "\"";
      return ss.str();
    });
    register_format<std::string>([](const Dictionary::Value& val){
      std::stringstream ss;
      ss << "\"" << any_cast<const std::string&>(val) << "\"";
      return ss.str();
    });

    register_format<std::vector<Value>>([](const Value& val){
      const auto& vec = any_cast<const std::vector<Value>&>(val);
      std::stringstream ss;
      ss << "[";
      for(std::size_t i = 0; i != vec.size(); ++i)
        ss << vec[i] << (i+1 == vec.size() ? "]" : ", ");
      return ss.str();
    });
  });

  std::set<Dictionary const *> refs;
  if (not documentation.empty()) {
    out << indent + "/* ";
    for (std::size_t i = 0; i < documentation.size(); i += 80){
      auto sub = documentation.substr(i, 80);
      if (i != 0) out << indent;
      out << sub;
      if (sub.size() == 80) out << "\n";
    }
    out << " */" << ((documentation.size() > 80) ? '\n' : ' ');
  }
  if (has_dictionary()) {
    out << indent + "  " + Dune::className<Dictionary>();
    if (   object.type() == typeid(std::weak_ptr<Dictionary>)
        or object.type() == typeid(std::shared_ptr<Dictionary>)
        or _dictionary_registry.contains(this)) {
      refs.insert(&as_dictionary());
      out << " = " << &as_dictionary();
    } else {
      out << " " << name << " ";
      refs.merge(as_dictionary().report(out, indent + "  "));
    }
  } else {
    if (not indent.empty() and not name.empty())
      out << indent + Impl::demangle(object.type().name()) + " " + name + " = ";
    if (object.has_value()) {
      if (auto it = _dictionary_format.find(object.type()); it != _dictionary_format.end())
        out <<  it->second(*this);
      else
        out << "<unregistered-formatter>";
    } else {
      out << "<empty>";
    }
  }
  return refs;
}

void Dictionary::Value::debug(std::ostream& out, const std::string& type) const {
  std::cerr << "=> Bad cast from 'std::any := "
            << Impl::demangle(object.type().name()) << "' to '" << type << "' <=\n";
  if (not documentation.empty())
    std::cerr << "\n************* Documentation *************\n"
              << documentation << "\n"
              << std::string(41, '*') + "\n\n";
}

std::ostream& operator<<(std::ostream& out, const Dictionary::Value& value) {
  value.report(out, "", "");
  return out;
}


Dictionary& Dictionary::sub(std::string_view key) {
  std::string_view::size_type dot = key.find(".");
  if (dot != std::string_view::npos)
    return sub(key.substr(0,dot)).sub(key.substr(dot+1));
  else
    return _param[std::string{key}].as_dictionary();
}

const Dictionary& Dictionary::sub(std::string_view key) const {
  std::string_view::size_type dot = key.find(".");
  if (dot != std::string_view::npos)
    return sub(key.substr(0,dot)).sub(key.substr(dot+1));
  else
    return _param.at(std::string{key}).as_dictionary();
}

bool Dictionary::hasKey(std::string_view key) const {
  std::string_view::size_type dot = key.find(".");
  if (dot != std::string_view::npos) {
    std::string_view prefix = key.substr(0,dot);
    if (not _param.contains(std::string{prefix}))
      return false;
    return sub(prefix).hasKey(key.substr(dot+1));
  } else
    return _param.contains(std::string{key});
}


bool Dictionary::hasSub(std::string_view key) const {
  if (not hasKey(key)) return false;
  return get(key).has_dictionary();
}

Dictionary::Value& Dictionary::get(std::string_view key) {
  std::string_view::size_type dot = key.find(".");
  if (dot != std::string_view::npos)
    return sub(key.substr(0,dot)).get(key.substr(dot+1));
  else
    return _param[std::string{key}];
}

const Dictionary::Value& Dictionary::get(std::string_view key) const {
  std::string_view::size_type dot = key.find(".");
  if (dot != std::string_view::npos)
    return sub(key.substr(0,dot)).get(key.substr(dot+1));
  else
    return _param.at(std::string{key});
}

Dictionary::Value& Dictionary::operator[] (std::string_view key) { return get(key); }
const Dictionary::Value& Dictionary::operator[] (std::string_view key) const { return get(key); }

std::vector<std::string_view> Dictionary::keys() const {
  std::vector<std::string_view> keys;
  for(const auto& [key, value] : _param) keys.emplace_back(key);
  return keys;
}

std::vector<std::string_view> Dictionary::subs() const {
  std::vector<std::string_view> keys;
  for(const auto& [key, value] : _param) if (hasSub(key)) keys.emplace_back(key);
  return keys;
}

void Dictionary::erase(std::string_view key) {
  std::string_view::size_type dot = key.find(".");
  if (dot != std::string_view::npos)
    sub(key.substr(0,dot)).erase(key.substr(dot+1));
  else
    _param.erase(std::string{key});
}


std::ostream& operator<<(std::ostream& out, const Dictionary& dict) {
  std::set<Dictionary const *> processed;
  std::set<Dictionary const *> tmp;
  auto refs = dict.report(out);
  processed.insert(&dict);
  while (not refs.empty()) {
    tmp.clear();
    for(auto& ref : refs) if (not processed.contains(ref)) tmp.merge(ref->report(out));
    processed.merge(refs);
    refs = std::move(tmp);
  }
  return out;
}

std::set<Dictionary const *> Dictionary::report(std::ostream& out, std::string indent) const {

  std::set<Dictionary const *> refs;
  out << "[" << this << "] {\n";
  for(const auto& [key, value] : _param) {
    refs.merge(value.report(out, key, indent + "  "));
    out << ";\n";
  }
  out << indent << "};\n";
  return refs;
}

void Dictionary::Value::register_format(const std::type_index& type, std::function<std::string(const Value&)> f) {
  _dictionary_format[type] = f;
}

void Dictionary::register_dictionary(const Value& value, Dictionary& dict, bool reference) {
  _dictionary_registry[&value] = {&dict, reference};
}

void Dictionary::unregister_dictionary(const Value& value) {
  _dictionary_registry.erase(&value);
}

} // namespace Dune::PDELab::inline Experimental
