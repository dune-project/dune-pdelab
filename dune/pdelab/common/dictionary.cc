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
inline static std::unordered_map<std::type_index, std::function<std::string(const std::any&)>> _dictionary_format;
inline static std::unordered_map<typename Dictionary::Value const *, std::pair<Dictionary*,bool>> _dictionary_registry;

Dictionary& Dictionary::sub(std::string_view key) {
  std::string_view::size_type dot = key.find(".");
  if (dot != std::string_view::npos)
    return sub(key.substr(0,dot)).sub(key.substr(dot+1));
  else {
    Dictionary * sub_view = nullptr;
    Value& sub_dict = _param[std::string{key}];
    if (not sub_dict.object.has_value())
      sub_dict = std::make_any<Dictionary>();

    if (sub_dict.object.type() == typeid(Dictionary)) {
      sub_view = &any_cast<Dictionary&>(sub_dict);
    } else if (sub_dict.object.type() == typeid(std::shared_ptr<Dictionary>)) {
      sub_view = any_cast<std::shared_ptr<Dictionary>&>(sub_dict).get();
    } else if (sub_dict.object.type() == typeid(std::weak_ptr<Dictionary>)) {
      if (auto observe_ptr = any_cast<std::weak_ptr<Dictionary>&>(sub_dict).lock())
        sub_view = observe_ptr.get();
      else
        DUNE_THROW(InvalidStateException, "Reference to \"" << key << "\" sub-dictionary does not exist anymore");
    } else if (auto it = _dictionary_registry.find(&sub_dict); it != _dictionary_registry.end()) {
      sub_view = it->second.first;
    } else {
      DUNE_THROW(RangeError, "Value (" << Impl::demangle(sub_dict.object.type().name()) << ") for '" << key << "' key is not a sub-dictionary");
    }
    assert(sub_view);
    return *sub_view;
  }
}

const Dictionary& Dictionary::sub(std::string_view key) const {
  std::string_view::size_type dot = key.find(".");
  if (dot != std::string_view::npos)
    return sub(key.substr(0,dot)).sub(key.substr(dot+1));
  else {
    Dictionary const * sub_view = nullptr;
    const Value& sub_dict = _param.at(std::string{key});
    if (not sub_dict.object.has_value())
      DUNE_THROW(RangeError, "Sub-dictionary \"" << key << "\" does not exist");

    if (sub_dict.object.type() == typeid(Dictionary)) {
      sub_view = &any_cast<const Dictionary&>(sub_dict);
    } else if (sub_dict.object.type() == typeid(std::shared_ptr<Dictionary>)) {
      sub_view = any_cast<const std::shared_ptr<Dictionary>&>(sub_dict).get();
    } else if (sub_dict.object.type() == typeid(std::weak_ptr<Dictionary>)) {
      if (auto observe_ptr = any_cast<const std::weak_ptr<Dictionary>&>(sub_dict).lock())
        sub_view = observe_ptr.get();
      else
        DUNE_THROW(InvalidStateException, "Sub-dictionary \"" << key << "\" does not exist");
    } else if (auto it = _dictionary_registry.find(&sub_dict); it != _dictionary_registry.end()) {
      sub_view = it->second.first;
    } else {
      DUNE_THROW(RangeError, "Key \"" << key << "\" is not a sub-dictionary");
    }
    assert(sub_view);
    return *sub_view;
  }
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
  const Value& value = get(key);
  const std::type_index& type = value.object.type();
  return   (type == typeid(Dictionary))
        or (type == typeid(std::shared_ptr<Dictionary>))
        or (type == typeid(std::weak_ptr<Dictionary>))
        or (_dictionary_registry.contains(&value));
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

  // register default formatters
  std::call_once(_dictionary_default_format_flag, []{
    auto trivial = std::tuple<std::string,
      double, bool, signed char, unsigned char, short int,
      unsigned short int, int, unsigned int, long int,
      unsigned long int, long long int, unsigned long long int>{};
    std::apply([]<class... T>(T...){
      (register_format<T>(), ...);
    }, trivial);
    register_format<void*>();
    register_format<char const*>([](const std::any& val){
      std::stringstream ss;
      ss << "\"" << std::any_cast<char const*>(val) << "\"";
      return ss.str();
    });
    register_format<std::string>([](const std::any& val){
      std::stringstream ss;
      ss << "\"" << std::any_cast<std::string>(val) << "\"";
      return ss.str();
    });
  });

  std::set<Dictionary const *> refs;
  out << "[" << this << "] {\n";
  for(const auto& [key, value] : _param) {
    const std::any& object = value.object;
    if (not value.documentation.empty()) {
      out << indent + "  /* ";
      for (std::size_t i = 0; i < value.documentation.size(); i += 80){
        auto sub = value.documentation.substr(i, 80);
        if (i != 0) out << indent + "  ";
        out << sub;
        if (sub.size() == 80) out << "\n";
      }
      out << " */\n";
    }
    if (hasSub(key)) {
      out << indent + "  " + Dune::className<Dictionary>();
      if (   object.type() == typeid(std::weak_ptr<Dictionary>)
          or object.type() == typeid(std::shared_ptr<Dictionary>)
          or _dictionary_registry.contains(&value)) {
        out << "& " << key;
        refs.insert(&sub(key));
        out << " = " << &sub(key) << ";\n";
      } else {
        out << " " << key << " ";
        refs.merge(sub(key).report(out, indent + "  "));
      }
    } else {
      out << indent + "  " + Impl::demangle(object.type().name()) + " " + key;
      if (object.has_value()) {
        if (auto it = _dictionary_format.find(object.type()); it != _dictionary_format.end())
          out << " = " + it->second(object);
        else
          out << " = <unregistered-formatter>";
      }
      out << ";\n";
    }
  }
  out << indent << "};\n";
  return refs;
}

void Dictionary::register_format(const std::type_index& type, std::function<std::string(const std::any&)> f) {
  _dictionary_format[type] = f;
}

void Dictionary::register_dictionary(const Value& value, Dictionary& dict, bool reference) {
  _dictionary_registry[&value] = {&dict, reference};
}

void Dictionary::unregister_dictionary(const Value& value) {
  _dictionary_registry.erase(&value);
}

} // namespace Dune::PDELab::inline Experimental
