#include <dune/pdelab/common/dictionary.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/classname.hh>

#include <cassert>
#include <algorithm>

namespace Dune::PDELab::inline Experimental {

Dictionary& Dictionary::sub(std::string_view key) {
  std::string_view::size_type dot = key.find(".");
  if (dot != std::string_view::npos)
    return sub(key.substr(0,dot)).sub(key.substr(dot+1));
  else {
    Dictionary * sub_view = nullptr;
    std::any& sub_param = _param[std::string{key}];
    if (not sub_param.has_value())
      sub_param = std::make_any<Dictionary>();

    if (sub_param.type() == typeid(Dictionary)) {
      sub_view = &std::any_cast<Dictionary&>(sub_param);
    } else if (sub_param.type() == typeid(std::shared_ptr<Dictionary>)) {
      sub_view = std::any_cast<std::shared_ptr<Dictionary>&>(sub_param).get();
    } else if (sub_param.type() == typeid(std::weak_ptr<Dictionary>)) {
      if (auto observe_ptr = std::any_cast<std::weak_ptr<Dictionary>&>(sub_param).lock())
        sub_view = observe_ptr.get();
      else
        DUNE_THROW(InvalidStateException, "Reference to \"" << key << "\" sub-dictionary does not exist anymore");
    } else {
      DUNE_THROW(RangeError, "Value (" << Impl::demangle(sub_param.type().name()) << ") for '" << key << "' key is not a sub-dictionary");
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
    const std::any& sub_param = _param.at(std::string{key});
    if (not sub_param.has_value())
      DUNE_THROW(RangeError, "Sub-dictionary \"" << key << "\" does not exist");

    if (sub_param.type() == typeid(Dictionary)) {
      sub_view = &std::any_cast<const Dictionary&>(sub_param);
    } else if (sub_param.type() == typeid(std::shared_ptr<Dictionary>)) {
      sub_view = std::any_cast<const std::shared_ptr<Dictionary>&>(sub_param).get();
    } else if (sub_param.type() == typeid(std::weak_ptr<Dictionary>)) {
      if (auto observe_ptr = std::any_cast<const std::weak_ptr<Dictionary>&>(sub_param).lock())
        sub_view = observe_ptr.get();
      else
        DUNE_THROW(InvalidStateException, "Sub-dictionary \"" << key << "\" does not exist");
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
  auto& sub_param = get(key);
  return   (sub_param.type() == typeid(Dictionary))
        or (sub_param.type() == typeid(std::shared_ptr<Dictionary>))
        or (sub_param.type() == typeid(std::weak_ptr<Dictionary>));
}

std::any& Dictionary::get(std::string_view key) {
  std::string_view::size_type dot = key.find(".");
  if (dot != std::string_view::npos)
    return sub(key.substr(0,dot)).get(key.substr(dot+1));
  else
    return _param[std::string{key}];
}

const std::any& Dictionary::get(std::string_view key) const {
  std::string_view::size_type dot = key.find(".");
  if (dot != std::string_view::npos)
    return sub(key.substr(0,dot)).get(key.substr(dot+1));
  else
    return _param.at(std::string{key});
}

std::any& Dictionary::operator[] (std::string_view key) { return get(key); }
const std::any& Dictionary::operator[] (std::string_view key) const { return get(key); }

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
  std::call_once(_default_format_flag, []{
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
    if (hasSub(key)) {
      out << indent + "  " + Dune::className<Dictionary>();
      if (value.type() == typeid(std::weak_ptr<Dictionary>) or value.type() == typeid(std::shared_ptr<Dictionary>)) {
        out << "& " << key;
        refs.insert(&sub(key));
        out << " = " << &sub(key) << ";\n";
      } else {
        out << " " << key << " ";
        refs.merge(sub(key).report(out, indent + "  "));
      }
    } else {
      out << indent + "  " + Impl::demangle(value.type().name()) + " " + key;
      if (value.has_value()) {
        if (auto it = _format.find(value.type()); it != _format.end())
          out << " = " + it->second(value);
        else
          out << " = <unregistered-formatter>";
      }
      out << ";\n";
    }
  }
  out << indent << "};\n";
  return refs;
}

} // namespace Dune::PDELab::inline Experimental
