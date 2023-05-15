#include <dune/pdelab/common/property_tree.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/classname.hh>

#include <unordered_set>
#include <mutex>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <sstream>

namespace Dune::PDELab::inline Experimental {

inline static std::once_flag _property_default_format_flag;
inline static std::unordered_map<std::type_index, std::function<std::string(const Property&)>> _property_format;
inline static std::unordered_map<void const *, std::pair<PropertyTree*,bool>> _ptree_registry;

void register_property_tree(void const * key, PropertyTree& ptree, bool reference) {
  _ptree_registry[key] = {&ptree, reference};
}

void unregister_property_tree(void const * key) {
  _ptree_registry.erase(key);
}

void Property::register_ptree(PropertyTree& ptree, bool reference) {
  register_property_tree(this, ptree, reference);
  _clean_up = [this]{ unregister_property_tree(this); };
}


bool Property::has_property_tree() const {
  const std::type_index& type = object.type();
  return   (type == typeid(PropertyTree))
        or (type == typeid(std::shared_ptr<PropertyTree>))
        or (type == typeid(std::weak_ptr<PropertyTree>))
        or (_ptree_registry.contains(this));
}

const PropertyTree& Property::as_property_tree() const {
  PropertyTree const * ptree = nullptr;
  if (object.type() == typeid(PropertyTree)) {
    ptree = &any_cast<const PropertyTree&>(*this);
  } else if (object.type() == typeid(std::shared_ptr<PropertyTree>)) {
    ptree = any_cast<const std::shared_ptr<PropertyTree>&>(*this).get();
  } else if (object.type() == typeid(std::weak_ptr<PropertyTree>)) {
    if (auto observe_ptr = any_cast<const std::weak_ptr<PropertyTree>&>(*this).lock())
      ptree = observe_ptr.get();
    else
      DUNE_THROW(InvalidStateException, "Reference to sub-ptreeionary does not exist anymore");
  } else if (auto it = _ptree_registry.find(this); it != _ptree_registry.end()) {
    ptree = it->second.first;
  } else {
    DUNE_THROW(RangeError, "Property (" << Impl::demangle(object.type().name()) << ") is not a ptreeionary");
  }
  assert(ptree);
  return *ptree;
}

PropertyTree& Property::as_property_tree() {
  PropertyTree* ptree = nullptr;
  if (not object.has_value())
    *this = std::make_any<PropertyTree>();

  if (object.type() == typeid(PropertyTree)) {
    ptree = &any_cast<PropertyTree&>(*this);
  } else if (object.type() == typeid(std::shared_ptr<PropertyTree>)) {
    ptree = any_cast<std::shared_ptr<PropertyTree>&>(*this).get();
  } else if (object.type() == typeid(std::weak_ptr<PropertyTree>)) {
    if (auto observe_ptr = any_cast<std::weak_ptr<PropertyTree>&>(*this).lock())
      ptree = observe_ptr.get();
    else
      DUNE_THROW(InvalidStateException, "Reference to sub-ptreeionary does not exist anymore");
  } else if (auto it = _ptree_registry.find(this); it != _ptree_registry.end()) {
    ptree = it->second.first;
  } else {
    DUNE_THROW(RangeError, "Property (" << Impl::demangle(object.type().name()) << ") is not a ptreeionary");
  }
  assert(ptree);
  return *ptree;
}

bool Property::has_array() const {
  const std::type_index& type = object.type();
  return   (type == typeid(std::vector<Property>));
}

const std::vector<Property>& Property::as_array() const {
  return any_cast<const std::vector<Property>&>(*this);
}

std::vector<Property>& Property::as_array() {
  return any_cast<std::vector<Property>&>(*this);
}

const Property& Property::operator[](std::size_t i) const {
  return as_array()[i];
}

Property& Property::operator[](std::size_t i) {
  if (not object.has_value())
    *this = std::vector<Property>(i+1);
  auto& vec = as_array();
  if (vec.size() <= i) vec.resize(i+1);
  return vec[i];
}

std::set<PropertyTree const *> Property::report(std::ostream& out, std::string name, std::string indent) const {

  // register default formatters
  std::call_once(_property_default_format_flag, []{
    auto trivial = std::tuple<std::string,
      double, signed char, unsigned char, short int,
      unsigned short int, int, unsigned int, long int,
      unsigned long int, long long int, unsigned long long int>{};
    std::apply([]<class... T>(T...){
      (register_format<T>(), ...);
    }, trivial);
    register_format<void*>();

    register_format<bool>([](const Property& val){
      std::stringstream ss;
      ss << std::boolalpha << any_cast<bool>(val);
      return ss.str();
    });

    register_format<char>([](const Property& val){
      std::stringstream ss;
      ss << "\'" << any_cast<char>(val) << "\'";
      return ss.str();
    });

    register_format<char const*>([](const Property& val){
      std::stringstream ss;
      ss << "\"" << any_cast<char const*>(val) << "\"";
      return ss.str();
    });
    register_format<std::string>([](const Property& val){
      std::stringstream ss;
      ss << "\"" << any_cast<const std::string&>(val) << "\"";
      return ss.str();
    });

    register_format<std::vector<Property>>([](const Property& val){
      const auto& vec = any_cast<const std::vector<Property>&>(val);
      std::stringstream ss;
      ss << "[";
      for(std::size_t i = 0; i != vec.size(); ++i)
        ss << vec[i] << (i+1 == vec.size() ? "]" : ", ");
      return ss.str();
    });
  });

  std::set<PropertyTree const *> refs;
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
  if (has_property_tree()) {
    out << indent + "  " + Dune::className<PropertyTree>();
    if (   object.type() == typeid(std::weak_ptr<PropertyTree>)
        or object.type() == typeid(std::shared_ptr<PropertyTree>)
        or _ptree_registry.contains(this)) {
      refs.insert(&as_property_tree());
      out << " = " << &as_property_tree();
    } else {
      out << " " << name << " ";
      refs.merge(as_property_tree().report(out, indent + "  "));
    }
  } else {
    if (not indent.empty() and not name.empty())
      out << indent + Impl::demangle(object.type().name()) + " " + name + " = ";
    if (object.has_value()) {
      if (auto it = _property_format.find(object.type()); it != _property_format.end())
        out <<  it->second(*this);
      else
        out << "<unregistered-formatter>";
    } else {
      out << "<empty>";
    }
  }
  return refs;
}

void Property::report_bad_cast(const std::string& type) const {
  std::cerr << "=> Bad cast from 'std::any := "
            << Impl::demangle(object.type().name()) << "' to '" << type << "' <=\n";
  if (not documentation.empty())
    std::cerr << "\n************* Documentation *************\n"
              << documentation << "\n"
              << std::string(41, '*') + "\n\n";
}

std::ostream& operator<<(std::ostream& out, const Property& ppt) {
  ppt.report(out, "", "");
  return out;
}

void Property::register_format(const std::type_index& type, std::function<std::string(const Property&)> f) {
  _property_format[type] = f;
}

} // namespace Dune::PDELab::inline Experimental
