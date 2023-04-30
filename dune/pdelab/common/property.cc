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
inline static std::unordered_map<Property const *, std::weak_ptr<PropertyTree>> _ptree_registry;

void Property::register_ptree(std::weak_ptr<PropertyTree> ptree) {
  _ptree_registry[this] = std::move(ptree);
  assert(not _clean_up);
  _clean_up = [this]{ _ptree_registry.erase(this); };
}

Property::Property(std::string_view ppt_name)
  : _name{ppt_name} {}

Property::Property(const Property& other) {
  *this = other;
}

Property::Property(Property&& other) {
  *this = std::move(other);
}

Property::~Property() {
  *this = nullptr;
}

Property& Property::operator=(const Property& other) {
  _name = other._name;
  documentation = other.documentation;
  getter = other.getter;
  setter = other.setter;
  _object = other._object;
  if (auto it = _ptree_registry.find(&other); it != _ptree_registry.end())
    register_ptree(it->second);
  return *this;
}

Property& Property::operator=(Property&& other) {
  _name = std::move(other._name);
  documentation = std::move(other.documentation);
  getter = std::move(other.getter);
  setter = std::move(other.setter);
  _object = std::move(other._object);
  if (auto it = _ptree_registry.find(&other); it != _ptree_registry.end())
    register_ptree(it->second);
  if (auto clean_up = std::move(other._clean_up)) clean_up();
  return *this;
}

bool Property::has_tree() const {
  const std::type_index& type = _object.type();
  return   (type == typeid(PropertyTree))
        or (type == typeid(std::shared_ptr<PropertyTree>))
        or (type == typeid(std::weak_ptr<PropertyTree>))
        or (_ptree_registry.contains(this));
}

const PropertyTree& Property::as_tree() const {
  PropertyTree const * ptree = nullptr;
  if (_object.type() == typeid(PropertyTree)) {
    ptree = &property_cast<const PropertyTree&>(*this);
  } else if (_object.type() == typeid(std::shared_ptr<PropertyTree>)) {
    ptree = property_cast<const std::shared_ptr<PropertyTree>&>(*this).get();
  } else if (_object.type() == typeid(std::weak_ptr<PropertyTree>)) {
    if (auto observe_ptr = property_cast<const std::weak_ptr<PropertyTree>&>(*this).lock())
      ptree = observe_ptr.get();
    else
      DUNE_THROW(InvalidStateException, "Reference to sub-PropertyTree does not exist anymore");
  } else if (auto it = _ptree_registry.find(this); it != _ptree_registry.end()) {
    if (auto observe_ptr = it->second.lock())
      ptree = observe_ptr.get();
    else
      DUNE_THROW(InvalidStateException, "Reference to sub-PropertyTree does not exist anymore");
  } else {
    DUNE_THROW(RangeError, "Property (" << Impl::demangle(_object.type().name()) << ") is not a Property Tree");
  }
  assert(ptree);
  return *ptree;
}

PropertyTree& Property::as_tree() {
  PropertyTree* ptree = nullptr;
  if (not _object.has_value())
    *this = PropertyTree{};

  if (_object.type() == typeid(PropertyTree)) {
    ptree = &property_cast<PropertyTree&>(*this);
  } else if (_object.type() == typeid(std::shared_ptr<PropertyTree>)) {
    ptree = property_cast<std::shared_ptr<PropertyTree>&>(*this).get();
  } else if (_object.type() == typeid(std::weak_ptr<PropertyTree>)) {
    if (auto observe_ptr = property_cast<std::weak_ptr<PropertyTree>&>(*this).lock())
      ptree = observe_ptr.get();
    else
      DUNE_THROW(InvalidStateException, "Reference to sub-PropertyTree does not exist anymore");
  } else if (auto it = _ptree_registry.find(this); it != _ptree_registry.end()) {
    if (auto observe_ptr = it->second.lock())
      ptree = observe_ptr.get();
    else
      DUNE_THROW(InvalidStateException, "Reference to sub-PropertyTree does not exist anymore");
  } else {
    DUNE_THROW(RangeError, "Property (" << Impl::demangle(_object.type().name()) << ") is not a PropertyTree");
  }
  assert(ptree);
  return *ptree;
}

bool Property::has_array() const {
  const std::type_index& type = _object.type();
  return   (type == typeid(std::vector<Property>));
}

const std::vector<Property>& Property::as_array() const {
  return property_cast<const std::vector<Property>&>(*this);
}

std::vector<Property>& Property::as_array() {
  return property_cast<std::vector<Property>&>(*this);
}

const Property& Property::operator[](std::size_t i) const {
  return as_array()[i];
}

Property& Property::operator[](std::size_t i) {
  if (not _object.has_value())
    *this = std::vector<Property>(i+1);
  auto& vec = as_array();
  if (vec.size() <= i) vec.resize(i+1);
  return vec[i];
}

std::set<PropertyTree const *> Property::report(std::ostream& out, std::string indent) const {

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
      ss << std::boolalpha << property_cast<bool>(val);
      return ss.str();
    });

    register_format<char>([](const Property& val){
      std::stringstream ss;
      ss << "\'" << property_cast<char>(val) << "\'";
      return ss.str();
    });

    register_format<char const*>([](const Property& val){
      std::stringstream ss;
      ss << "\"" << property_cast<char const*>(val) << "\"";
      return ss.str();
    });
    register_format<std::string>([](const Property& val){
      std::stringstream ss;
      ss << "\"" << property_cast<const std::string&>(val) << "\"";
      return ss.str();
    });

    register_format<std::vector<Property>>([](const Property& val){
      const auto& vec = property_cast<const std::vector<Property>&>(val);
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
  if (has_tree()) {
    out << indent + Dune::className<PropertyTree>() + " " + _name + " ";
    if (   _object.type() == typeid(std::weak_ptr<PropertyTree>)
        or _object.type() == typeid(std::shared_ptr<PropertyTree>)
        or _ptree_registry.contains(this)) {
      refs.insert(&as_tree());
      out << "= " << &as_tree();
      out << ";\n";
    } else {
      refs.merge(as_tree().report(out, indent + "  "));
    }
  } else {
    if (not indent.empty() and not _name.empty())
      out << indent + Impl::demangle(_object.type().name()) + " " + _name + " = ";
    if (_object.has_value()) {
      if (auto it = _property_format.find(_object.type()); it != _property_format.end())
        out <<  it->second(*this);
      else
        out << "<unregistered-formatter>";
    } else {
      out << "<empty>";
    }
    if (not indent.empty() and not _name.empty())
      out << ";\n";
  }
  return refs;
}

BadPropertyCast Property::bad_cast_exception(const std::string& type) const {
  std::stringstream ss;
  ss << "\n===========> Bad property cast <===========\n"
     << "Property name          := '" << _name << "'\n"
     << "Source type (std::any) := '" << Impl::demangle(_object.type().name()) << "'\n"
     << "Target type            := '" << type << "'\n";
  if (not documentation.empty())
    ss << "\n************* Documentation *************\n"
       << documentation << "\n"
       << std::string(41, '*') + "\n\n";
  BadPropertyCast exception{};
  exception.message(ss.str());
  return exception;
}

std::ostream& operator<<(std::ostream& out, const Property& ppt) {
  ppt.report(out, "");
  return out;
}

void Property::register_format(const std::type_info& type, std::function<std::string(const Property&)> f) {
  _property_format[type] = f;
}

} // namespace Dune::PDELab::inline Experimental
