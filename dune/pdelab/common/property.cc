#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/common/property_tree.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/classname.hh>

#include <unordered_set>
#include <mutex>
#include <shared_mutex>
#include <cassert>
#include <algorithm>
#include <utility>
#include <iostream>
#include <sstream>
#include <unordered_map>

#if __cpp_lib_stacktrace >= 202011L
#include <stacktrace>
#endif

namespace Dune::PDELab::inline Experimental {

inline static std::once_flag _property_default_format_flag;
inline static std::unordered_map<std::type_index, std::function<std::string(const Property&)>> _property_format;

class PropertyTreeRegistry {
public:
  template<class F>
  decltype(auto) unique(F apply) {
    std::unique_lock guard{_mutex};
    return apply(_map);
  }

  template<class F>
  decltype(auto) shared(F apply) {
    std::shared_lock guard{_mutex};
    return apply(std::as_const(_map));
  }

  std::optional<Impl::WeakPropertyTree> find(Property const * key) {
    return shared([key](const auto& map) -> std::optional<Impl::WeakPropertyTree> {
      if (auto it = map.find(key); it != map.end())
        return it->second;
      else
        return std::nullopt;
    });
  }
private:
  std::unordered_map<Property const *, Impl::WeakPropertyTree> _map;
  std::shared_mutex _mutex;
};
inline static PropertyTreeRegistry _ptree_registry;

void Property::register_derived_ptree(std::optional<Impl::WeakPropertyTree> ptree) {
  if (not ptree) return;
  assert(not _clean_up);
  _ptree_registry.unique([&](auto& map){ map[this] = ptree.value(); });
  _clean_up = [this]{
    _ptree_registry.unique([&](auto& map){ map.erase(this); });
  };
}

Property::Property(std::string name_, std::string documentation_, Setter setter_, Getter getter_)
  : name{std::move(name_)}
  , documentation(std::move(documentation_))
  , setter(std::move(setter_))
  , getter(std::move(getter_))
{}

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
  name = other.name;
  documentation = other.documentation;
  getter = other.getter;
  setter = other.setter;
  _object = other._object;
  register_derived_ptree(_ptree_registry.find(&other));
  return *this;
}

Property& Property::operator=(Property& other) {
  return *this = std::as_const(other);
}


Property& Property::operator=(Property&& other) {
  name = std::move(other.name);
  documentation = std::move(other.documentation);
  getter = std::move(other.getter);
  setter = std::move(other.setter);
  _object = std::move(other._object);
  register_derived_ptree(_ptree_registry.find(&other));
  if (auto clean_up = std::move(other._clean_up)) clean_up();
  return *this;
}

Property& Property::operator=(std::nullptr_t) {
  if (auto clean_up = std::move(_clean_up)) clean_up();
  _clean_up = {};
  _object.reset();
  return *this;
}

bool Property::is_tree() const {
  const std::type_index& type = _object.type();
  return    (type == typeid(PropertyTree))
         or (type == typeid(std::shared_ptr<PropertyTree>))
         or (type == typeid(std::weak_ptr<PropertyTree>))
         or (_ptree_registry.find(this));
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
  } else if (auto ptree_opt = _ptree_registry.find(this)) {
    std::visit([&](auto&& weak_ptr){
      if (auto observe_ptr = weak_ptr.lock())
        ptree = observe_ptr.get();
      else
        DUNE_THROW(InvalidStateException, "Reference to sub-PropertyTree does not exist anymore");
    }, ptree_opt.value());
  }
  if (not ptree)
    throw bad_cast_exception(typeid(const PropertyTree&), "Property is not a const Property Tree");
  return *ptree;
}

PropertyTree& Property::as_tree() {
  PropertyTree* ptree = nullptr;
  if (not _object.has_value())
    *this = PropertyTree{};

  if (_object.type() == typeid(PropertyTree)) {
    ptree = &property_cast<PropertyTree&>(*this);
  } else if (_object.type() == typeid(std::shared_ptr<PropertyTree>)) {
    ptree = property_cast<std::shared_ptr<PropertyTree>>(*this).get();
  } else if (_object.type() == typeid(std::weak_ptr<PropertyTree>)) {
    if (auto observe_ptr = property_cast<std::weak_ptr<PropertyTree>>(*this).lock())
      ptree = observe_ptr.get();
    else
      DUNE_THROW(InvalidStateException, "Reference to sub-PropertyTree does not exist anymore");
  } else if (auto ptree_opt = _ptree_registry.find(this)) {
    std::visit([&ptree]<class T>(const T& weak_ptr){
      if constexpr (std::same_as<T, std::weak_ptr<PropertyTree>>) {
        if (auto observe_ptr = weak_ptr.lock())
          ptree = observe_ptr.get();
        else
          DUNE_THROW(InvalidStateException, "Reference to sub-PropertyTree does not exist anymore");
      }
    }, ptree_opt.value());
  }
  if (not ptree)
    throw bad_cast_exception(typeid(PropertyTree&), "Property is not a const Property Tree");
  return *ptree;
}

bool Property::is_vector() const {
  const std::type_index& type = _object.type();
  return   (type == typeid(std::vector<Property>));
}

const std::vector<Property>& Property::as_vector() const & {
  return property_cast<const std::vector<Property>&>(*this);
}

std::vector<Property>& Property::as_vector() & {
  return property_cast(*this, std::vector<Property>{});
}

std::vector<Property>&& Property::as_vector() && {
  return std::move(property_cast(*this, std::vector<Property>{}));
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
  if (is_tree()) {
    out << indent + Dune::className<PropertyTree>() + " " + name + " ";
    if (   _object.type() == typeid(std::weak_ptr<PropertyTree>)
        or _object.type() == typeid(std::shared_ptr<PropertyTree>)
        or _ptree_registry.find(this)) {
      refs.insert(&as_tree());
      out << "= " << &as_tree();
      out << ";\n";
    } else {
      refs.merge(as_tree().report(out, indent + "  "));
    }
  } else {
    if (not indent.empty() and not name.empty())
      out << indent + Dune::Impl::demangle(_object.type().name()) + " " + name + " = ";
    if (_object.has_value()) {
      if (auto it = _property_format.find(_object.type()); it != _property_format.end())
        out <<  it->second(*this);
      else
        out << "<unregistered-formatter>";
    } else {
      out << "<empty>";
    }
    if (not indent.empty() and not name.empty())
      out << ";\n";
  }
  return refs;
}

BadPropertyCast Property::bad_cast_exception(const std::type_info& type, std::optional<std::string> other_info) const {
  std::stringstream ss;
  ss << "\n===========> Bad property cast <===========\n"
     << "Property name          := '" << name << "'\n"
     << "Source type (std::any) := '" << Dune::Impl::demangle(_object.type().name()) << "'\n"
     << "Target type            := '" << Dune::Impl::demangle(type.name()) << "'\n";
  if (other_info)
    ss << "\n--------------- Other Info ---------------\n"
      << other_info.value() << '\n';
#if __cpp_lib_stacktrace >= 202011L
  ss << "\n--------------- StackTrace ---------------\n"
     << std::stacktrace::current() << '\n';
#endif
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
