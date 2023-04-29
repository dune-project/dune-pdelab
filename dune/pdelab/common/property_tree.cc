#include <dune/pdelab/common/property_tree.hh>

#include <cassert>
#include <algorithm>
#include <iostream>
#include <sstream>

namespace Dune::PDELab::inline Experimental {

PropertyTree& PropertyTree::sub(std::string_view key) {
  std::string_view::size_type dot = key.find(".");
  if (dot != std::string_view::npos)
    return sub(key.substr(0,dot)).sub(key.substr(dot+1));
  else
    return _param[std::string{key}].as_property_tree();
}

const PropertyTree& PropertyTree::sub(std::string_view key) const {
  std::string_view::size_type dot = key.find(".");
  if (dot != std::string_view::npos)
    return sub(key.substr(0,dot)).sub(key.substr(dot+1));
  else
    return _param.at(std::string{key}).as_property_tree();
}

bool PropertyTree::hasKey(std::string_view key) const {
  std::string_view::size_type dot = key.find(".");
  if (dot != std::string_view::npos) {
    std::string_view prefix = key.substr(0,dot);
    if (not _param.contains(std::string{prefix}))
      return false;
    return sub(prefix).hasKey(key.substr(dot+1));
  } else
    return _param.contains(std::string{key});
}

bool PropertyTree::hasSub(std::string_view key) const {
  if (not hasKey(key)) return false;
  return get(key).has_property_tree();
}

Property& PropertyTree::get(std::string_view key) {
  std::string_view::size_type dot = key.find(".");
  if (dot != std::string_view::npos)
    return sub(key.substr(0,dot)).get(key.substr(dot+1));
  else
    return _param[std::string{key}];
}

const Property& PropertyTree::get(std::string_view key) const {
  std::string_view::size_type dot = key.find(".");
  if (dot != std::string_view::npos)
    return sub(key.substr(0,dot)).get(key.substr(dot+1));
  else
    return _param.at(std::string{key});
}

Property& PropertyTree::operator[] (std::string_view key) { return get(key); }
const Property& PropertyTree::operator[] (std::string_view key) const { return get(key); }

std::vector<std::string_view> PropertyTree::keys() const {
  std::vector<std::string_view> keys;
  for(const auto& [key, ppt] : _param) keys.emplace_back(key);
  return keys;
}

std::vector<std::string_view> PropertyTree::subs() const {
  std::vector<std::string_view> keys;
  for(const auto& [key, ppt] : _param) if (hasSub(key)) keys.emplace_back(key);
  return keys;
}

void PropertyTree::erase(std::string_view key) {
  std::string_view::size_type dot = key.find(".");
  if (dot != std::string_view::npos)
    sub(key.substr(0,dot)).erase(key.substr(dot+1));
  else
    _param.erase(std::string{key});
}


std::ostream& operator<<(std::ostream& out, const PropertyTree& ptree) {
  std::set<PropertyTree const *> processed;
  std::set<PropertyTree const *> tmp;
  auto refs = ptree.report(out);
  processed.insert(&ptree);
  while (not refs.empty()) {
    tmp.clear();
    for(auto& ref : refs) if (not processed.contains(ref)) tmp.merge(ref->report(out));
    processed.merge(refs);
    refs = std::move(tmp);
  }
  return out;
}

std::set<PropertyTree const *> PropertyTree::report(std::ostream& out, std::string indent) const {

  std::set<PropertyTree const *> refs;
  out << "[" << this << "] {\n";
  for(const auto& [key, ppt] : _param) {
    refs.merge(ppt.report(out, key, indent + "  "));
    out << ";\n";
  }
  out << indent << "};\n";
  return refs;
}

void PropertyTree::report_bad_key(std::string_view key) const {
  std::cerr << "=> Bad cast of key '" << key << "' <=\n";
}

} // namespace Dune::PDELab::inline Experimental
