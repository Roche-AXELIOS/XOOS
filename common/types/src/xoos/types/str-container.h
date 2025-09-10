#pragma once

#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>

// This file contains template implementations of transparent functions that utilize the std::string key
// to enhance performance
// https://rules.sonarsource.com/cpp/tag/since-c++14/RSPEC-6045/

namespace xoos {

// unordered containers with transparent functions
template <typename Tp, typename Alloc = std::allocator<Tp>>
using UnorderedSet = std::unordered_set<Tp, std::hash<Tp>, std::equal_to<>, Alloc>;
using StrUnorderedSet = std::unordered_set<std::string, std::hash<std::string>, std::equal_to<>>;

template <typename Tp, typename Alloc = std::allocator<std::pair<const std::string, Tp>>>
using StrUnorderedMap = std::unordered_map<std::string, Tp, std::hash<std::string>, std::equal_to<std::string>, Alloc>;

// ordered containers with transparent functions
using StrSet = std::set<std::string, std::less<>>;

template <typename Tp, typename Alloc = std::allocator<std::pair<const std::string, Tp>>>
using StrMap = std::map<std::string, Tp, std::less<>, Alloc>;

}  // namespace xoos
