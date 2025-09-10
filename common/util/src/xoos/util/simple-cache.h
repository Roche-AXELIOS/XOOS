#pragma once

#include <functional>
#include <map>
#include <memory>
#include <mutex>

namespace xoos::util {

/**
 * A thread-safe cache that stores weak references to shared_ptr values.
 *
 * The cache automatically loads values using a provided loader function when
 * cache misses occur. Values are stored as weak_ptr to allow automatic cleanup
 * when no external references exist.
 *
 * @tparam K The type of keys used for lookups
 * @tparam V The type of cached values (stored as shared_ptr<V>)
 * @tparam K0 The type used for internal storage keys (defaults to K)
 */
template <typename K, typename V, typename K0 = K>
class SimpleCache {
 public:
  /**
   * Constructor with loader function only.
   * Uses identity transformation for keys (K -> K).
   *
   * @param loader Function to load values when cache misses occur
   */
  explicit SimpleCache(std::function<std::shared_ptr<V>(const K&)> loader) : _loader(std::move(loader)) {
  }

  /**
   * Constructor with key transformation and loader functions.
   *
   * @param key_transform Function to transform lookup keys to storage keys
   * @param loader Function to load values when cache misses occur
   */
  SimpleCache(std::function<K0(K)> key_transform, std::function<std::shared_ptr<V>(const K&)> loader)
      : _key_transform(std::move(key_transform)), _loader(std::move(loader)) {
  }

  /**
   * Retrieves a value from the cache or loads it if not present.
   *
   * This method is thread-safe. If the value exists in cache and is still
   * referenced elsewhere, returns the cached value. If the weak_ptr has
   * expired or the key is not found, loads a new value using the loader.
   *
   * @param key The key to lookup
   * @return Shared pointer to the cached or newly loaded value
   */
  std::shared_ptr<V> Get(const K& key) {
    std::scoped_lock lock(_mutex);

    // Transform the lookup key to storage key
    auto k = _key_transform(key);

    // Check if value exists in cache and is still valid
    if (auto it = _cache.find(k); it != _cache.end()) {
      if (auto ptr = it->second.lock()) {
        return ptr;  // Cache hit - return existing value
      }
    }

    // Cache miss or expired weak_ptr - load new value
    auto value = _loader(key);
    _cache[k] = value;  // Store as weak_ptr
    return value;
  }

 private:
  // Key transformation function (defaults to identity)
  std::function<K0(K)> _key_transform = [](const K& k) { return k; };

  // Value loader function
  std::function<std::shared_ptr<V>(const K&)> _loader;

  // Internal cache storage using weak_ptr for automatic cleanup
  std::map<K0, std::weak_ptr<V>> _cache;

  // Mutex for thread safety when updating the cache
  std::mutex _mutex{};
};

}  // namespace xoos::util
