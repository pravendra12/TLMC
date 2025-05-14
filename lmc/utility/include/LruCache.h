#ifndef LMC_UITILITY_INCLUDE_LRUCACHE_H_
#define LMC_UITILITY_INCLUDE_LRUCACHE_H_

#include <thread>
#include <mutex>
#include <shared_mutex>
#include <list>
#include <unordered_map>

using namespace std;

template <typename K, class V>
class LruCache
{
public:
  // Constructor
  explicit LruCache(size_t capacity);
  
  // Setter
  void Add(const K key, const V &value);
  
  // Getter
  [[nodiscard]] bool Get(const K key, V &value);
  
  // Strictly only use in a destructor
  [[nodiscard]] size_t GetSizeOfCacheList() const;


private:
  size_t capacity_;
  // cache list of key and value
  list<pair<K, V>> cache_list_{};
  // hashmap of key and iterator of cache list
  unordered_map<K, decltype(cache_list_.begin())> cache_hashmap_{};
  shared_mutex mu_{};
};

#include "LruCache.inl"

#endif // LMC_UITILITY_INCLUDE_LRUCACHE_H_

