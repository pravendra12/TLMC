#include "LruCache.h"

// Implementation of LruCache class

template <typename K, class V>
LruCache<K, V>::LruCache(size_t capacity) : capacity_(capacity)
{}

template <typename K, class V>
void LruCache<K, V>::Add(const K key, const V &value)
{
  unique_lock<shared_mutex> lock(mu_);
  auto it = cache_hashmap_.find(key);

  if (it != cache_hashmap_.end())
  {
    cache_list_.erase(it->second);
  }
  cache_list_.push_front(make_pair(key, value));
  cache_hashmap_[key] = cache_list_.begin();

  if (cache_hashmap_.size() > capacity_)
  {
    auto last = cache_list_.rbegin()->first;
    cache_list_.pop_back();
    cache_hashmap_.erase(last);
  }
}

template <typename K, class V>
bool LruCache<K, V>::Get(const K key, V &value)
{
  unique_lock<shared_mutex> lock(mu_);
  auto it = cache_hashmap_.find(key);

  if (it == cache_hashmap_.end())
  {
    return false;
  }
  
  cache_list_.splice(cache_list_.begin(), cache_list_, it->second);
  value = it->second->second;
  return true;
}

template <typename K, class V>
inline size_t LruCache<K, V>::GetSizeOfCacheList() const
{
  return cache_list_.size();
}
