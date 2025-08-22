/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.                                
 * @Author: Zhucong Xi                                                          
 * @Date: 6/14/22 12:36 PM                                                      
 * @Last Modified by: pravendra12                                               
 * @Last Modified: 2025-06-01                                       
 ******************************************************************************/

/*! @file LruCache.h
    @brief File contains declaration of Least Recently Used Cache template class
*/

#ifndef LMC_PRED_INCLUDE_LRUCACHE_HPP_
#define LMC_PRED_INCLUDE_LRUCACHE_HPP_

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
  /*! @brief Constructor for LruCache
      @param capacity Capacity of the LruCache
  */
  explicit LruCache(size_t capacity) : capacity_(capacity)
  {
  }

  /*! @brief Adds the key, value pair
   */
  void Add(const K key, const V &value)
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

  /*! @brief Returns the value give key
   */
  [[nodiscard]] bool Get(const K key, V &value)
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

private:
  /*! @brief Total capacity of the LRU Cache
   */
  size_t capacity_;

  /*! @brief Cache list of key and value
   */
  list<pair<K, V>> cache_list_{};

  /*! @brief Hashmap of key and iterator of cache list
   */
  unordered_map<K, decltype(cache_list_.begin())> cache_hashmap_{};
  shared_mutex mu_{};
};
#endif // LMC_PRED_INCLUDE_LRUCACHE_HPP_