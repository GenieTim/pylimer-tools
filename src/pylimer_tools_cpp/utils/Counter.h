#pragma once

#include "stddef.h"
#include <unordered_map>
#include <vector>

namespace pylimer_tools {
namespace utils {
  class IndexCounter
  {
  private:
    std::vector<int> counts = {};

  public:
    explicit IndexCounter(const int size = 0) { this->counts.resize(size, 0); }

    void increment(const int index)
    {
      if (index >= this->counts.size()) {
        this->counts.resize(index + 1, 0);
      }

      this->counts[index] += 1;
    }

    std::unordered_map<int, int> asMap() const
    {
      std::unordered_map<int, int> result;
      for (size_t i = 0; i < this->counts.size(); ++i) {
        if (this->counts[i] > 0) {
          result[i] = this->counts[i];
        }
      }
      return result;
    }
  };

  struct IntDefaultZero
  {
    int i = 0;

    operator int() const { return this->i; };
  };

  template<class T>
  class Counter
  {
  private:
    std::unordered_map<T, IntDefaultZero> counts;

  public:
    explicit Counter(int size = 0) { this->counts.reserve(size); }

    void increment(const T& value) { this->counts[value].i += 1; }

    std::unordered_map<T, int> asMap() const
    {
      std::unordered_map<T, int> result;
      result.reserve(this->counts.size());
      for (const auto& pair : this->counts) {
        if (pair.second.i > 0) {
          result[pair.first] = pair.second.i;
        }
      }
      return result;
    }
  };
}
}
