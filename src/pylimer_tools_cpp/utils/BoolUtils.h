#pragma once
#include <cstdarg>
#include <iostream>

template<typename T>
bool
all_equal(const int count, ...)
{
  va_list args;
  va_start(args, count);

  bool result = true;
  T first = va_arg(args, T);
  for (int i = 1; i < count; ++i) {
    T next = va_arg(args, T);
    if (first != next) {
      std::cout << "Discrepancy at i = " << i << ": " << first << " vs. "
                << next << std::endl;
      result = false;
      break;
    }
  }

  va_end(args);

  return result;
}
