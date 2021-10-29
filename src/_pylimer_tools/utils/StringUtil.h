#ifndef STRING_UTIL_H
#define STRING_UTIL_H

#include <string>

#include <cstring>

namespace pylimer_tools
{
  namespace utils
  {

    static bool isUpper(std::string str)
    {
      for (int i = 0; i < str.length(); ++i)
      {
        if (!std::isupper(str[i]))
        {
          return false;
        }
      }
      return true;
    }

    static inline bool contains(const std::string *haystack, const std::string needle) {
      return haystack->find(needle) != std::string::npos;
    }
  }

}

#endif
