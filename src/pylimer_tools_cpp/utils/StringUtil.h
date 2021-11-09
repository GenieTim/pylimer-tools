#ifndef STRING_UTIL_H
#define STRING_UTIL_H

#include <string>
#include <iostream>
#include <cstring>
// #include <ranges>
// #include <string_view>
#include <boost/tokenizer.hpp>

namespace std
{
  static std::string to_string(std::string input)
  {
    return input;
  }
}

namespace pylimer_tools
{
  namespace utils
  {

    static bool isUpper(std::string str)
    {
      for (size_t i = 0; i < str.length(); ++i)
      {
        if (!std::isupper(str[i]))
        {
          return false;
        }
      }
      return true;
    }

    static inline bool contains(const std::string *haystack, const std::string needle)
    {
      return haystack->find(needle) != std::string::npos;
    }

    template <class T, class A>
    T join(const A &begin, const A &end, const T &t)
    {
      T result;
      for (A it = begin;
           it != end;
           it++)
      {
        if (!result.empty())
          result.append(t);
        result.append(std::to_string(*it));
      }
      // std::cout << result << std::endl;
      return result;
    }

    class CsvTokenizer
    {
      typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

    public:
      CsvTokenizer(std::string subject)
      {
        // Either use C++ 20 implementation, if <ranges> is available
        // std::vector<std::string> results;
        // constexpr std::string_view words{subject};
        // constexpr std::string_view delim{" ,;\t\n"};
        // for (auto word : std::views::split(words, delim))
        // {
        //   results.push_back(word);
        // }

        // otherwise, use boost
        boost::char_separator<char> sep{" ,;\t\n"};
        tokenizer tok{subject, sep};
        for (tokenizer::iterator it = tok.begin(); it != tok.end(); ++it)
        {
          // we are only interested in the first value.
          // also, we cannot cast the strings that follow
          this->results.push_back(*it);
        }
      }

      int getLength() const { return this->results.size(); }

      template <typename OUT>
      OUT get(size_t index) const { return (OUT)this->results[index]; }

      template <>
      std::string get<std::string>(size_t index) const
      {
        return this->results[index];
      }

#define MAKE_GET(TYPE, METHOD) \
  template <>                  \
  TYPE get<TYPE>(size_t index) const { return METHOD(this->results[index]); }

      MAKE_GET(double, std::stod)
      MAKE_GET(long double, std::stold)
      MAKE_GET(float, std::stof)
      MAKE_GET(int, std::stoi)
      MAKE_GET(long int, std::stol)
      MAKE_GET(unsigned int, std::stoul)
      MAKE_GET(unsigned long int, std::stoull)

    private:
      std::vector<std::string> results;
    };
  }

}

#endif
