#ifndef STRING_UTIL_H
#define STRING_UTIL_H

#include <string>
#include <iostream>
#include <cstring>
#include <algorithm>
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

    static inline bool contains(const std::string haystack, const std::string needle)
    {
      return haystack.find(needle) != std::string::npos;
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

    static const std::string WHITESPACE = " \n\r\t\f\v";

    static inline std::string ltrim(const std::string &s)
    {
      size_t start = s.find_first_not_of(WHITESPACE);
      return (start == std::string::npos) ? "" : s.substr(start);
    }

    static inline std::string rtrim(const std::string &s)
    {
      size_t end = s.find_last_not_of(WHITESPACE);
      return (end == std::string::npos) ? "" : s.substr(0, end + 1);
    }

    static inline std::string trim(const std::string &s)
    {
      return rtrim(ltrim(s));
    }

    static inline std::string rstrip(std::string haystack, const std::string needle)
    {
      auto pos = haystack.find(needle);
      if (pos != std::string::npos)
      {
        haystack.erase(pos);
      }
      return haystack;
    }

    static inline bool startsWith(const std::string haystack, const std::string needle)
    {
      return haystack.compare(0, needle.size(), needle) == 0;
    }

    static inline std::string trimLineOmitComment(std::string line)
    {
      line = pylimer_tools::utils::ltrim(line);
      // trim comments
      line = pylimer_tools::utils::rstrip(line, "#");
      return line;
    }

    static inline std::string trimLineOmitComment(char *line)
    {
      std::string tempString = std::string(line);
      return pylimer_tools::utils::trimLineOmitComment(tempString);
    }

    class CsvTokenizer
    {
      typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

    public:
      CsvTokenizer(std::string subject)
      {
        this->source = subject;
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

      CsvTokenizer(std::string subject, size_t maxNrToRead)
      {
        this->source = subject;
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
        size_t iteration = 0;
        for (tokenizer::iterator it = tok.begin(); it != tok.end(); ++it)
        {
          // we are only interested in the first value.
          // also, we cannot cast the strings that follow
          this->results.push_back(*it);
          iteration++;
          if (iteration == maxNrToRead)
          {
            break;
          }
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

#define MAKE_GET(TYPE, METHOD)                                                                    \
  template <>                                                                                     \
  TYPE get<TYPE>(size_t index) const                                                              \
  {                                                                                               \
    if (this->results.size() <= index)                                                            \
    {                                                                                             \
      throw std::runtime_error("Index out of bounds when parsing string '" + this->source + "'"); \
    }                                                                                             \
    try                                                                                           \
    {                                                                                             \
      return METHOD(this->results[index]);                                                        \
    }                                                                                             \
    catch (std::invalid_argument e)                                                               \
    {                                                                                             \
      throw std::runtime_error("Failed to convert string " + this->results[index]);               \
    }                                                                                             \
  }

      MAKE_GET(double, std::stod)
      MAKE_GET(long double, std::stold)
      MAKE_GET(float, std::stof)
      MAKE_GET(int, std::stoi)
      MAKE_GET(long int, std::stol)
      MAKE_GET(unsigned int, std::stoul)
      MAKE_GET(unsigned long int, std::stoull)

    private:
      std::string source;
      std::vector<std::string> results;
    };
  }

}

#endif
