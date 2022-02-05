#ifndef STRING_UTIL_H
#define STRING_UTIL_H

#include <algorithm>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
// #include <ranges>
// #include <string_view>

namespace std {
static std::string to_string(std::string input) { return input; }
} // namespace std

namespace pylimer_tools {
namespace utils {

static bool isUpper(std::string str) {
  for (size_t i = 0; i < str.length(); ++i) {
    if (!std::isupper(str[i])) {
      return false;
    }
  }
  return true;
}

static inline bool contains(const std::string haystack,
                            const std::string needle) {
  return haystack.find(needle) != std::string::npos;
}

template <class T, class A> T join(const A &begin, const A &end, const T &t) {
  T result;
  for (A it = begin; it != end; it++) {
    if (!result.empty())
      result.append(t);
    result.append(std::to_string(*it));
  }
  // std::cout << result << std::endl;
  return result;
}

static const std::string WHITESPACE = " \n\r\t\f\v";

static inline std::string ltrim(const std::string &s) {
  size_t start = s.find_first_not_of(WHITESPACE);
  return (start == std::string::npos) ? "" : s.substr(start);
}

static inline std::string rtrim(const std::string &s) {
  size_t end = s.find_last_not_of(WHITESPACE);
  return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

static inline std::string trim(const std::string &s) { return rtrim(ltrim(s)); }

static inline std::string rstrip(std::string haystack,
                                 const std::string needle) {
  auto pos = haystack.find(needle);
  if (pos != std::string::npos) {
    haystack.erase(pos);
  }
  return haystack;
}

static inline bool startsWith(const std::string haystack,
                              const std::string needle) {
  return haystack.compare(0, needle.size(), needle) == 0;
}

static inline std::string trimLineOmitComment(std::string line) {
  line = pylimer_tools::utils::ltrim(line);
  // trim comments
  line = pylimer_tools::utils::rstrip(line, "#");
  return line;
}

static inline std::string trimLineOmitComment(char *line) {
  std::string tempString = std::string(line);
  return pylimer_tools::utils::trimLineOmitComment(tempString);
}

class CsvTokenizer {

public:
  CsvTokenizer(const std::string subject) {
    this->source = subject;
    // Either use C++ 20 implementation, if <ranges> is available
    // std::vector<std::string> results;
    // constexpr std::string_view words{subject};
    // constexpr std::string_view delim{" ,;\t\n"};
    // for (auto word : std::views::split(words, delim))
    // {
    //   results.push_back(word);
    // }

    // or the "manual" one below
    std::string separators = " ,;\t\n";
    std::string text = subject;
    size_t start = text.find_first_not_of(separators);
    size_t end;
    do {
      end = text.find_first_of(separators, start);

      if (end == std::string::npos) {
        std::string token = text.substr(start);
        start = end;
        this->results.push_back(token);
        break;
      }

      std::string token = text.substr(start, end - start);
      this->results.push_back(token);
      start = text.find_first_not_of(separators, end + 1);
    } while (start != std::string::npos);
  }

  CsvTokenizer(const std::string subject, const size_t maxNrToRead) {
    this->source = subject;
    // Either use C++ 20 implementation, if <ranges> is available
    // std::vector<std::string> results;
    // constexpr std::string_view words{subject};
    // constexpr std::string_view delim{" ,;\t\n"};
    // for (auto word : std::views::split(words, delim))
    // {
    //   results.push_back(word);
    // }

    // or the "manual" one below
    this->results.reserve(maxNrToRead);
    std::string separators = " ,;\t\n";
    std::string text = subject;
    size_t start = text.find_first_not_of(separators);
    size_t end;
    size_t iteration = 0;
    do {
      end = text.find_first_of(separators, start);

      if (end == std::string::npos) {
        std::string token = text.substr(start);
        start = end;
        // this->results[iteration] = token;
        this->results.push_back(token);
        break;
      }

      std::string token = text.substr(start, end - start);
      // this->results[iteration] = token;
      this->results.push_back(token);
      iteration++;
      if (iteration == maxNrToRead) {
        break;
      }
      start = text.find_first_not_of(separators, end + 1);
    } while (start != std::string::npos);
  }

  size_t getLength() const { return this->results.size(); }

  template <typename OUT> inline OUT get(size_t index) const {
    return dynamic_cast<OUT>(this->results[index]);
  };

private:
  std::string source;
  std::vector<std::string> results;
};

  template <> inline std::string CsvTokenizer::get<std::string>(size_t index) const {
    return this->results[index];
  };

#define MAKE_GET(TYPE, METHOD)                                                 \
  template <> inline TYPE CsvTokenizer::get<TYPE>(size_t index) const {                             \
    if (this->results.size() <= index) {                                       \
      throw std::runtime_error("Index " + std::to_string(index) +              \
                               " out of bounds when parsing string '" +        \
                               this->source + "'");                            \
    }                                                                          \
    try {                                                                      \
      return METHOD(this->results[index]);                                     \
    } catch (std::invalid_argument e) {                                        \
      throw std::runtime_error("Failed to convert string '" +                  \
                               this->results[index] + "' to " #TYPE);          \
    }                                                                          \
  }

  MAKE_GET(double, std::stod)
  MAKE_GET(long double, std::stold)
  MAKE_GET(float, std::stof)
  MAKE_GET(int, std::stoi)
  MAKE_GET(long int, std::stol)
  MAKE_GET(unsigned int, std::stoul)
  MAKE_GET(unsigned long int, std::stoull)
} // namespace utils

} // namespace pylimer_tools

#endif
