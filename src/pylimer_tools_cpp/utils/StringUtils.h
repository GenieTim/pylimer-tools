#pragma once

// #include <ranges>
// #include <string_view>
#include <Eigen/Dense>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace std {
inline static std::string
to_string(std::string input)
{
  return input;
}

inline static std::string
to_string(Eigen::Vector3d input)
{
  return std::to_string(input[0]) + ", " + std::to_string(input[1]) + ", " +
         std::to_string(input[2]);
}

template<typename T>
static std::string
to_string(std::pair<T, T> pair)
{
  return "(" + std::to_string(pair.first) + ", " + std::to_string(pair.second) +
         ")";
}

template<typename Duration>
std::string
duration_to_string(Duration duration)
{
  auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
  duration -= hours;
  auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
  duration -= minutes;
  auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
  duration -= seconds;
  auto milliseconds =
    std::chrono::duration_cast<std::chrono::milliseconds>(duration);

  std::ostringstream oss;
  oss << hours.count() << ":" << std::setw(2) << std::setfill('0')
      << minutes.count() << ":" << std::setw(2) << std::setfill('0')
      << seconds.count() << "." << std::setw(3) << std::setfill('0')
      << milliseconds.count();

  return oss.str();
}
} // namespace std

namespace pylimer_tools::utils {

inline static bool
isUpper(const std::string& str)
{
  for (size_t i = 0; i < str.length(); ++i) {
    if (!std::isupper(str[i])) {
      return false;
    }
  }
  return true;
}

static inline std::string
removeAll(std::string& haystack, std::string needle)
{
  std::string::size_type i = haystack.find(needle);
  while (i != std::string::npos) {
    haystack.erase(i, needle.length());
    i = haystack.find(needle, i);
  }

  return haystack;
}

static inline std::string
removeAllRegex(const std::string& haystack,
               std::string pattern,
               const bool caseInsensitive = false)
{
  const std::regex reg(pattern,
                       caseInsensitive ? std::regex_constants::icase |
                                           std::regex_constants::ECMAScript
                                       : std::regex_constants::ECMAScript);
  return std::regex_replace(haystack, reg, "");
}

static inline bool
contains(const std::string& haystack, const std::string& needle)
{
  return haystack.find(needle) != std::string::npos;
}

template<class T, class A>
T
join(const A& begin, const A& end, const T& t)
{
  T result;
  for (A it = begin; it != end; it++) {
    if (!result.empty()) {
      result.append(t);
    }
    result.append(std::to_string(*it));
  }
  // std::cout << result << std::endl;
  return result;
}

static const std::string WHITESPACE = " \n\r\t\f\v";

static inline std::string
ltrim(const std::string& s)
{
  const size_t start = s.find_first_not_of(WHITESPACE);
  return (start == std::string::npos) ? "" : s.substr(start);
}

static inline std::string
rtrim(const std::string& s)
{
  const size_t end = s.find_last_not_of(WHITESPACE);
  return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

static inline std::string
trim(const std::string& s)
{
  return rtrim(ltrim(s));
}

static inline std::string
rstrip(std::string haystack, const std::string& needle)
{
  const auto pos = haystack.find(needle);
  if (pos != std::string::npos) {
    haystack.erase(pos);
  }
  return haystack;
}

static inline bool
startsWith(const std::string& haystack, const std::string& needle)
{
  return haystack.compare(0, needle.size(), needle) == 0;
}

static inline std::string
trimLineOmitComment(std::string line)
{
  line = pylimer_tools::utils::ltrim(line);
  // trim comments
  line = pylimer_tools::utils::rstrip(line, "#");
  return line;
}

static inline std::string
trimLineOmitComment(char* line)
{
  const std::string tempString = std::string(line);
  return pylimer_tools::utils::trimLineOmitComment(tempString);
}

static inline std::vector<std::string>
split(const std::string& s, const char delim)
{
  std::vector<std::string> result;
  std::stringstream ss(s);
  std::string item;

  while (getline(ss, item, delim)) {
    result.push_back(item);
  }

  return result;
}

static inline int
split(std::vector<std::string>& res,
      std::string s,
      const std::string& delimiter)
{
  size_t pos_start = 0, pos_end, delim_len = delimiter.length();
  std::string token;

  size_t res_idx = 0;
  while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
    token = s.substr(pos_start, pos_end - pos_start);
    pos_start = pos_end + delim_len;
    if (!token.empty()) {
      if (res_idx >= res.size()) {
        res.push_back(token);
      } else {
        res[res_idx] = token;
      }
      res_idx += 1;
    }
  }
  if (res_idx >= res.size()) {
    res.push_back(s.substr(pos_start));
  } else {
    res[res_idx] = s.substr(pos_start);
  }
  res_idx += 1;
  return static_cast<int>(res_idx);
}

class CsvTokenizer
{

public:
  explicit CsvTokenizer(const std::string& subject)
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

    // or the "manual" one below
    const std::string separators = " ,;\t\n";
    const std::string& text = subject;
    size_t start = text.find_first_not_of(separators);
    do {
      const size_t end = text.find_first_of(separators, start);

      if (end == std::string::npos) {
        const std::string token = text.substr(start);
        start = end;
        this->results.push_back(token);
        break;
      }

      std::string token = text.substr(start, end - start);
      this->results.push_back(token);
      start = text.find_first_not_of(separators, end + 1);
    } while (start != std::string::npos);
  }

  CsvTokenizer(const std::string& subject, const size_t maxNrToRead)
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

    // or the "manual" one below
    this->results.reserve(maxNrToRead);
    const std::string separators = " ,;\t\n";
    const std::string& text = subject;
    size_t start = text.find_first_not_of(separators);
    size_t iteration = 0;
    do {
      const size_t end = text.find_first_of(separators, start);

      if (end == std::string::npos) {
        const std::string token = text.substr(start);
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

  [[nodiscard]] size_t getLength() const { return this->results.size(); }

  template<typename OUT>
  inline OUT get(const size_t index) const
  {
    return dynamic_cast<OUT>(this->results[index]);
  };

private:
  std::string source;
  std::vector<std::string> results;
};

template<>
inline std::string
CsvTokenizer::get<std::string>(const size_t index) const
{
  return this->results[index];
};

#define MAKE_GET(TYPE, METHOD)                                                 \
  template<>                                                                   \
  inline TYPE CsvTokenizer::get<TYPE>(const size_t index) const                \
  {                                                                            \
    if (this->results.size() <= index) {                                       \
      throw std::runtime_error("Index " + std::to_string(index) +              \
                               " out of bounds when parsing string '" +        \
                               this->source + "'");                            \
    }                                                                          \
    try {                                                                      \
      return METHOD(this->results[index]);                                     \
    } catch (std::invalid_argument & e) {                                      \
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
} // namespace pylimer_tools::utils
