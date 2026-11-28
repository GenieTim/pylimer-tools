#pragma once

#include <Eigen/Dense>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
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
  const std::chrono::hours hours =
    std::chrono::duration_cast<std::chrono::hours>(duration);
  duration -= hours;
  const std::chrono::minutes minutes =
    std::chrono::duration_cast<std::chrono::minutes>(duration);
  duration -= minutes;
  const std::chrono::seconds seconds =
    std::chrono::duration_cast<std::chrono::seconds>(duration);
  duration -= seconds;
  const std::chrono::milliseconds milliseconds =
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
    if (!std::isupper(static_cast<unsigned char>(str[i]))) {
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

inline constexpr std::string_view WHITESPACE_VIEW = " \n\r\t\f\v";
static const std::string WHITESPACE = " \n\r\t\f\v";

static inline std::string
ltrim(const std::string& s)
{
  const size_t start = s.find_first_not_of(WHITESPACE_VIEW);
  return (start == std::string::npos) ? "" : s.substr(start);
}

static inline std::string
rtrim(const std::string& s)
{
  const size_t end = s.find_last_not_of(WHITESPACE_VIEW);
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
  const size_t pos = haystack.find(needle);
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

static inline std::string_view
trimWhitespaceView(std::string_view sv)
{
  const size_t start = sv.find_first_not_of(WHITESPACE_VIEW);
  if (start == std::string_view::npos) {
    return std::string_view();
  }
  const size_t end = sv.find_last_not_of(WHITESPACE_VIEW);
  return sv.substr(start, end - start + 1);
}

static inline std::vector<std::string_view>
splitWhitespaceView(std::string_view sv)
{
  std::vector<std::string_view> tokens;
  size_t pos = 0;
  while (pos < sv.size()) {
    const size_t start = sv.find_first_not_of(WHITESPACE_VIEW, pos);
    if (start == std::string_view::npos) {
      break;
    }
    const size_t end = sv.find_first_of(WHITESPACE_VIEW, start);
    if (end == std::string_view::npos) {
      tokens.push_back(sv.substr(start));
      break;
    }
    tokens.push_back(sv.substr(start, end - start));
    pos = end + 1;
  }
  return tokens;
}

static inline std::string_view
trimLineOmitCommentView(std::string_view line)
{
  const size_t start = line.find_first_not_of(WHITESPACE_VIEW);
  if (start == std::string_view::npos) {
    return std::string_view();
  }
  line.remove_prefix(start);
  const size_t commentPos = line.find('#');
  if (commentPos != std::string_view::npos) {
    line = line.substr(0, commentPos);
  }
  while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) {
    line.remove_suffix(1);
  }
  return line;
}

static inline std::string
trimLineOmitComment(std::string_view line)
{
  return std::string(trimLineOmitCommentView(line));
}

static inline std::string
trimLineOmitComment(const char* line)
{
  return std::string(trimLineOmitCommentView(std::string_view(line)));
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
  explicit CsvTokenizer(const std::string& subject,
                        const size_t maxNrToRead = 0)
    : source(subject)
  {
    if (maxNrToRead > 0) {
      this->results.reserve(maxNrToRead);
    }
    constexpr std::string_view separators = " ,;\t\r\n";
    const std::string_view text(this->source);
    size_t start = text.find_first_not_of(separators);
    while (start != std::string_view::npos) {
      const size_t end = text.find_first_of(separators, start);
      if (end == std::string_view::npos) {
        this->results.emplace_back(text.substr(start));
        break;
      }
      this->results.emplace_back(text.substr(start, end - start));
      if (maxNrToRead > 0 && this->results.size() >= maxNrToRead) {
        break;
      }
      start = text.find_first_not_of(separators, end + 1);
    }
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
  if (this->results.size() <= index) {
    throw std::runtime_error("Index " + std::to_string(index) +
                             " out of bounds when parsing string '" +
                             this->source + "'");
  }
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
      return static_cast<TYPE>(METHOD(this->results[index]));                  \
    } catch (const std::invalid_argument&) {                                   \
      throw std::runtime_error("Failed to convert string '" +                  \
                               this->results[index] + "' to " #TYPE);          \
    } catch (const std::out_of_range&) {                                       \
      throw std::runtime_error("Value '" + this->results[index] +              \
                               "' out of range for " #TYPE);                   \
    }                                                                          \
  }

MAKE_GET(double, std::stod)
MAKE_GET(long double, std::stold)
MAKE_GET(float, std::stof)
MAKE_GET(int, std::stoi)
MAKE_GET(long int, std::stol)
MAKE_GET(unsigned int, std::stoul)
MAKE_GET(unsigned long int, std::stoull)

inline bool
isNumericTokenBoundary(const char c)
{
  return c == '\0' || c == ' ' || c == '\t' || c == '\r' || c == '\n' ||
         c == '\f' || c == '\v' || c == '#';
}

template<typename T>
inline bool
parseNumericToken(const char*& p, T& out)
{
  char* endptr = nullptr;
  if constexpr (std::is_same_v<T, float>) {
    const float val = std::strtof(p, &endptr);
    if (endptr == p || !isNumericTokenBoundary(*endptr)) {
      return false;
    }
    out = val;
  } else if constexpr (std::is_same_v<T, long double>) {
    const long double val = std::strtold(p, &endptr);
    if (endptr == p || !isNumericTokenBoundary(*endptr)) {
      return false;
    }
    out = val;
  } else if constexpr (std::is_floating_point_v<T>) {
    const double val = std::strtod(p, &endptr);
    if (endptr == p || !isNumericTokenBoundary(*endptr)) {
      return false;
    }
    out = static_cast<T>(val);
  } else if constexpr (std::is_signed_v<T>) {
    const long long val = std::strtoll(p, &endptr, 10);
    if (endptr == p || !isNumericTokenBoundary(*endptr) ||
        val < static_cast<long long>(std::numeric_limits<T>::min()) ||
        val > static_cast<long long>(std::numeric_limits<T>::max())) {
      return false;
    }
    out = static_cast<T>(val);
  } else {
    const char* scan = p;
    while (*scan == ' ' || *scan == '\t' || *scan == '\r' || *scan == '\n' ||
           *scan == '\f' || *scan == '\v') {
      ++scan;
    }
    if (*scan == '-') {
      return false;
    }
    const unsigned long long val = std::strtoull(scan, &endptr, 10);
    if (endptr == scan || !isNumericTokenBoundary(*endptr) ||
        val > static_cast<unsigned long long>(std::numeric_limits<T>::max())) {
      return false;
    }
    out = static_cast<T>(val);
  }
  p = endptr;
  return true;
}

template<typename OUT>
inline std::vector<OUT>
parseTypesInLine(const std::string& line, const int nToRead = -1)
{
  if (nToRead == 0) {
    return {};
  }
  std::vector<OUT> resultnumbers;
  const CsvTokenizer tokenizer =
    (nToRead > 0) ? CsvTokenizer(line, static_cast<size_t>(nToRead))
                  : CsvTokenizer(line);
  resultnumbers.reserve(tokenizer.getLength());
  for (size_t i = 0; i < tokenizer.getLength(); ++i) {
    resultnumbers.push_back(tokenizer.get<OUT>(i));
  }
  return resultnumbers;
}
} // namespace pylimer_tools::utils
