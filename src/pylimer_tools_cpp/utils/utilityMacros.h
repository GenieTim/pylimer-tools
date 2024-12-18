#ifndef UTILITY_MACROS_H
#define UTILITY_MACROS_H

// #include <iostream>

// to string, without macro expansion
#define STRINGINFY(s) #s
// to string, with macro expansion
#define XSTRINGINFY(s) STRINGINFY(s)

// raise exceptions under condition – similar to assert, but kept when compiling
// for any optimisation
#define INVALIDARG_EXP_IFN(condition, message)                                 \
  if (!(condition)) {                                                          \
  //std::cerr << "Argument error: " << message << std::endl;                    \
    throw std::invalid_argument(                                               \
      std::string(message) + std::string("\nFailed condition: " #condition));  \
  }

#define RUNTIME_EXP_IFN(condition, message)                                    \
  if (!(condition)) {                                                          \
  //std::cerr << "Runtime error: " << message << std::endl;                    \
    throw std::runtime_error(std::string(message) +                            \
                             std::string("\nFailed condition: " #condition));  \
  }

#define REQUIRE_IGRAPH_SUCCESS(igraph_call)                                    \
  if (igraph_call) {                                                           \
    throw std::runtime_error("Failure when calling igraph: " #igraph_call);    \
  }

// mathematical closeness
#define APPROX_EQUAL(value1, value2, eps)                                      \
  (((value1 + eps) >= value2) && ((value1 - eps) <= value2))

// NOTE: the "<=" is needed (not "<") to account for 0.
#define APPROX_REL_EQUAL(value1, value2, eps)                                  \
  (std::abs(value1 - value2) <=                                                \
   (eps * std::max(std::abs(value1), std::abs(value2))))

#define APPROX_WITHIN(value1, lo, hi, eps)                                     \
  (((value1 + eps) >= lo) && ((value1 - eps) <= hi))

#define SQUARE(expr) ((expr) * (expr))

#define XOR(value1, value2) !(value1) != !(value2)

// enum as flags
#define MAKE_FLAGS_ENUM(TEnum, TUnder)                                         \
  constexpr TEnum operator~(TEnum a)                                           \
  {                                                                            \
    return static_cast<TEnum>(~static_cast<TUnder>(a));                        \
  }                                                                            \
  constexpr TEnum operator|(TEnum a, TEnum b)                                  \
  {                                                                            \
    return static_cast<TEnum>(static_cast<TUnder>(a) |                         \
                              static_cast<TUnder>(b));                         \
  }                                                                            \
  constexpr TEnum operator&(TEnum a, TEnum b)                                  \
  {                                                                            \
    return static_cast<TEnum>(static_cast<TUnder>(a) &                         \
                              static_cast<TUnder>(b));                         \
  }                                                                            \
  constexpr TEnum operator^(TEnum a, TEnum b)                                  \
  {                                                                            \
    return static_cast<TEnum>(static_cast<TUnder>(a) ^                         \
                              static_cast<TUnder>(b));                         \
  }                                                                            \
  constexpr TEnum& operator|=(TEnum& a, TEnum b)                               \
  {                                                                            \
    a = static_cast<TEnum>(static_cast<TUnder>(a) | static_cast<TUnder>(b));   \
    return a;                                                                  \
  }                                                                            \
  constexpr TEnum& operator&=(TEnum& a, TEnum b)                               \
  {                                                                            \
    a = static_cast<TEnum>(static_cast<TUnder>(a) & static_cast<TUnder>(b));   \
    return a;                                                                  \
  }                                                                            \
  constexpr TEnum& operator^=(TEnum& a, TEnum b)                               \
  {                                                                            \
    a = static_cast<TEnum>(static_cast<TUnder>(a) ^ static_cast<TUnder>(b));   \
    return a;                                                                  \
  }                                                                            \
  constexpr bool operator==(const TEnum& a, const TEnum& b)                    \
  {                                                                            \
    return static_cast<TUnder>(a) == static_cast<TUnder>(b);                   \
  }                                                                            \
  constexpr bool operator==(const TEnum& a, const TUnder& b)                   \
  {                                                                            \
    return static_cast<TUnder>(a) == b;                                        \
  }

#endif
