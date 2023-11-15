#ifndef UTILITY_MACROS_H
#define UTILITY_MACROS_H

// to string, without macro expansion
#define STRINGINFY(s) #s
// to string, with macro expansion
#define XSTRINGINFY(s) STRINGINFY(s)

#define INVALIDARG_EXP_IFN(condition, message)                                 \
  if (!(condition)) {                                                          \
    throw std::invalid_argument(message "\nFailed condition: " #condition);    \
  }

#define RUNTIME_EXP_IFN(condition, message)                                    \
  if (!(condition)) {                                                          \
    throw std::runtime_error(message "\nFailed condition: " #condition);       \
  }

#define REQUIRE_IGRAPH_SUCCESS(igraph_call)                                    \
  if (igraph_call) {                                                           \
    throw std::runtime_error("Failure when calling igraph: " #igraph_call);    \
  }

#define APPROX_EQUAL(value1, value2, eps)                                      \
  value1 + eps >= value2&& value1 - eps <= value2

#define APPROX_WITHIN(value1, lo, hi, eps)                                     \
  value1 + eps >= lo&& value1 - eps <= hi

#define XOR(value1, value2) !(value1) != !(value2)

#define SQUARE(expr) ((expr)*(expr))

#endif
