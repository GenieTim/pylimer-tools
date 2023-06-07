
#include <random>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <format>
#include <iostream>
#include <sstream>
#include <string>

typedef int64_t msec_t;
#if defined(__WIN32__)
#include <windows.h>
msec_t
currentTimeInMillis(void)
{
  return timeGetTime();
}
#else
#include <sys/time.h>
msec_t
currentTimeInMillis(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (msec_t)tv.tv_sec * 1000 + tv.tv_usec / 1000;
}
#endif

int
main(int argc, char* argv[])
{
  msec_t start;
  msec_t stop;

  long MAXCOUNT = 10000000;

  std::srand(123);

  int* numbers = new int[MAXCOUNT];
  for (int i = 0; i < MAXCOUNT; i++) {
    numbers[i] = std::rand();
  }

  {
    std::string result;
    result.reserve(MAXCOUNT * 100);
    start = currentTimeInMillis();
    for (int i = 0; i < MAXCOUNT; i++) {
      std::format_to(
        std::back_inserter(result), "Number {} is great!", numbers[i]);
    }
    stop = currentTimeInMillis();
    assert(result.size() > MAXCOUNT);
    std::cout << "timing fmt: " << stop - start
              << " ms   /   string length: " << result.size() << std::endl;
  }

  {
    std::string result;
    result.reserve(MAXCOUNT * 100);
    start = currentTimeInMillis();
    for (int i = 0; i < MAXCOUNT; i++) {
      result += "Number " + std::to_string(numbers[i]) + " is great!";
    }
    stop = currentTimeInMillis();
    assert(result.size() > MAXCOUNT);
    std::cout << "timing to_string: " << stop - start
              << " ms   /   string length: " << result.size() << std::endl;
  }

  {
    std::string result;
    result.reserve(MAXCOUNT * 100);
    std::ostringstream ss;
    start = currentTimeInMillis();
    for (int i = 0; i < MAXCOUNT; i++) {
      ss << "Number " << numbers[i] << " is great!";
    }
    result = ss.str();
    stop = currentTimeInMillis();
    assert(result.size() > MAXCOUNT);
    std::cout << "timing stds: " << stop - start
              << " ms   /   string length: " << result.size() << std::endl;
  }
}