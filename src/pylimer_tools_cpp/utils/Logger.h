#pragma once
#include <iostream>
#include <string>
#include <vector>

namespace pylimer_tools {
namespace utils {
  class Logger
  {
  public:
    static void log(std::string message)
    {
      if (Logger::logEnabled) {
        std::cout << message << std::endl;
      }
    }
    static void enableLog(const bool enable = true)
    {
      Logger::logEnabled = enable;
    }

  private:
    Logger() { Logger::logEnabled = false; }
    static bool logEnabled = false;
  };
} // namespace utils
} // namespace pylimer_tools
