#pragma once

#include "./utilityMacros.h"
#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

#ifdef OPENMP_FOUND
#include <omp.h>
#endif

namespace pylimer_tools::utils {

template<unsigned int S>
class PerformanceTimer
{

private:
  std::array<std::string, S> sectionNames;
  std::array<double, S> sectionMean_mus;
  std::array<double, S> sectionVariance_mus;
  std::array<int, S> sectionNumMeasurements;
  int currentSection = -1;
  std::chrono::time_point<std::chrono::high_resolution_clock> startTime;

  const int col1width = 30;
  const int col2width = 10;
  const int col3width = 15;
  const int col4width = 15;
  const int col5width = 12;
  const int unitwidth = 5;

  template<typename T>
  void printNumberRightWithWidth(T t,
                                 const std::string& unit,
                                 const int width,
                                 const int precision) const
  {
    std::cout << std::right << std::setw(width - unitwidth) << std::setfill(' ')
              << std::setprecision(precision) << t;
    // add unit
    std::cout << std::left << std::setw(unitwidth) << std::setfill(' ')
              << (" [" + unit + "]");
  }

  template<typename T>
  void printLeftWithWidth(T t, const int width) const
  {
    std::cout << std::left << std::setw(width) << std::setfill(' ') << t;
  }

  void printLine() const
  {
    std::cout << std::left
              << std::setw(col1width + col2width + col3width + col4width +
                           col5width)
              << std::setfill('-') << "-";
    std::cout << std::endl;
  }

  void printHeader() const
  {
    this->printLeftWithWidth("Section", col1width);
    this->printLeftWithWidth("Num", col2width);
    this->printLeftWithWidth("Mean time", col3width);
    this->printLeftWithWidth("Std time", col4width);
    this->printLeftWithWidth("[%] of total", col5width);
    std::cout << std::endl;
    this->printLine();
  }

  std::pair<std::string, double> formatTime(const double us_time) const
  {
    std::string unit = "µs";
    double conversionFactor = 1.;
    if (us_time > 1e9) {
      unit = "h";
      conversionFactor = 1. / (1e6 * 60 * 60);
    } else if (us_time > 1e6) {
      unit = "s";
      conversionFactor = 1. / 1e6;
    } else if (us_time > 1e3) {
      unit = "ms";
      conversionFactor = 1. / 1e3;
    }

    return std::make_pair(unit, conversionFactor);
  }

  void printSection(const unsigned int idx, const double total) const
  {
    const std::pair<std::string, double> timeConversion =
      this->formatTime(total);
    std::string unit = timeConversion.first;
    double conversionFactor = timeConversion.second;

    this->printLeftWithWidth(this->sectionNames[idx], this->col1width);
    this->printLeftWithWidth(this->sectionNumMeasurements[idx],
                             this->col2width);
    this->printNumberRightWithWidth(
      this->sectionMean_mus[idx] * conversionFactor, unit, this->col3width, 5);
    this->printNumberRightWithWidth(this->sectionVariance_mus[idx] *
                                      conversionFactor,
                                    unit,
                                    this->col4width,
                                    5);
    this->printNumberRightWithWidth(
      100 * ((static_cast<double>(this->sectionNumMeasurements[idx]) *
              this->sectionMean_mus[idx]) /
             total),
      "%",
      this->col5width,
      2);
    std::cout << std::endl;
  }

public:
  PerformanceTimer()
  {
    this->sectionMean_mus.fill(0.);
    this->sectionVariance_mus.fill(0.);
    this->sectionNumMeasurements.fill(0);
    this->sectionNames.fill("");
  }

  /**
   * @brief Set the section name
   *
   * @param idx
   * @param name
   */
  void registerSection(const unsigned int idx, std::string name)
  {
    static_assert(idx < S);
    this->sectionNames[idx] = name;
  }

  /**
   * @brief Set all section names
   *
   * @param sections
   */
  void registerSections(const std::array<std::string, S>& sections)
  {
    this->sectionNames = sections;
  }

  /**
   * @brief Switch measurement context to the specified section
   *
   * @param idx
   */
  void section(const unsigned int idx)
  {
    this->stop();
    this->start(idx);
  }

  /**
   * @brief Stop and record the currently running measurement.
   *
   * Does nothing if no measurement is currently running
   */
  void stop()
  {
    if (this->currentSection == -1) {
      // first need to start measurement
      return;
    }
    std::chrono::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration duration = end - this->startTime;
    long int us =
      (std::chrono::duration_cast<std::chrono::microseconds>(duration)).count();
    const double numMeasurementsBefore =
      static_cast<double>(this->sectionNumMeasurements[this->currentSection]);
    this->sectionNumMeasurements[this->currentSection] += 1;
    double denominator =
      1. /
      static_cast<double>(this->sectionNumMeasurements[this->currentSection]);
    // record the new variance/standard deviation
    // see e.g.
    // https://math.stackexchange.com/questions/102978/incremental-computation-of-standard-deviation#comment241843_103025
    this->sectionVariance_mus[this->currentSection] =
      (numMeasurementsBefore > 0
         ? (((numMeasurementsBefore - 1) / (numMeasurementsBefore)) *
            this->sectionVariance_mus[this->currentSection])
         : 0) +
      denominator * SQUARE(us - this->sectionMean_mus[this->currentSection]);

    // record the new average microseconds
    this->sectionMean_mus[this->currentSection] =
      (this->sectionMean_mus[this->currentSection]) *
        (numMeasurementsBefore * denominator) +
      us * denominator;

    this->currentSection = -1;
  }

  /**
   * @brief Start measuring, discarding anything that may have been already
   * started
   *
   * @param sectionIdx the section, if any
   */
  void start(const unsigned int sectionIdx = 0)
  {
    // static_assert(sectionIdx < S);
    this->currentSection = static_cast<int>(sectionIdx);
    this->startTime = std::chrono::high_resolution_clock::now();
  }

  /**
   * @brief Output the results of all the measurements
   *
   */
  void output()
  {
    this->printHeader();

    double total = 0.0;
    long int numMeasurements = 0;
    for (unsigned int i = 0; i < S; ++i) {
      total += this->sectionMean_mus[i] *
               static_cast<double>(this->sectionNumMeasurements[i]);
      numMeasurements += this->sectionNumMeasurements[i];
    }

    for (size_t i = 0; i < S; ++i) {
      this->printSection(i, total);
    }
    this->printLine();

    std::cout << "\n";

    const std::pair<std::string, double> totalTimeConversion =
      this->formatTime(total);
    // std::pair<double, std::string> meanTimeConversion =
    // this->formatTime(total/(static_cast<double>(numMeasurements)));

#ifdef OPENMP_FOUND
    std::cout << "OpenMP: " << omp_get_max_threads() << " threads."
              << std::endl;
#endif
    std::cout << "Total: " << (total * totalTimeConversion.second) << " "
              << totalTimeConversion.first << " for " << numMeasurements
              << " measurements." << std::endl;
  }
};
}
