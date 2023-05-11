#include "../../src/pylimer_tools_cpp/entities/EigenNeighbourList.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <random>
#include <vector>

namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;

template<typename... Args>
std::string
string_format(const std::string& format, Args... args)
{
  int size_s = std::snprintf(nullptr, 0, format.c_str(), args...) +
               1; // Extra space for '\0'
  if (size_s <= 0) {
    throw std::runtime_error("Error during formatting.");
  }
  auto size = static_cast<size_t>(size_s);
  std::unique_ptr<char[]> buf(new char[size]);
  std::snprintf(buf.get(), size, format.c_str(), args...);
  return std::string(buf.get(),
                     buf.get() + size - 1); // We don't want the '\0' inside
}

void
runBenchmarkWithConfig(double neighbourBinSize = 1.0,
                       double cutoff = 1.0,
                       double density = 3.,
                       double scatter = 5.)
{
  double boxLen = 12.5;
  pe::Box box = pe::Box(boxLen, boxLen, boxLen);
  int numAtoms = density * box.getVolume();
  std::uniform_real_distribution<double> dist(0, boxLen);
  std::random_device rd;
  std::mt19937 gen = std::mt19937(rd());
  // setup the neighbour list
  Eigen::VectorXd coordinates =
    Eigen::VectorXd::Random(3 * numAtoms) * boxLen * scatter;
  // random [?] access atom list
  // std::cout << "Coordinates set up" << std::endl;
  std::vector<size_t> ids;
  ids.reserve(numAtoms);
  for (size_t i = 0; i < numAtoms; ++i) {
    ids.push_back(i);
  }
  // std::cout << "Ids set up" << std::endl;
  std::shuffle(ids.begin(), ids.end(), gen);
  // std::cout << "Ids shuffled" << std::endl;

  // actually benchmark
  // char benchmarkName[75];
  // int nchar = sprintf(benchmarkName,
  //                     "EigenNeighList_con_%lf_query_%lf_rho_%lf_s_%lf",
  //                     neighbourBinSize,
  //                     cutoff,
  //                     density,
  //                     scatter);
  std::string benchmarkName =
    string_format("EigenNeighList_con_%lf_query_%lf_rho_%lf_s_%lf",
                  neighbourBinSize,
                  cutoff,
                  density,
                  scatter);
  BENCHMARK(benchmarkName.c_str())
  {
    // std::cout << "Setting up Neighlist" << std::endl;
    pe::EigenNeighbourList neighbourlist =
      pe::EigenNeighbourList(coordinates, box, neighbourBinSize);
    // std::cout << "Neighlist set up" << std::endl;

    double meanNrOfNeighs = 0.0;
    // pre-allocate the neighbor indices array
    Eigen::ArrayXi neighbors = Eigen::ArrayXi(static_cast<int>(
      numAtoms * (std::ceil((3.1 * cutoff) * (3.1 * cutoff) * (3.1 * cutoff)) /
                  box.getVolume())));
    for (int step = 0; step < 50; ++step) {
      // std::cout << "Step " << step << std::endl;
      Eigen::Vector3d pairdistance;
      for (size_t i : ids) {
        int numNeighbors = neighbourlist.getIndicesCloseToCoordinates(
          neighbors, coordinates.segment(3 * i, 3), cutoff);
        int actualNumNeighs = 0;
        // pair forces
        for (size_t neigh_idx = 0; neigh_idx < numNeighbors; ++neigh_idx) {
          const size_t j = neighbors[neigh_idx];
          if (j <= i) {
            continue;
          }
          pairdistance =
            coordinates.segment(3 * i, 3) - coordinates.segment(3 * j, 3);
          box.handlePBC(pairdistance);
          const double rNorm = pairdistance.norm();
          if (rNorm >= cutoff || rNorm < 1e-12) {
            continue;
          }

          actualNumNeighs += 1;
        }
        meanNrOfNeighs +=
          static_cast<double>(actualNumNeighs) / static_cast<double>(numAtoms);
      }
      coordinates = Eigen::VectorXd::Random(numAtoms * 3) * boxLen * scatter;
      neighbourlist.resetCoordinates(coordinates);
    }
    return meanNrOfNeighs;
  };
}

TEST_CASE("Eigen Neighbourlist Benchmark",
          "[benchmark][EigenNeighbourList][long]")
{
  for (double outerCutoff = 0.5; outerCutoff <= 2.5; outerCutoff += 0.5) {
    for (double innerCutoff = 0.5 * outerCutoff; innerCutoff <= 2 * outerCutoff;
         innerCutoff +=
         0.2 * 1.5 * outerCutoff) { // 5 iterations per outer cut-off
      runBenchmarkWithConfig(innerCutoff, outerCutoff, 3., 5.);
    }
  }
  REQUIRE(true);
}
