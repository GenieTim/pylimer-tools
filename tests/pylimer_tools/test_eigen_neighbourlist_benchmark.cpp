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

void
runBenchmarkWithConfig(double neighbourBinSize = 1.0,
                       double cutoff = 1.0,
                       double density = 3.,
                       double scatter = 5.)
{
  double boxLen = 20.;
  pe::Box box = pe::Box(boxLen, boxLen, boxLen);
  int numAtoms = density * box.getVolume();
  std::uniform_real_distribution<double> dist(0, boxLen);
  std::random_device rd;
  std::mt19937 gen = std::mt19937(rd());
  // setup the neighbour list
  Eigen::VectorXd coordinates =
    Eigen::VectorXd::Random(3*numAtoms) * boxLen * scatter;
  // random [?] access atom list
  std::vector<size_t> ids; ids.reserve(numAtoms);
  for (size_t i = 0; i < numAtoms; ++i) {
    ids.push_back(i);
  }
  std::shuffle(ids.begin(), ids.end(), gen);

  // actually benchmark
  BENCHMARK("EigenNeighbourList, construct" + std::to_string(neighbourBinSize) +
            "query" + std::to_string(cutoff) + "rho" + std::to_string(density) +
            "s" + std::to_string(scatter))
  {
    pe::EigenNeighbourList neighbourlist =
      pe::EigenNeighbourList(coordinates, box, neighbourBinSize);

    double meanNrOfNeighs = 0.0;
    // pre-allocate the neighbor indices array
    Eigen::ArrayXi neighbors = Eigen::ArrayXi(static_cast<int>(
      numAtoms * (std::ceil((3.1 * cutoff) * (3.1 * cutoff) * (3.1 * cutoff)) /
                  box.getVolume())));
    for (int step = 0; step < 250; ++step) {
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
      coordinates = Eigen::VectorXd::Random(numAtoms*3) * boxLen * scatter;
      neighbourlist.resetCoordinates(coordinates);
    }
    return meanNrOfNeighs;
  };
}

TEST_CASE("Eigen Neighbourlist Benchmark", "[benchmark][EigenNeighbourList][long]")
{
  for (double outerCutoff = 0.5; outerCutoff <= 2.5; outerCutoff += 0.5) {
    for (double innerCutoff = 0.1; innerCutoff <= 2 * outerCutoff;
         innerCutoff += 0.1) {
      runBenchmarkWithConfig(innerCutoff, outerCutoff, 3., 5.);
    }
  }
}