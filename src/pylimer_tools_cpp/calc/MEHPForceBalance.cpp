#include "MEHPForceBalance.h"
#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <map>
#include <nlopt.hpp>
#include <string>
#include <tuple>
#include <vector>

namespace pylimer_tools {
namespace calc {
  namespace mehp {

    /**
     * FORCE RELAXATION
     */
    void MEHPForceBalance::runForceRelaxation(
      long int maxNrOfSteps, // default: 10000
      double xtol,
      long int innerMaxNrOfSteps,
      double innerXtol,
      double innerAlphaTol)
    {
      this->simulationHasRun = true;
      double stress[3][3];

      for (size_t j = 0; j < 3; j++) {
        for (size_t k = 0; k < 3; k++) {
          stress[j][k] = 0.;
        }
      }

      ForceBalanceNetwork net = this->initialConfig;
      const int M = this->universe.getMolecules(crosslinkerType).size();
      const int N = this->universe.getMeanStrandLength(crosslinkerType) + 1;
      const double bM = this->universe.computeMeanBondLength();
      const int f =
        this->universe.determineFunctionalityPerType()[crosslinkerType];
      bool is2D = this->is2D;

      /* array allocation */
      Eigen::VectorXd u = Eigen::VectorXd::Zero(3 * net.nrOfLinks);

      /* force relaxation */
      double maxDistanceMoved = 0.0;
      size_t iterationsDone = 0;
      size_t totalInnerIterationsDone = 0;
      Eigen::VectorXd springPartitions = this->currentSpringPartitionsVec;
      // std::vector<Eigen::ArrayXi> independentVertexSets =
      //   getHeuristicallyIndependentCoordinateSets(&net);
      // this->getIndependentCoordinateSets(&net);
      // { net.springPartCoordinateIndexA, net.springPartCoordinateIndexB };
      std::cout << "Starting force balance procedure "
                // "with " << independentVertexSets.size() << "vertex sets."
                << std::endl;
      do {
        // std::vector<Eigen::ArrayXi> independentVertexSets =
        //   getRandomCoordinateSets(&net);
        maxDistanceMoved = 0.0;
        // first, place cross-links
        // Eigen::VectorXd oneOverSpringPartitions =
        //   this->assembleOneOverSpringPartition(&net, springPartitions);
        // for (Eigen::ArrayXi vertexSet : independentVertexSets) {
        //   maxDistanceMoved =
        //     std::max(maxDistanceMoved,
        //              this->displaceLinksToMeanPosition(
        //                &net, u, oneOverSpringPartitions, vertexSet, 1.0));
        // }

        // maxDistanceMoved =
        //   this->displaceLinksToMeanPosition(&net, u, springPartitions, 1.0);
        // maxDistanceMoved = std::max(
        //   maxDistanceMoved,
        //   this->displaceLinksToMeanPosition(&net, u, springPartitions, 0.5));
        for (size_t link_idx = 0; link_idx < net.nrOfNodes; ++link_idx) {
          double distanceMoved =
            this->displaceToMeanPosition(&net, u, springPartitions, link_idx);
          maxDistanceMoved = std::max(maxDistanceMoved, distanceMoved);
        }
        // then, place slip-link
        for (size_t link_idx = net.nrOfNodes; link_idx < net.nrOfLinks;
             ++link_idx) {
          size_t innerIterationsDone = 0;
          double displacementDone = 0.0;
          double parametrisationChange = 0.0;
          do {
            parametrisationChange =
              this->updateSpringPartition(&net, u, springPartitions, link_idx);
            displacementDone =
              this->displaceToMeanPosition(&net, u, springPartitions, link_idx);
            innerIterationsDone += 1;
          } while (displacementDone > innerXtol &&
                   innerIterationsDone < innerMaxNrOfSteps &&
                   parametrisationChange > innerAlphaTol);
          totalInnerIterationsDone += innerIterationsDone;
        }
        iterationsDone += 1;
        if (iterationsDone % 50 == 0) {
          std::cout << "Iteration " << iterationsDone << " " << maxDistanceMoved
                    << std::endl;
        }
      } while (maxDistanceMoved > xtol && iterationsDone < maxNrOfSteps);

      // query solution & exit reason
      assert(u.size() == 3 * net.nrOfLinks);
      this->currentDisplacements = u;
      this->currentSpringPartitionsVec = springPartitions;
      this->currentSpringDistances =
        this->evaluateSpringDistances(&net, this->currentDisplacements, is2D);

      this->exitReason = iterationsDone == maxNrOfSteps
                           ? ExitReason::MAX_STEPS
                           : ExitReason::X_TOLERANCE;
      this->nrOfStepsDone += iterationsDone;
      std::cout << iterationsDone << " steps done, " << totalInnerIterationsDone
                << " inner iterations. Last max distance moved: "
                << maxDistanceMoved;
    }

    /**
     * @brief Estimate some sets of random vertices
     *
     * This is particularly useful (used) for parallelising the displacements.
     * Note that the returned Eigen::ArrayXi will contain the indices to the
     * coordinates rather than the actual vertex indices.
     *
     * @param net
     * @return std::vector<Eigen::ArrayXi>
     */
    std::vector<Eigen::ArrayXi> MEHPForceBalance::getRandomCoordinateSets(
      ForceBalanceNetwork* net) const
    {
      std::vector<Eigen::ArrayXi> results;
      results.reserve(2);
      Eigen::ArrayXf randomFloats = Eigen::ArrayXf::Random(net->nrOfLinks);
      size_t nrOfPositiveFloats = (randomFloats > 0).count();
      // TODO: implement something better
      Eigen::ArrayXi links1 = Eigen::ArrayXi(nrOfPositiveFloats * 3);
      size_t links1Idx = 0;
      Eigen::ArrayXi links2 =
        Eigen::ArrayXi((net->nrOfLinks - nrOfPositiveFloats) * 3);
      size_t links2Idx = 0;
      for (size_t i = 0; i < net->nrOfLinks; ++i) {
        if (randomFloats[i] > 0) {
          links1[links1Idx * 3] = 3 * i;
          links1[links1Idx * 3 + 1] = 3 * i + 1;
          links1[links1Idx * 3 + 2] = 3 * i + 2;
          links1Idx += 1;
        } else {
          links2[links2Idx * 3] = 3 * i;
          links2[links2Idx * 3 + 1] = 3 * i + 1;
          links2[links2Idx * 3 + 2] = 3 * i + 2;
          links2Idx += 1;
        }
      }
      results.push_back(links1);
      results.push_back(links2);
      return results;
    }

    /**
     * @brief Estimate some sets of independent vertices
     *
     * This is useful (used) for parallelising the displacements.
     * Note that the returned Eigen::ArrayXi will contain the indices to the
     * coordinates rather than the actual vertex indices.
     *
     * Time complexity: ca. O(|v||e|)
     *
     * @param net
     * @return std::vector<Eigen::ArrayXi>
     */
    std::vector<Eigen::ArrayXi>
    MEHPForceBalance::getHeuristicallyIndependentCoordinateSets(
      ForceBalanceNetwork* net) const
    {
      // std::vector<Eigen::ArrayXi> results = { net->springCoordinateIndexA,
      //                                         net->springCoordinateIndexB };
      std::vector<Eigen::ArrayXi> results;
      // global block list: block all indices that are added to a result already
      ArrayXb globalBlocked = ArrayXb::Constant(net->nrOfLinks, false);
      size_t remainingLinks = net->nrOfLinks;
      size_t globalStartingIdx = 0;
      while (remainingLinks > 0) {
        ArrayXb localBlocked = globalBlocked;
        std::vector<size_t> localIndexList;

        size_t localStartingIndex = globalStartingIdx;
        while (localStartingIndex < net->nrOfLinks) {
          while (localBlocked[localStartingIndex]) {
            localStartingIndex += 1;
            if (!(localStartingIndex < net->nrOfLinks)) {
              goto while2exit;
            }
          }
          // add this link to the current results list
          localIndexList.push_back(localStartingIndex);
          remainingLinks -= 1;
          localBlocked[localStartingIndex] = true;
          globalBlocked[localStartingIndex] = true;
          // block the neighbours
          std::vector<size_t> connections =
            net->springIndicesOfLinks[localStartingIndex];
          for (size_t springIdx : connections) {
            std::vector<size_t> springPartners =
              net->linkIndicesOfSprings[springIdx];
            for (size_t i = 0; i < springPartners.size(); i++) {
              if (springPartners[i] == localStartingIndex) {
                if (i > 0) {
                  localBlocked[springPartners[i - 1]] = true;
                }
                if (i < springPartners.size() - 1) {
                  localBlocked[springPartners[i + 1]] = true;
                }
                break;
              }
            }
          }
        }
      while2exit:

        while (globalStartingIdx < net->nrOfLinks &&
               globalBlocked[globalStartingIdx]) {
          globalStartingIdx += 1;
        }

        // translate the localIndexList to the results
        Eigen::ArrayXi localRes =
          Eigen::ArrayXi::Zero(3 * localIndexList.size());
        for (int i = 0; i < localIndexList.size(); i++) {
          localRes.segment(3 * i, 3) << 3 * localIndexList[i],
            3 * localIndexList[i] + 1, 3 * localIndexList[i] + 2;
        }
        results.push_back(localRes);
      }
      // TODO: implement something better
      return results;
    }

    /**
     * @brief Build a graph of the current configuration and find all the sets
     * of independent vertices
     *
     * This is particularly useful (used) for parallelising the displacements.
     * Note that the returned Eigen::ArrayXi will contain the indices to the
     * coordinates rather than the actual vertex indices.
     *
     * @param net
     * @return std::vector<Eigen::ArrayXi>
     */
    std::vector<Eigen::ArrayXi> MEHPForceBalance::getIndependentCoordinateSets(
      ForceBalanceNetwork* net) const
    {
      std::vector<Eigen::ArrayXi> results;
      // build graph
      igraph_t graph;
      igraph_empty(&graph, net->nrOfLinks, IGRAPH_UNDIRECTED);
      igraph_vector_int_t edges;
      igraph_vector_int_init(&edges, net->nrOfPartialSprings * 2);
      for (int i = 0; i < net->nrOfPartialSprings; ++i) {
        // igraph_vector_int_push_back(&edges, net->springPartIndexA[i]);
        // igraph_vector_int_push_back(&edges, net->springPartIndexB[i]);
        igraph_vector_int_set(&edges, 2 * i, net->springPartIndexA[i]);
        igraph_vector_int_set(&edges, 2 * i + 1, net->springPartIndexB[i]);
      }
      assert(igraph_vector_int_size(&edges) == net->nrOfPartialSprings * 2);
      igraph_add_edges(&graph, &edges, nullptr);
      igraph_vector_int_destroy(&edges);

      // find dependencies in graph
      igraph_vector_int_list_t dependencies;
      igraph_vector_int_list_init(&dependencies, 0);
      igraph_independent_vertex_sets(&graph, &dependencies, -1, -1);

      // assemble to results
      results.reserve(igraph_vector_int_list_size(&dependencies));
      for (size_t i = 0; i < igraph_vector_int_list_size(&dependencies); ++i) {
        igraph_vector_int_t* depsI =
          igraph_vector_int_list_get_ptr(&dependencies, i);
        Eigen::ArrayXi result = Eigen::ArrayXi(igraph_vector_int_size(depsI));
        for (size_t j = 0; j < igraph_vector_int_size(depsI); ++j) {
          igraph_integer_t resI = igraph_vector_int_get(depsI, j);
          for (size_t dir = 0; dir < 3; ++dir) {
            result[3 * j + dir] = 3 * resI + dir;
          }
        }
        results.push_back(result);
      }
      igraph_vector_int_list_destroy(&dependencies);

      return results;
    };

    /**
     * @brief Updates the partition/parametrisation of a spring around one link
     *
     */
    double MEHPForceBalance::updateSpringPartition(
      const ForceBalanceNetwork* net,
      const Eigen::VectorXd& u,
      Eigen::VectorXd& springPartitions, /* gives the parametrisation of N */
      const size_t linkIdx) const
    {
      INVALIDARG_EXP_IFN(linkIdx < net->springIndicesOfLinks.size(),
                         "Link to update needs to be in the list");
      INVALIDARG_EXP_IFN(net->linkIsSliplink[linkIdx],
                         "Only slip-links may slip along a spring");
      std::vector<size_t> springIndices = net->springIndicesOfLinks[linkIdx];
      double maxDiff = 0.0;
      for (size_t springIndex : springIndices) {
        std::vector<size_t> springsPartners =
          net->linkIndicesOfSprings[springIndex];
        for (size_t partner_idx = 1; partner_idx < springsPartners.size() - 1;
             ++partner_idx) {
          if (springsPartners[partner_idx] == linkIdx) {
            // found position of this link in this spring
            // want to find the ideal value for
            // net->springPartitions[springIndex][partner_idx-1]
            double distanceBack = ((MEHPForceBalance::evaluateDistanceBetween(
                                      net,
                                      u,
                                      springsPartners[partner_idx],
                                      springsPartners[partner_idx - 1],
                                      this->is2D))
                                     .squaredNorm());
            double distanceForward =
              ((MEHPForceBalance::evaluateDistanceBetween(
                  net,
                  u,
                  springsPartners[partner_idx + 1],
                  springsPartners[partner_idx],
                  this->is2D))
                 .squaredNorm());
            double idealValue =
              1. / (1. + sqrt(distanceForward /
                              distanceBack)); // TODO: above, below?
            if (distanceBack < 1e-9) {
              idealValue = 0.0; // TODO: really?
            }
            size_t currentSpringGlobalIdx =
              net->localToGlobalSpringIndex.at(springIndex)[partner_idx - 1];
            size_t neighbourSpringGlobalIdx =
              net->localToGlobalSpringIndex.at(springIndex)[partner_idx];
            double currentS = springPartitions[currentSpringGlobalIdx];
            double nextS = springPartitions[neighbourSpringGlobalIdx];
            double newS = idealValue * (nextS + currentS);
            maxDiff = std::max(currentS - newS, maxDiff);
            springPartitions[currentSpringGlobalIdx] = newS;
            springPartitions[neighbourSpringGlobalIdx] =
              (1. - idealValue) * (nextS + currentS);
          }
        }
      }
      return maxDiff;
    }

    /**
     * @brief Displace one link to the mean of all connected neighbours
     *
     * @param net the force balance network
     * @param u the current displacements, wherein the resulting coordinates
     * shall be stored
     * @param linkIdx the idx of the link to displace
     * @return double, the distance (squared norm) displaced
     */
    double MEHPForceBalance::displaceToMeanPosition(
      const ForceBalanceNetwork* net,
      Eigen::VectorXd& u,
      const Eigen::VectorXd& springPartitions,
      const size_t linkIdx) const
    {
      std::vector<size_t> springIndices = net->springIndicesOfLinks[linkIdx];
      // Eigen::Vector3d currentDisplacement = u.segment(3 * linkIdx, 3);
      Eigen::Vector3d objectiveDisplacement =
        Eigen::Vector3d::Zero(); // = remainingDisplacement.array();
      double objectiveDisplacementContributors = 0.0;
      for (size_t spring_index = 0; spring_index < springIndices.size();
           ++spring_index) {
        // compute partial distances & total distance of this spring
        std::vector<size_t> springsPartners =
          net->linkIndicesOfSprings[springIndices[spring_index]];
        for (size_t partner_idx = 0; partner_idx < springsPartners.size() - 1;
             ++partner_idx) {
          if (springsPartners[partner_idx] == linkIdx ||
              springsPartners[partner_idx + 1] == linkIdx) {
            Eigen::Vector3d partialDistance;
            // add partial distance to the total distance
            if (springsPartners[partner_idx] == linkIdx) {
              partialDistance = MEHPForceBalance::evaluateDistanceBetween(
                net,
                u,
                springsPartners[partner_idx + 1],
                springsPartners[partner_idx],
                this->is2D);
            } else {
              assert(springsPartners[partner_idx + 1] == linkIdx);
              partialDistance = MEHPForceBalance::evaluateDistanceBetween(
                net,
                u,
                springsPartners[partner_idx],
                springsPartners[partner_idx + 1],
                this->is2D);
            }
            // add to displacement
            double contourLengthFraction =
              springPartitions[net->localToGlobalSpringIndex.at(
                springIndices[spring_index])[partner_idx]];
            if (contourLengthFraction > 1e-12) {
              double oneOverContourLengthFraction =
                1.0 / (net->springsContourLength[springIndices[spring_index]] *
                       contourLengthFraction);
              objectiveDisplacement +=
                (partialDistance)*oneOverContourLengthFraction; // /
              // totalDistance.array());
              objectiveDisplacementContributors += oneOverContourLengthFraction;
            } else {
              objectiveDisplacement = 1e12 * (partialDistance);
              objectiveDisplacementContributors += 1e12;
            }
          }
        }
      }
      // take mean for displacement
      u.segment(3 * linkIdx, 3) +=
        objectiveDisplacement / objectiveDisplacementContributors;

      double dist = objectiveDisplacement.squaredNorm();
      // if (dist > 500) {
      //   std::cout << "Moving " << linkIdx << " for " << dist
      //             << " with displacements " << currentDisplacement[0] << ",
      //             "
      //             << currentDisplacement[1] << ", " <<
      //             currentDisplacement[2]
      //             << std::endl;
      //   std::cout << "For objective displacements "
      //             << objectiveDisplacement[0] << ", "
      //             << objectiveDisplacement[1] << ", "
      //             << objectiveDisplacement[2] << ", for "
      //             << objectiveDisplacementContributors << "." << std::endl;
      // }
      return dist;
    }

    Eigen::Vector3d MEHPForceBalance::evaluateDistanceBetween(
      const ForceBalanceNetwork* net,
      const Eigen::VectorXd& u,
      const size_t linkIndexA,
      const size_t linkIndexB,
      const bool is2D)
    {
      Eigen::Vector3d distances = net->coordinates.segment(3 * linkIndexA, 3) +
                                  u.segment(3 * linkIndexA, 3) -
                                  (net->coordinates.segment(3 * linkIndexB, 3) +
                                   u.segment(3 * linkIndexB, 3));

      // Possibly improvable PBC
      for (size_t j = 0; j < 3; ++j) {
        int iterations = 0;
        assert(!std::isinf(distances[j]) && !std::isnan(distances[j]));
        while (distances[j] > net->boxHalfs[j % 3]) {
          distances[j] -= net->L[j % 3];
          iterations++;
          if (iterations > 10) {
            throw std::runtime_error(
              "Too many iterations in PBC from " + std::to_string(linkIndexA) +
              " to " + std::to_string(linkIndexB) + ", currently at " +
              std::to_string(distances[j]) + " of " +
              std::to_string(net->boxHalfs[j % 3]) + " after " +
              std::to_string(iterations) + " iterations");
          }
        }
        iterations = 0;
        while (distances[j] < -net->boxHalfs[j % 3]) {
          distances[j] += net->L[j % 3];
          iterations++;
          if (iterations > 10) {
            throw std::runtime_error(
              "Too many iterations in PBC from " + std::to_string(linkIndexA) +
              " to " + std::to_string(linkIndexB) + ", currently at " +
              std::to_string(distances[j]) + " of " +
              std::to_string(net->boxHalfs[j % 3]) + " after " +
              std::to_string(iterations) + " iterations");
          }
        }
      }

      if (is2D) {
        distances[2] = 0.0;
      }

      return distances;
    }

    Eigen::VectorXd MEHPForceBalance::evaluateSpringDistances(
      const ForceBalanceNetwork* net,
      const Eigen::VectorXd& u,
      const bool is2D)
    {
      // first, the distances
      assert(u.size() == net->coordinates.size());
      Eigen::VectorXd springDistances =
        Eigen::VectorXd::Zero(3 * net->nrOfSprings);
      Eigen::VectorXd actualCoordinates = net->coordinates + u;

      for (size_t i = 0; i < net->nrOfSprings; ++i) {
        std::vector<size_t> springsPartners = net->linkIndicesOfSprings[i];
        for (size_t partner_idx = 0; partner_idx < springsPartners.size() - 1;
             ++partner_idx) {
          // add partial distance to the total distance
          Eigen::Vector3d partialDistance =
            MEHPForceBalance::evaluateDistanceBetween(
              net,
              u,
              springsPartners[partner_idx],
              springsPartners[partner_idx + 1],
              is2D);
          for (size_t distance_idx = 0; distance_idx < 3; ++distance_idx) {
            springDistances[i * 3 + distance_idx] +=
              partialDistance[distance_idx];
            assert(!isnan(partialDistance[distance_idx]));
          }
        }
      }

      if (is2D) {
        // springDistances(Eigen::seq(2, Eigen::last, Eigen::fix<3>)) =
        //   Eigen::VectorXd::Zero(net->nrOfSprings / 3);
        for (size_t i = 2; i < 3 * net->nrOfSprings; i += 3) {
          springDistances[i] = 0.0;
        }
      }
      assert(springDistances.size() == net->nrOfSprings * 3);

      return springDistances;
    }

    /**
     * FORCE RELAXATION DATA ACCESS
     */
    pylimer_tools::entities::Universe MEHPForceBalance::getCrosslinkerVerse(
      int newCrosslinkerType) const
    {
      // convert nodes & springs back to a universe
      pylimer_tools::entities::Universe xlinkUniverse =
        pylimer_tools::entities::Universe(this->universe.getBox());
      std::vector<long int> ids;
      std::vector<int> types = pylimer_tools::utils::initializeWithValue(
        this->initialConfig.nrOfNodes, crosslinkerType);
      std::vector<double> x;
      std::vector<double> y;
      std::vector<double> z;
      std::vector<int> zeros = pylimer_tools::utils::initializeWithValue(
        this->initialConfig.nrOfNodes, 0);
      ids.reserve(this->initialConfig.nrOfNodes);
      x.reserve(this->initialConfig.nrOfNodes);
      y.reserve(this->initialConfig.nrOfNodes);
      z.reserve(this->initialConfig.nrOfNodes);
      for (int i = 0; i < this->initialConfig.nrOfNodes; ++i) {
        x.push_back(this->initialConfig.coordinates[3 * i + 0] +
                    this->currentDisplacements[3 * i + 0]);
        y.push_back(this->initialConfig.coordinates[3 * i + 1] +
                    this->currentDisplacements[3 * i + 1]);
        z.push_back(this->initialConfig.coordinates[3 * i + 2] +
                    this->currentDisplacements[3 * i + 2]);
        ids.push_back(this->initialConfig.oldAtomIds[i]);
      }
      xlinkUniverse.addAtoms(ids, types, x, y, z, zeros, zeros, zeros);
      std::vector<long int> bondFrom;
      std::vector<long int> bondTo;
      bondFrom.reserve(this->initialConfig.nrOfSprings);
      bondTo.reserve(this->initialConfig.nrOfSprings);
      for (int i = 0; i < this->initialConfig.nrOfSprings; ++i) {
        bondFrom.push_back(
          this->initialConfig.oldAtomIds[this->initialConfig.springIndexA[i]]);
        bondTo.push_back(
          this->initialConfig.oldAtomIds[this->initialConfig.springIndexB[i]]);
      }
      xlinkUniverse.addBonds(
        bondFrom.size(),
        bondFrom,
        bondTo,
        pylimer_tools::utils::initializeWithValue(bondFrom.size(), 1),
        false,
        false); // disable simplify to keep the self-loops etc.
      return xlinkUniverse;
    }

    void MEHPForceBalance::addSlipLinks(const std::vector<size_t> strandIdx1,
                                        const std::vector<size_t> strandIdx2,
                                        const std::vector<double> x,
                                        const std::vector<double> y,
                                        const std::vector<double> z,
                                        const std::vector<double> alpha1,
                                        const std::vector<double> alpha2)
    {
      size_t additionalLen = strandIdx1.size();
      size_t currentNrOfLinks = this->initialConfig.nrOfLinks;
      size_t currentNrOfPartialSprings = this->initialConfig.nrOfPartialSprings;
      if (additionalLen != x.size() || additionalLen != y.size() ||
          additionalLen != z.size()) {
        throw std::invalid_argument("x, y and z must have the same dimensions");
      }
      if (additionalLen != strandIdx2.size() ||
          additionalLen != alpha1.size() || additionalLen != alpha2.size()) {
        throw std::invalid_argument(
          "Strand indices and alpha estimates must have the same length");
      }
      // actually start adding them
      this->initialConfig.nrOfLinks += additionalLen;
      // but first, indicate the resize
      this->initialConfig.coordinates.conservativeResize(
        3 * this->initialConfig.nrOfLinks);
      this->initialConfig.springIndicesOfLinks.reserve(
        this->initialConfig.nrOfLinks);
      this->initialConfig.linkIsSliplink.conservativeResize(
        this->initialConfig.nrOfLinks);
      this->initialConfig.springPartCoordinateIndexA.conservativeResize(
        3 * (currentNrOfPartialSprings + 2 * additionalLen));
      this->initialConfig.springPartCoordinateIndexB.conservativeResize(
        3 * (currentNrOfPartialSprings + 2 * additionalLen));
      this->initialConfig.springPartIndexA.conservativeResize(
        currentNrOfPartialSprings + 2 * additionalLen);
      this->initialConfig.springPartIndexB.conservativeResize(
        currentNrOfPartialSprings + 2 * additionalLen);
      this->currentSpringPartitionsVec.conservativeResize(
        currentNrOfPartialSprings + 2 * additionalLen);
      this->currentDisplacements.conservativeResize(
        3 * this->initialConfig.nrOfLinks);
      size_t partialSpringsAdded = 0;
      // then, loop the slip-links to add
      for (size_t i = 0; i < additionalLen; ++i) {
        // add the info that is straight-forward to add
        this->initialConfig.coordinates[3 * currentNrOfLinks + 3 * i] = x[i];
        this->initialConfig.coordinates[3 * currentNrOfLinks + 3 * i + 1] =
          y[i];
        this->initialConfig.coordinates[3 * currentNrOfLinks + 3 * i + 2] =
          z[i];
        this->initialConfig.linkIsSliplink[currentNrOfLinks + i] = true;
        std::vector<size_t> springIndices{ strandIdx1[i], strandIdx2[i] };
        this->initialConfig.springIndicesOfLinks.push_back(springIndices);
        // add to the springs
        int springIndexIndex = 0;
        for (size_t springIndex : springIndices) {
          std::vector<size_t> springParticipants =
            this->initialConfig.linkIndicesOfSprings[springIndex];
          double alpha = (springIndexIndex == 0 ? alpha1[i] : alpha2[i]);
          // detect the position in the spring
          std::vector<double> partitionsStrand;
          partitionsStrand.reserve(springParticipants.size() - 1);
          for (size_t i = 0; i < springParticipants.size() - 1; ++i) {
            partitionsStrand.push_back(
              this->currentSpringPartitionsVec[this->initialConfig
                                                 .localToGlobalSpringIndex.at(
                                                   springIndex)[i]]);
          }

          bool wasAdded = false;
          size_t targetIndexInSpring = 0;
          double cumulativePartition = 0.0;
          for (size_t p_idx = 0; p_idx < partitionsStrand.size(); ++p_idx) {
            cumulativePartition += partitionsStrand[p_idx];
            if (cumulativePartition > alpha) {
              targetIndexInSpring = p_idx;
              if (p_idx > 0) {
                alpha = alpha - (cumulativePartition - partitionsStrand[p_idx]);
              }
              wasAdded = true;
              break;
            }
          }
          if (!wasAdded) {
            targetIndexInSpring = springParticipants.size() - 1;
            if (partitionsStrand.size() > 0) {
              alpha = alpha - cumulativePartition;
            }
          }

          // TODO: we should not need to be doing this.
          if (alpha < 0. || alpha > 1.) {
            std::cout << "WARNING: alpha = " << alpha << " for spring "
                      << springIndex << std::endl;
            alpha = std::clamp(alpha, 0., 1.);
          }
          // have to adjust the existing springs, too!
          size_t springPartner1 = springParticipants[targetIndexInSpring];
          size_t springPartner2 = springParticipants[targetIndexInSpring + 1];
          size_t newNodeIdx = currentNrOfLinks + i;

          // update connectivity
          size_t lastSpringIndex =
            this->initialConfig.localToGlobalSpringIndex.at(springIndex)
              .at(targetIndexInSpring);
          size_t newSpringIndex =
            currentNrOfPartialSprings + partialSpringsAdded;

          this->initialConfig.localToGlobalSpringIndex.at(springIndex)
            .insert(this->initialConfig.localToGlobalSpringIndex.at(springIndex)
                        .begin() +
                      targetIndexInSpring + 1,
                    newSpringIndex);

          // adjust also the coordinates
          this->currentDisplacements[newNodeIdx] = 0.0;
          if (this->initialConfig.springPartIndexA[lastSpringIndex] ==
              springPartner1) {
            this->initialConfig.springPartIndexB[lastSpringIndex] = newNodeIdx;
            for (size_t offset = 0; offset < 3; ++offset) {
              this->initialConfig
                .springPartCoordinateIndexB[3 * lastSpringIndex + offset] =
                3 * newNodeIdx + offset;
            }
          } else {
            assert(this->initialConfig.springPartIndexA[lastSpringIndex] ==
                   springPartner2);
            this->initialConfig.springPartIndexA[lastSpringIndex] = newNodeIdx;
            for (size_t offset = 0; offset < 3; ++offset) {
              this->initialConfig
                .springPartCoordinateIndexA[3 * lastSpringIndex + offset] =
                3 * newNodeIdx + offset;
            }
          }
          // add the new one
          this->initialConfig.springPartIndexA[newSpringIndex] = newNodeIdx;
          this->initialConfig.springPartIndexB[newSpringIndex] = springPartner2;
          for (size_t offset = 0; offset < 3; ++offset) {
            this->initialConfig
              .springPartCoordinateIndexA[3 * newSpringIndex + offset] =
              3 * newNodeIdx + offset;
            this->initialConfig
              .springPartCoordinateIndexB[3 * newSpringIndex + offset] =
              3 * springPartner2 + offset;
          }

          this->currentSpringPartitionsVec[lastSpringIndex] -= alpha;
          // TODO: we should not need to be doing this.
          if (this->currentSpringPartitionsVec[lastSpringIndex] < 0. ||
              this->currentSpringPartitionsVec[lastSpringIndex] > 1.) {
            std::cout << "WARNING: alpha = "
                      << this->currentSpringPartitionsVec[lastSpringIndex]
                      << " for spring " << lastSpringIndex
                      << " after subtracting alpha = " << alpha << std::endl;
            this->currentSpringPartitionsVec[lastSpringIndex] = std::clamp(
              this->currentSpringPartitionsVec[lastSpringIndex], 0., 1.);
          }
          this->currentSpringPartitionsVec[newSpringIndex] = alpha;
          this->initialConfig.linkIndicesOfSprings[springIndex].insert(
            this->initialConfig.linkIndicesOfSprings[springIndex].begin() +
              targetIndexInSpring +
              1, // + 1 to compensate for the first cross-link
            newNodeIdx);

          partialSpringsAdded += 1;
          springIndexIndex += 1;
        }
      }
      this->initialConfig.nrOfPartialSprings += partialSpringsAdded;
      // do we really want to?
      this->validateNetwork();
      assert(partialSpringsAdded == 2 * additionalLen);
    };

    /**
     * @brief Get the Average Spring Length at the current step
     *
     * @return double
     */
    double MEHPForceBalance::getAverageSpringLength() const
    {
      double r2 = 0.0;
      for (int i = 0; i < this->initialConfig.nrOfSprings; i++) {
        double r2local = 0.0;
        for (int j = 0; j < 3; ++j) {
          r2local += this->currentSpringDistances[i * 3 + j] *
                     this->currentSpringDistances[i * 3 + j];
        }
        r2 += sqrt(r2local);
      }
      return r2 / this->initialConfig.nrOfSprings;
    }

    /**
     * @brief Compute the stress tensor
     *
     * @param net
     * @param u
     * @return std::array<std::array<double, 3>, 3>
     */
    std::array<std::array<double, 3>, 3> MEHPForceBalance::evaluateStressTensor(
      const Eigen::VectorXd& springDistances,
      const double volume) const
    {
      std::array<std::array<double, 3>, 3> stress;

      for (size_t i = 0; i < springDistances.size() / 3; ++i) {
        double s[3] = { springDistances[3 * i + 0],
                        springDistances[3 * i + 1],
                        springDistances[3 * i + 2] };
        /* spring contribution to the overall stress tensor */
        for (size_t j = 0; j < 3; j++) {
          for (size_t k = 0; k < 3; k++) {
            double contribution = this->kappa * s[j] * s[k];
            stress[j][k] += contribution;
          }
        }
      }

      for (size_t j = 0; j < 3; j++) {
        for (size_t k = 0; k < 3; k++) {
          stress[j][k] /= volume;
        }
      }

      return stress;
    }

    /**
     * @brief Compute the stress tensor
     *
     * @param net
     * @param u
     * @param loopTol
     * @return std::array<std::array<double, 3>, 3>
     */
    std::array<std::array<double, 3>, 3> MEHPForceBalance::evaluateStressTensor(
      ForceBalanceNetwork* net,
      const Eigen::VectorXd& u,
      const double loopTol) const
    {
      Eigen::VectorXd springDistances =
        this->evaluateSpringDistances(net, u, this->is2D);

      return this->evaluateStressTensor(springDistances, net->vol);
    }

    std::array<std::array<double, 3>, 3> MEHPForceBalance::getStressTensor()
      const
    {
      return this->evaluateStressTensor(this->currentSpringDistances,
                                        this->initialConfig.vol);
    }

    /**
     * @brief Get the Effective Functionality Of each node
     *
     * Returns the number of active springs connected to each atom, atomId
     * used as index
     *
     * @param tolerance the tolerance: springs under a certain length are
     * considered inactive
     * @return std::unordered_map<long int, int>
     */
    std::unordered_map<long int, int>
    MEHPForceBalance::getEffectiveFunctionalityOfAtoms(double tolerance) const
    {
      std::unordered_map<long int, int> results;
      results.reserve(this->initialConfig.nrOfNodes);

      Eigen::VectorXi nrOfActiveSpringsConnected =
        this->getNrOfActiveSpringsConnected(tolerance);
      for (size_t i = 0; i < this->initialConfig.nrOfNodes; i++) {
        results.emplace(this->initialConfig.oldAtomIds[i],
                        nrOfActiveSpringsConnected[i]);
      }
      return results;
    }

    /**
     * @brief Get the Ids Of active Nodes
     *
     * @param tolerance the tolerance: springs under a certain length are
     * considered inactive
     * @param minimumNrOfActiveConnections the number of active springs
     * required for this node to qualify as active
     * @return std::vector<long int> the atom ids
     */
    std::vector<long int> MEHPForceBalance::getIdsOfActiveNodes(
      double tolerance,
      int minimumNrOfActiveConnections,
      int maximumNrOfActiveConnections) const
    {
      std::vector<long int> results;
      results.reserve(this->initialConfig.nrOfNodes);

      Eigen::VectorXi nrOfActiveSpringsConnected =
        this->getNrOfActiveSpringsConnected(tolerance);
      for (size_t i = 0; i < this->initialConfig.nrOfNodes; i++) {
        if (nrOfActiveSpringsConnected[i] >= minimumNrOfActiveConnections &&
            (maximumNrOfActiveConnections < 0 ||
             maximumNrOfActiveConnections >= nrOfActiveSpringsConnected[i])) {
          results.push_back(this->initialConfig.oldAtomIds[i]);
        }
      }

      return results;
    }

    /**
     * @brief Get the Nr Of Active Springs connected to each node
     *
     * @param tolerance the tolerance: springs under a certain length are
     * considered inactive
     * @return Eigen::VectorXi
     */
    Eigen::VectorXi MEHPForceBalance::getNrOfActiveSpringsConnected(
      double tolerance) const
    {
      Eigen::VectorXi nrOfActiveSpringsConnected =
        Eigen::VectorXi::Zero(this->initialConfig.nrOfNodes);
      ArrayXb springIsActive =
        this->findActiveSprings(this->currentSpringDistances, tolerance);
      for (size_t i = 0; i < this->initialConfig.nrOfSprings; i++) {
        if (springIsActive[i] == true) /* active spring */
        {
          int a = this->initialConfig.springIndexA[i];
          int b = this->initialConfig.springIndexB[i];
          ++(nrOfActiveSpringsConnected[a]);
          ++(nrOfActiveSpringsConnected[b]);
        }
      }
      return nrOfActiveSpringsConnected;
    }

    /**
     * @brief Get the Gamma Factor at the current step
     *
     * @param r02 the melt <R_0^2>, for phantom = Nb^2
     * @param nrOfChains the nr of chains to average over (can be different
     * from the nr of springs thanks to omitted free chains or primary loops)
     * @return double
     */
    double MEHPForceBalance::getGammaFactor(double r02, int nrOfChains) const
    {
      if (r02 < 0) {
        r02 = this->defaultR0Squared;
      }
      if (nrOfChains < 1) {
        nrOfChains = this->defaultNrOfChains;
      }

      return this->evaluateGammaFactor(
        this->currentSpringDistances, r02, nrOfChains);
    }

    bool MEHPForceBalance::validateNetwork(const ForceBalanceNetwork* net)
    {
      if (net == nullptr) {
        net = &this->initialConfig;
      }
      RUNTIME_EXP_IFN(!std::isinf(net->L[0]) && !std::isnan(net->L[0]),
                      "Box direction x must be scalar");
      RUNTIME_EXP_IFN(!std::isinf(net->L[1]) && !std::isnan(net->L[1]),
                      "Box direction y must be scalar");
      RUNTIME_EXP_IFN(!std::isinf(net->L[2]) && !std::isnan(net->L[2]),
                      "Box direction z must be scalar");
      RUNTIME_EXP_IFN(net->coordinates.size() == net->nrOfLinks * 3,
                      "Invalid size of coordinates");
      RUNTIME_EXP_IFN(this->currentDisplacements.size() == net->nrOfLinks * 3,
                      "Invalid size of current displacement");
      RUNTIME_EXP_IFN(net->localToGlobalSpringIndex.size() == net->nrOfSprings,
                      "Invalid size of connectivity map");
      RUNTIME_EXP_IFN(net->springsContourLength.size() == net->nrOfSprings,
                      "Invalid size of contour lengths");
      RUNTIME_EXP_IFN(net->springIndicesOfLinks.size() == net->nrOfLinks,
                      "Invalid size of spring indices of links");
      RUNTIME_EXP_IFN(net->linkIndicesOfSprings.size() == net->nrOfSprings,
                      "Invalid size of link indices of springs");
      RUNTIME_EXP_IFN(net->linkIsSliplink.size() == net->nrOfLinks,
                      "Invalid size of link is sliplink");
      RUNTIME_EXP_IFN(
        net->linkIsSliplink.count() == net->nrOfLinks - net->nrOfNodes,
        "Nr of nodes plus nr of slp-links should give the total nr of links");
      RUNTIME_EXP_IFN(net->oldAtomIds.size() == net->nrOfNodes,
                      "Invalid size of old atom ids");
      RUNTIME_EXP_IFN(net->springCoordinateIndexA.size() ==
                        net->nrOfSprings * 3,
                      "Invalid size of springCoordinateIndexA");
      RUNTIME_EXP_IFN(net->springCoordinateIndexB.size() ==
                        net->nrOfSprings * 3,
                      "Invalid size of springCoordinateIndexB");
      RUNTIME_EXP_IFN(net->springPartCoordinateIndexA.size() ==
                        net->nrOfPartialSprings * 3,
                      "Invalid size of springPartCoordinateIndexA");
      RUNTIME_EXP_IFN(net->springPartCoordinateIndexB.size() ==
                        net->nrOfPartialSprings * 3,
                      "Invalid size of springPartCoordinateIndexB");
      RUNTIME_EXP_IFN(net->springIndexA.size() == net->nrOfSprings,
                      "Invalid size of springIndexA");
      RUNTIME_EXP_IFN(net->springIndexB.size() == net->nrOfSprings,
                      "Invalid size of springIndexB");
      RUNTIME_EXP_IFN(net->springPartIndexA.size() == net->nrOfPartialSprings,
                      "Invalid size of springPartIndexA");
      RUNTIME_EXP_IFN(net->springPartIndexB.size() == net->nrOfPartialSprings,
                      "Invalid size of springPartIndexB");
      RUNTIME_EXP_IFN(net->springIsActive.size() == net->nrOfSprings,
                      "Invalid size of springIsActive");
      RUNTIME_EXP_IFN(this->currentSpringPartitionsVec.size() ==
                        net->nrOfPartialSprings,
                      "Invalid size of currentSpringPartitionsVec");
      RUNTIME_EXP_IFN(APPROX_EQUAL(this->currentSpringPartitionsVec.sum(),
                                   net->nrOfSprings,
                                   1e-3),
                      "Spring partitions should sum to 1 per spring");
      for (size_t i = 0; i < this->currentSpringPartitionsVec.size(); i++) {
        RUNTIME_EXP_IFN(this->currentSpringPartitionsVec[i] <= 1.0 &&
                          this->currentSpringPartitionsVec[i] >= 0.0,
                        "Spring partitions must be between 0 & 1, got " +
                          std::to_string(this->currentSpringPartitionsVec[i]) +
                          " at i = " + std::to_string(i) + ".");
      }
      for (size_t i = 0; i < net->nrOfSprings; ++i) {
        RUNTIME_EXP_IFN(net->linkIndicesOfSprings[i].size() >= 2,
                        "Each spring requires at least two links, got " +
                          std::to_string(net->linkIndicesOfSprings[i].size()) +
                          " at i = " + std::to_string(i) + ".");
        RUNTIME_EXP_IFN(net->localToGlobalSpringIndex.at(i).size() ==
                          net->linkIndicesOfSprings[i].size() - 1,
                        "Require a global index for each local one");
        RUNTIME_EXP_IFN(
          net->linkIndicesOfSprings[i][0] <=
            net->linkIndicesOfSprings[i]
                                     [net->linkIndicesOfSprings[i].size() - 1],
          "Springs must have increasing end-point indices");
      }
      return true;
    }
  }
}
}
