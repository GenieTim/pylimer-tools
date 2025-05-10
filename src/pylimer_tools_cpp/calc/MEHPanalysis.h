#ifndef MEHP_ANALYSIS_H
#define MEHP_ANALYSIS_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Molecule.h"
#include "../entities/Universe.h"
#include "../entities/UniverseSequence.h"
#include <array>
#include <map>
#include <string>

namespace pylimer_tools::calc::mehp {

typedef std::array<double, 3> position_vec_t;

/*
Compute the end to end vectors between each pair of (indirectly) connected
crossLinker

Arguments:
  - network: The polymer network to do the computation for
  - crossLinkerType: The atom type to compute the in-between vectors for

Returns:
  - endToEndVectors (map): a map with key: "{molecule.key}"
          and value: Their difference vector
*/
std::map<std::string, position_vec_t>
computeEndToEndVectors(pylimer_tools::entities::Universe network,
                       const int crossLinkerType)
{
  std::map<std::string, position_vec_t> result;
  std::vector<pylimer_tools::entities::Molecule> molecules =
    network.getChainsWithCrosslinker(crossLinkerType);
  pylimer_tools::entities::Box box = network.getBox();

  for (pylimer_tools::entities::Molecule chain : molecules) {
    std::vector<pylimer_tools::entities::Atom> crossLinkers =
      chain.getAtomsOfType(crossLinkerType);

    if (crossLinkers.size() != 2 ||
        chain.getType() ==
          pylimer_tools::entities::MoleculeType::PRIMARY_LOOP ||
        chain.getType() ==
          pylimer_tools::entities::MoleculeType::DANGLING_CHAIN) {
      continue;
    }

    // use id to keep direction of the vector constant
    Eigen::Vector3d distanceVec;
    if (crossLinkers[0].getId() > crossLinkers[1].getId()) {
      distanceVec = crossLinkers[0].vectorTo(crossLinkers[1], box);
    } else {
      distanceVec = crossLinkers[1].vectorTo(crossLinkers[0], box);
    }
    position_vec_t distanceVecT;
    std::copy_n(std::begin(distanceVec), 3, std::begin(distanceVecT));
    result.insert_or_assign(chain.getKey(), distanceVecT);
  }

  return result;
}

/*
Compute the mean end to end vectors between each pair of (indirectly)
connected crossLinker

Arguments:
  - networks: The different configurations of the polymer network to do the
computation for
  - crossLinkerType: The atom type to compute the in-between vectors for

Returns:
  - endToEndVectors (map): a dictionary with key: "{chain.key}"
          and value: Their mean distance difference vector
*/
std::map<std::string, position_vec_t>
computeMeanEndToEndVectors(pylimer_tools::entities::UniverseSequence networks,
                           const int crossLinkerType)
{
  std::map<std::string, position_vec_t> result;

  if (networks.getLength() == 0) {
    return result;
  }

  double multiplier = 1.0 / networks.getLength();

  for (int i = 0; i < networks.getLength(); ++i) {
    pylimer_tools::entities::Universe network = networks.atIndex(i);
    std::map<std::string, position_vec_t> currentEndToEndVectors =
      computeEndToEndVectors(network, crossLinkerType);

    for (auto const& [key, vec] : currentEndToEndVectors) {
      if (!pylimer_tools::utils::map_has_key(result, key)) {
        position_vec_t zeroPosition;
        zeroPosition.fill(0.0);
        result.insert_or_assign(key, zeroPosition);
      }
      for (int j = 0; j < 3; j++) {
        result[key][j] += vec[j] * multiplier;
      }
    }
  }

  return result;
}

/*
Compute the mean end to end distance between each pair of (indirectly)
connected crossLinker

Arguments:
- networks: The different configurations of the polymer network to do the
computation for
- crossLinkerType: The atom type to compute the in-between vectors for

Returns:
- endToEndDistances (dict): a dictionary with key:
"{atom1.name}+{atom2.name}" and value: The norm of the mean difference
vector
*/
std::map<std::string, double>
computeMeanEndToEndDistances(pylimer_tools::entities::UniverseSequence networks,
                             const int crossLinkerType)
{
  std::map<std::string, position_vec_t> distanceVectors =
    computeMeanEndToEndVectors(networks, crossLinkerType);
  std::map<std::string, double> results;

  for (auto const& [key, vec] : distanceVectors) {
    results.insert_or_assign(
      key, sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]));
  }

  return results;
}

}

#endif
