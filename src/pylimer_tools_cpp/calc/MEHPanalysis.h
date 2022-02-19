#ifndef MEHP_ANALYSIS_H
#define MEHP_ANALYSIS_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Molecule.h"
#include "../entities/Universe.h"
#include "../entities/UniverseSequence.h"
#include <algorithm>
#include <array>
#include <map>
#include <string>

namespace pylimer_tools {
namespace calc {
namespace mehp {

typedef std::array<double, 3> position_vec_t;

/*
Compute the end to end vectors between each pair of (indirectly) connected
crosslinker

Arguments:
  - network: the polymer network to do the computation for
  - crosslinkerType: the atom type to compute the in-between vectors for

Returns:
  - endToEndVectors (map): a map with key: "{molecule.key}"
      and value: their difference vector
*/
std::map<std::string, position_vec_t>
computeEndToEndVectors(pylimer_tools::entities::Universe network,
                       int crosslinkerType) {
  std::map<std::string, position_vec_t> result;
  std::vector<pylimer_tools::entities::Molecule> molecules =
      network.getChainsWithCrosslinker(crosslinkerType);
  pylimer_tools::entities::Box box = network.getBox();

  for (pylimer_tools::entities::Molecule chain : molecules) {
    std::vector<pylimer_tools::entities::Atom> crosslinkers =
        chain.getAtomsWithType(crosslinkerType);

    if (crosslinkers.size() != 2 ||
        chain.getType() ==
            pylimer_tools::entities::MoleculeType::PRIMARY_LOOP ||
        chain.getType() ==
            pylimer_tools::entities::MoleculeType::DANGLING_CHAIN) {
      continue;
    }

    // use id to keep direction of the vector constant
    double distanceVec[3];
    if (crosslinkers[0].getId() > crosslinkers[1].getId()) {
      crosslinkers[0].vectorTo(crosslinkers[1], &box, distanceVec);
    } else {
      crosslinkers[1].vectorTo(crosslinkers[0], &box, distanceVec);
    }
    position_vec_t distanceVecT;
    std::copy_n(std::begin(distanceVec), 3, std::begin(distanceVecT));
    result.insert_or_assign(chain.getKey(), distanceVecT);
  }

  return result;
}

/*
Compute the mean end to end vectors between each pair of (indirectly) connected
crosslinker

Arguments:
  - networks: the different configurations of the polymer network to do the
computation for
  - crosslinkerType: the atom type to compute the in-between vectors for

Returns:
  - endToEndVectors (map): a dictionary with key: "{chain.key}"
      and value: their mean distance difference vector
*/
std::map<std::string, position_vec_t>
computeMeanEndToEndVectors(pylimer_tools::entities::UniverseSequence networks,
                           int crosslinkerType) {
  std::map<std::string, position_vec_t> result;

  if (networks.getLength() == 0) {
    return result;
  }

  double multiplier = 1.0 / networks.getLength();

  for (int i = 0; i < networks.getLength(); ++i) {
    pylimer_tools::entities::Universe network = networks.atIndex(i);
    std::map<std::string, position_vec_t> currentEndToEndVectors =
        computeEndToEndVectors(network, crosslinkerType);

    for (auto const &[key, vec] : currentEndToEndVectors) {
      if (!result.contains(key)) {
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
Compute the mean end to end distance between each pair of (indirectly) connected
crosslinker

Arguments:
- networks: the different configurations of the polymer network to do the
computation for
- crosslinkerType: the atom type to compute the in-between vectors for

Returns:
- endToEndDistances (dict): a dictionary with key: "{atom1.name}+{atom2.name}"
    and value: the norm of the mean difference vector
*/
std::map<std::string, double>
computeMeanEndToEndDistances(pylimer_tools::entities::UniverseSequence networks,
                             int crosslinkerType) {
  std::map<std::string, position_vec_t> distanceVectors =
      computeMeanEndToEndVectors(networks, crosslinkerType);
  std::map<std::string, double> results;

  for (auto const &[key, vec] : distanceVectors) {
    results.insert_or_assign(
        key, sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]));
  }

  return results;
}

} // namespace mehp
} // namespace calc
} // namespace pylimer_tools

#endif
