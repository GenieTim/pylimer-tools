#include "Molecule.h"
#include "../utils/GraphUtils.h"
#include "../utils/StringUtils.h"
#include "Atom.h"
#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <set>
extern "C"
{
#include <igraph/igraph.h>
}
#ifdef OPENMP_FOUND
#include <omp.h>
#endif

namespace pylimer_tools {
namespace entities {

  Molecule::Molecule(const Box& parent,
                     const igraph_t* inGraph,
                     const MoleculeType type,
                     const std::map<int, double>& massPerType)
  {
    this->parent = parent;
    this->initializeFromGraph(inGraph);
    this->typeOfThisMolecule = type;
    this->massPerType = massPerType;
  };

  /**
   * @brief Utility method to initialise the molecule, reading all info from the
   * graph
   *
   * @param inGraph
   */
  void Molecule::initializeFromGraph(const igraph_t* inGraph)
  {
    igraph_copy(&this->graph, inGraph);
    this->size = igraph_vcount(&this->graph);
    // construct a key for this molecule: a concatenation of all ids in this
    // molecule
    if (!igraph_cattribute_has_attr(
          &this->graph, IGRAPH_ATTRIBUTE_VERTEX, "id")) {
      throw std::runtime_error("Molecule's graph does not have attribute id");
    }
    igraph_vector_t allIds;
    igraph_vector_init(&allIds, this->size);
    VANV(&this->graph, "id", &allIds);
    if (igraph_cattribute_VANV(&this->graph, "id", igraph_vss_all(), &allIds)) {
      throw std::runtime_error(
        "Molecule's graph's attribute id is not accessible.");
    };
    std::vector<int> ids;
    pylimer_tools::utils::igraphVectorTToStdVector(&allIds, ids);
    if (ids.size() == 0 && this->size > 0) {
      throw std::runtime_error(
        "Molecule's graph's attribute id was not queried.");
    }
    this->atomIdToVertexIdx.reserve(ids.size());
    for (int i = 0; i < ids.size(); ++i) {
      this->atomIdToVertexIdx[ids[i]] = i;
    }
    // TODO: possibly use the lined up atoms instead? Check performance.
    std::sort(ids.begin(), ids.end());
    this->key =
      pylimer_tools::utils::join(ids.begin(), ids.end(), std::string("-"));
    igraph_vector_destroy(&allIds);
    // detect whether the graph has more than the standard attributes
    this->atomsHaveCustomAttributes = this->checkIfAtomsHaveCustomAttributes();
  };

  // rule of three:
  // 1. destructor (to destroy the graph)
  Molecule::~Molecule()
  {
    // in addition to basic fields being deleted, we need to clean up the graph
    // as is done in parent
    igraph_destroy(&this->graph);
  };
  // 2. copy constructor
  Molecule::Molecule(const Molecule& src)
    : Molecule(src.parent,
               &src.graph,
               src.typeOfThisMolecule,
               src.massPerType) {};
  // 3. copy assignment operator
  Molecule& Molecule::operator=(Molecule src)
  {
    std::swap(this->parent, src.parent);
    std::swap(this->typeOfThisMolecule, src.typeOfThisMolecule);
    std::swap(this->atomsHaveCustomAttributes, src.atomsHaveCustomAttributes);
    std::swap(this->size, src.size);
    std::swap(this->key, src.key);
    std::swap(this->massPerType, src.massPerType);
    std::swap(this->graph, src.graph);

    return *this;
  };

  // other operators
  bool Molecule::operator==(const Molecule& ref) const
  {
    igraph_t difference;
    igraph_difference(&difference, &this->graph, &ref.graph);
    igraph_integer_t ecount = igraph_ecount(&difference);
    igraph_destroy(&difference);
    if (ecount != 0) {
      return false;
    }
    std::vector<Atom> thisAtoms = this->getAtoms();
    std::vector<Atom> otherAtoms = ref.getAtoms();
    if (thisAtoms.size() != otherAtoms.size()) {
      return false;
    }
    for (size_t i = 0; i < otherAtoms.size(); ++i) {
      if (thisAtoms[i] != otherAtoms[i]) {
        return false;
      }
    }
    return this->massPerType == ref.massPerType &&
           this->typeOfThisMolecule == ref.typeOfThisMolecule &&
           this->parent == ref.parent;
  }

  std::vector<double> Molecule::computeBondLengths()
  {
    return AtomGraphParent::computeBondLengths(this->parent);
  }

  double Molecule::computeTotalLength()
  {
    std::vector<double> bondLens = this->computeBondLengths();
    return std::accumulate(bondLens.begin(), bondLens.end(), 0.0);
  }

  double Molecule::computeEndToEndDistance() const
  {
    return this->computeEndToEndVector().norm();
  }

  Eigen::Vector3d Molecule::computeEndToEndVector() const
  {
    if (this->getNrOfAtoms() < 2) {
      return Eigen::Vector3d::Zero();
    }

    std::vector<Atom> endNodes = this->getChainEnds();

    // we only compute an end-to-end distance if we have exactly two ends.
    // this is clearly not optimal, but at least unambiguous
    if (endNodes.size() == 2) {
      // TODO: this is more intensive than needed
      // check whether the compiler optimizes this or not
      Atom atom1 = endNodes[0];
      Atom atom2 = endNodes[1];
      return atom1.vectorToUnwrapped(atom2, this->parent);
    }

    throw std::runtime_error(
      "Cannot compute end-to-end vector for Molecule with " +
      std::to_string(endNodes.size()) + " end(s) in molecule " +
      this->getKey() + ".");
  }

  double Molecule::computeEndToEndDistanceWithDerivedImageFlags() const
  {
    return this->computeEndToEndVectorWithDerivedImageFlags().norm();
  }

  Eigen::Vector3d Molecule::computeEndToEndVectorWithDerivedImageFlags() const
  {
    const std::vector<long int> vertices = this->getVerticesLinedUp();
    if (vertices.size() == 0 || vertices.size() == 1) {
      return Eigen::Vector3d::Zero();
    }

    Eigen::VectorXd coordinates = Eigen::VectorXd(vertices.size() * 3);
    this->getAssumedVertexCoordinates(coordinates, this->parent, vertices);

    Eigen::Vector3d distance =
      coordinates.segment(3 * (vertices.size() - 1), 3) -
      coordinates.segment(0, 3);

    return distance;
  }

  /**
   * @brief compute the weight of this molecule
   *
   * @return double the total weight
   */
  double Molecule::computeTotalMass()
  {
    std::vector<int> presentTypes = this->getPropertyValues<int>("type");
    double totalWeight =
      std::accumulate(presentTypes.begin(),
                      presentTypes.end(),
                      0.0,
                      [&massPerType = this->massPerType](double val, int type) {
                        return val + massPerType[type];
                      });
    return totalWeight;
  }

  long int Molecule::getAtomIdByIdx(const int vertexId) const
  {
    return VAN(&this->graph, "id", vertexId);
  };

  long int Molecule::getIdxByAtomId(const int atomId) const
  {
    if (!pylimer_tools::utils::map_has_key(this->atomIdToVertexIdx, atomId)) {
      throw std::invalid_argument("Molecule cannot return vertex idx of this "
                                  "atom: an atom with this id (" +
                                  std::to_string(atomId) + ") does not exist");
    }
    return this->atomIdToVertexIdx.at(atomId);
  };

  bool Molecule::containsAtom(const Atom& atom) const
  {
    if (!pylimer_tools::utils::map_has_key(this->atomIdToVertexIdx,
                                           atom.getId())) {
      return false;
    }
    return this->getAtomByVertexIdx(this->getIdxByAtomId(atom.getId())) == atom;
  }

  /**
   * @brief Get the nr of atoms in the molecule
   *
   * @return int
   */
  int Molecule::getLength() const
  {
    return this->size;
  };

  /**
   * @brief Get the nr of atoms in the molecule
   *
   * @return int
   */
  int Molecule::getNrOfAtoms() const
  {
    return this->size;
    // return this->getNrOfVertices();
  }

  int Molecule::getNrOfBonds() const
  {
    return this->getNrOfEdges();
  }

  /**
   * @brief Get the type of the molecule
   *
   * @return MoleculeType
   */
  MoleculeType Molecule::getType() const
  {
    return this->typeOfThisMolecule;
  };

  const Box& Molecule::getBox() const
  {
    return this->parent;
  }

  double Molecule::computeRadiusOfGyration()
  {
    double meanX = 0.0, meanY = 0.0, meanZ = 0.0;
    // would be faster to just query the attributes.
    // But the OOP interface is just too tempting
    // as long as there are no external additional performance demands
    std::vector<Atom> allAtoms = this->getAtoms();
    if (allAtoms.size() == 0) {
      return 0.0;
    }
    double multiplier = 1. / static_cast<double>(allAtoms.size());
    double totalMass = 0.0;

    if (this->massPerType.empty()) {
      throw std::runtime_error(
        "Cannot compute radius of gyration without masses.");
    }

    // TODO: might want to use the raw values, use std::accumulate or
    // std::reduce
#pragma omp parallel for reduction(+ : meanX, meanY, meanZ)
    for (Atom a : allAtoms) {
      meanX +=
        this->massPerType.at(a.getType()) * a.getUnwrappedX(this->parent);
      // meanNx += a.getNX();
      meanY +=
        this->massPerType.at(a.getType()) * a.getUnwrappedY(this->parent);
      // meanNy += a.getNY();
      meanZ +=
        this->massPerType.at(a.getType()) * a.getUnwrappedZ(this->parent);
      // meanNz += a.getNZ();
      totalMass += this->massPerType.at(a.getType());
    }

    Atom virtualCenterAtom = Atom(0,
                                  0,
                                  meanX * multiplier,
                                  meanY * multiplier,
                                  meanZ * multiplier,
                                  0,
                                  0,
                                  0);

    double correctingFactor = 1. / totalMass;

    // reduce to the mean
    auto innerReduction = [&virtualCenterAtom,
                           correctingFactor,
                           &massPerType = this->massPerType,
                           box = this->parent](double val,
                                               const Atom a) -> double {
      double dist = a.distanceToUnwrapped(virtualCenterAtom, box);
      return val +
             (correctingFactor * massPerType.at(a.getType()) * dist * dist);
    };
    double Rg2 =
      std::accumulate(allAtoms.begin(), allAtoms.end(), 0., innerReduction);

    return Rg2;
  }

  double Molecule::computeRadiusOfGyrationWithDerivedImageFlags() const
  {
    const std::vector<long int> vertices = this->getVerticesLinedUp();
    if (vertices.size() == 0 || vertices.size() == 1) {
      return 0.0;
    }

    if (this->massPerType.empty()) {
      throw std::runtime_error(
        "Cannot compute radius of gyration without masses.");
    }

    double multiplier = 1. / (static_cast<double>(vertices.size()));

    // compute the mean position based on the
    // image flags of the first atom
    // Atom lastAtom = this->getAtomByVertexIdx(vertices[0]);
    std::vector<int> atomTypes = this->getPropertyValues<int>("type", vertices);
    Eigen::VectorXd assumedCoordinates =
      Eigen::VectorXd::Zero(vertices.size() * 3);
    this->getAssumedVertexCoordinates<Eigen::VectorXd>(
      assumedCoordinates, this->parent, vertices);
    double totalMass = 0.0;
    Eigen::Vector3d meanCoords = Eigen::Vector3d::Zero();
    // find mean position
    for (size_t i = 0; i < vertices.size(); ++i) {
      double localMultiplier = multiplier * this->massPerType.at(atomTypes[i]);

      totalMass += this->massPerType.at(atomTypes[i]);
      meanCoords += localMultiplier * assumedCoordinates.segment(3 * i, 3);
    }

    // use it to compute the r_g
    double Rg2 = 0.0;

    multiplier = 1. / totalMass;
    for (size_t i = 0; i < vertices.size(); ++i) {
      double localMultiplier = multiplier * this->massPerType.at(atomTypes[i]);
      Eigen::Vector3d distanceFromMean =
        assumedCoordinates.segment(3 * i, 3) - meanCoords;
      Rg2 += localMultiplier * distanceFromMean.squaredNorm();
    }

    return Rg2;
  };

  std::string Molecule::getKey() const
  {
    return this->key;
  }

  std::vector<Atom> Molecule::getAtoms() const
  {
    std::vector<Atom> results;
    size_t nrOfAtoms = this->getNrOfAtoms();
    results.reserve(nrOfAtoms);

    // #pragma omp declare reduction (merge : std::vector<Atom> :
    // omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
    // #pragma omp parallel for reduction(merge: results)
    for (size_t i = 0; i < nrOfAtoms; ++i) {
      results.push_back(this->getAtomByVertexIdx(i));
    }

    return results;
  };

  /**
   * @brief Get the atoms at the end of each chain, and/or, in case of a primary
   * loop, two times the cross-link.
   *
   * This method does not handle branched structures.
   *
   * @param crossLinkerType
   * @return std::vector<Atom, Atom>
   */
  std::vector<Atom> Molecule::getChainEnds(const int crossLinkerType,
                                           bool closePrimaryLoop) const
  {
    if (this->getNrOfAtoms() == 0) {
      return {};
    }
    if (this->getNrOfAtoms() == 1) {
      if (this->getNrOfBonds() > 0 && closePrimaryLoop) {
        return { this->getAtomByVertexIdx(0), this->getAtomByVertexIdx(0) };
      } else {
        return { this->getAtomByVertexIdx(0) };
      }
    }
    std::vector<pylimer_tools::entities::Atom> endsOfChain =
      this->getAtomsOfDegree(1);

    if (this->getType() == MoleculeType::PRIMARY_LOOP ||
        ((endsOfChain.size() == 1 || endsOfChain.size() == 0) &&
         this->getType() == MoleculeType::UNDEFINED)) {
      std::vector<pylimer_tools::entities::Atom> crossLinks =
        this->getAtomsOfType(crossLinkerType);
      if (crossLinks.size() == 0) {
        // can happen e.g. in "free" primary loops.
        // we could alternatively select a random atom,
        // e.g. the one with the lowest id
        return {};
      }
      // otherwise, primary loop,
      // but with cross-link, which will be the start of the loop.
      RUNTIME_EXP_IFN(
        crossLinks.size() == 1 ||
          (crossLinks.size() == 2 &&
           crossLinks[0].getId() == crossLinks[1].getId()),
        "For a primary loop, expected exactly one cross-link, got " +
          std::to_string(crossLinks.size()) + " in molecule " + this->getKey() +
          ".");
      return closePrimaryLoop
               ? std::vector<Atom>({ crossLinks[0], crossLinks[0] })
               : std::vector<Atom>({ crossLinks[0] });
    }
    RUNTIME_EXP_IFN(
      endsOfChain.size() == 2,
      "Expected to find the two ends of the chain of type " +
        std::to_string(this->getType()) + " as atoms of degree 1, but found " +
        std::to_string(endsOfChain.size()) +
        " atoms with degree 1 in molecule " + this->getKey() + ".");
    if (endsOfChain[0].getId() > endsOfChain[1].getId()) {
      return { endsOfChain[1], endsOfChain[0] };
    }
    return { endsOfChain[0], endsOfChain[1] };
  }

  /**
   * @brief Get the overall offset in terms of boxes (for PBC)
   *
   * The offset is computed as if computing the vector of the first to the last
   * atom (coords of last minus coords of first).
   *
   * NOTE: even for primary loops, it is possible that this is not equal to
   * zero.
   * @param crossLinkerType
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d Molecule::getOverallBondSum(const int crossLinkerType,
                                              const bool closeLoop) const
  {
    std::vector<long int> alignedVertices =
      this->getVerticesLinedUp(crossLinkerType, closeLoop);
    Eigen::VectorXd alignedCoordinates =
      Eigen::VectorXd::Zero(3 * alignedVertices.size());
    this->getAssumedVertexCoordinates(
      alignedCoordinates, this->parent, alignedVertices);
    Eigen::Vector3d result = Eigen::Vector3d::Zero();
    for (size_t i = 1; i < alignedVertices.size(); ++i) {
      Eigen::Vector3d distance = alignedCoordinates.segment((i) * 3, 3) -
                                 alignedCoordinates.segment((i - 1) * 3, 3);
      this->parent.handlePBC(distance);
      result += distance;
    }
    return result;
  };

  /**
   * @brief Get the overall offset in terms of boxes (for PBC)
   *
   * The offset is computed as if computing the vector of the first to the last
   * atom (coords of last minus coords of first).
   *
   * NOTE: even for primary loops, it is possible that this is not equal to
   * zero.
   * @param atomIdFrom
   * @param atomIdTo
   * @param crossLinkerType
   * @param requireOrder
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d Molecule::getOverallBondSumFromTo(size_t atomIdFrom,
                                                    size_t atomIdTo,
                                                    const int crossLinkerType,
                                                    bool requireOrder) const
  {
    std::vector<long int> alignedVertices =
      this->getVerticesLinedUp(crossLinkerType, true);
    Eigen::VectorXd alignedCoordinates =
      Eigen::VectorXd::Zero(3 * alignedVertices.size());
    this->getAssumedVertexCoordinates(
      alignedCoordinates, this->parent, alignedVertices);
    Eigen::Vector3d result = Eigen::Vector3d::Zero();
    size_t vertexIdFrom = this->atomIdToVertexIdx.at(atomIdFrom);
    size_t vertexIdTo = this->atomIdToVertexIdx.at(atomIdTo);
    bool recording = false;
    for (size_t i = 0; i < alignedVertices.size(); ++i) {
      if (recording) {
        Eigen::Vector3d distance = alignedCoordinates.segment((i) * 3, 3) -
                                   alignedCoordinates.segment((i - 1) * 3, 3);
        this->parent.handlePBC(distance);
        result += distance;

        if ((alignedVertices[i] == vertexIdFrom && !requireOrder) ||
            alignedVertices[i] == vertexIdTo) {
          return result * (alignedVertices[i] == vertexIdFrom ? -1. : 1.);
        }
      }

      if (alignedVertices[i] == vertexIdFrom ||
          (alignedVertices[i] == vertexIdTo && !requireOrder)) {
        recording = true;
      }
    }
    throw std::runtime_error(
      "Did not find both vertices to compute overall bond sum for");
  };

  size_t Molecule::getNrOfBondsFromTo(size_t atomIdFrom,
                                      size_t atomIdTo,
                                      const int crossLinkerType,
                                      bool requireOrder) const
  {
    std::vector<long int> alignedVertices =
      this->getVerticesLinedUp(crossLinkerType, true);
    size_t vertexIdFrom = this->atomIdToVertexIdx.at(atomIdFrom);
    size_t vertexIdTo = this->atomIdToVertexIdx.at(atomIdTo);
    bool recording = false;
    size_t result = 0;
    for (size_t i = 0; i < alignedVertices.size(); ++i) {
      if (recording) {
        result += 1;

        if ((alignedVertices[i] == vertexIdFrom && !requireOrder) ||
            alignedVertices[i] == vertexIdTo) {
          return result;
        }
      }

      if (alignedVertices[i] == vertexIdFrom ||
          (alignedVertices[i] == vertexIdTo && !requireOrder)) {
        recording = true;
      }
    }
    throw std::runtime_error(
      "Did not find both vertices (" + std::to_string(atomIdTo) + " & " +
      std::to_string(atomIdFrom) + ") to count atoms between");
  }

  /**
   * @brief Get the ids of the vertices in order of the chain, starting from one
   * end to the other
   *
   * @param crossLinkerType
   * @return std::vector<long int>
   */
  std::vector<long int> Molecule::getVerticesLinedUp(int crossLinkerType,
                                                     bool closeLoop) const
  {
    std::vector<long int> results;
    size_t nrOfAtoms = this->getNrOfAtoms();
    results.reserve(nrOfAtoms);

    if (nrOfAtoms == 1) {
      results.push_back(0);
      return results;
    }
    if (nrOfAtoms == 0) {
      return results;
    }

    std::vector<Atom> chainEnds = this->getChainEnds(crossLinkerType);
    long int vertexIdToStartWith = 0;
    // chainEnds may be empty e.g. for "free" primary loops
    if (chainEnds.size() > 0) {
      vertexIdToStartWith = this->getIdxByAtomId(chainEnds[0].getId());
    }

    std::vector<long int> connections =
      this->getVertexIdxsConnectedTo(vertexIdToStartWith);
    results.push_back(vertexIdToStartWith);
    bool loopFound = false;
    for (long int connection : connections) {
      long int currentCenter = connection;
      results.push_back(currentCenter);
      long int lastCenter = vertexIdToStartWith;
      std::vector<long int> subConnections =
        this->getVertexIdxsConnectedTo(currentCenter);
      while (subConnections.size() > 0) {
        if (subConnections.size() == 1) {
          break;
        }
        // we assume a functionality of 2 for ordinary strands
        if (subConnections.size() != 2) {
          throw std::runtime_error(
            "Failed to align all atoms on one strand, as a functionality of " +
            std::to_string(subConnections.size()) +
            " was found and 1 or 2 expected.");
        }
        int subConnectionDirection = (subConnections[0] == lastCenter) ? 1 : 0;
        if (subConnections[subConnectionDirection] == vertexIdToStartWith) {
          loopFound = true;
          break;
        }
        lastCenter = currentCenter;
        currentCenter = subConnections[subConnectionDirection];
        results.push_back(currentCenter);
        subConnections = this->getVertexIdxsConnectedTo(currentCenter);
      }
      if (loopFound) {
        assert(this->getType() == MoleculeType::UNDEFINED ||
               this->getType() == MoleculeType::PRIMARY_LOOP);
        break;
      }
    }

    if (results.size() != this->getNrOfAtoms()) {
      throw std::runtime_error(
        "Failed to align all atoms on one strand: Lined up " +
        std::to_string(results.size()) + " instead of " +
        std::to_string(this->getNrOfAtoms()) + " atoms.");
    }

    if (closeLoop && loopFound) {
      results.push_back(results[0]);
    }

    return results;
  }

  std::vector<Atom> Molecule::getAtomsLinedUp(int crossLinkerType,
                                              bool assumedCoordinates,
                                              bool closeLoop) const
  {
    std::vector<long int> vertices =
      this->getVerticesLinedUp(crossLinkerType, closeLoop);
    if (!assumedCoordinates) {
      return this->verticesToAtoms(vertices);
    } else {
      std::vector<pylimer_tools::entities::Atom> results;
      results.reserve(vertices.size());
      std::vector<double> coordinates =
        pylimer_tools::utils::initializeWithValue(vertices.size() * 3, 0.0);
      this->getAssumedVertexCoordinates<std::vector<double>>(
        coordinates, this->parent, vertices);
      for (size_t i = 0; i < vertices.size(); ++i) {
        long int vertex = vertices[i];
        results.push_back(pylimer_tools::entities::Atom(
          this->getAtomIdByIdx(vertex),
          this->getPropertyValue<int>("type", vertex),
          coordinates[3 * i],
          coordinates[3 * i + 1],
          coordinates[3 * i + 2],
          0,
          0,
          0));
      }
      return results;
    }
  };

} // namespace entities
} // namespace pylimer_tools
