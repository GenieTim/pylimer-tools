#pragma once

extern "C"
{
#include <igraph/igraph.h>
}
#include "Atom.h"
#include "AtomGraphParent.h"
#include "Box.h"
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace pylimer_tools::entities {

enum MoleculeType
{
  UNDEFINED,
  NETWORK_STRAND,
  PRIMARY_LOOP,
  DANGLING_CHAIN,
  FREE_CHAIN
};

class Molecule final : public AtomGraphParent
{
public:
  Molecule(const Box& parent,
           const igraph_t* graph,
           MoleculeType type,
           const std::map<int, double>& massPerType);

  // rule of three:
  // 1. destructor (to destroy the graph)
  ~Molecule() override;
  // 2. copy constructor
  Molecule(const Molecule& src);
  // 3. copy assignment operator
  Molecule& operator=(Molecule src);

  // other operators
  bool operator==(const Molecule& ref) const;

  // getters
  [[nodiscard]] int getLength() const;
  [[nodiscard]] MoleculeType getType() const;
  [[nodiscard]] std::vector<Atom> getAtoms() const;
  // std::map<std::string, std::vector<long int>> getBonds() const;
  [[nodiscard]] std::vector<Atom> getAtomsLinedUp(
    int crossLinkerType = 2,
    bool assumedCoordinates = false,
    bool closeLoop = false) const;
  [[nodiscard]] std::vector<igraph_integer_t> getVerticesLinedUp(
    int crossLinkerType = 2,
    bool closeLoop = false) const;
  [[nodiscard]] int getNrOfAtoms() const;
  [[nodiscard]] int getNrOfBonds() const;
  [[nodiscard]] const Box& getBox() const;
  [[nodiscard]] std::string getKey() const;
  [[nodiscard]] std::vector<int> getAtomTypes() const
  {
    return this->getPropertyValues<int>("type");
  }
  [[nodiscard]] long int getAtomIdByIdx(
    igraph_integer_t vertexId) const override;
  [[nodiscard]] igraph_integer_t getIdxByAtomId(long int atomId) const override;
  [[nodiscard]] bool containsAtom(const Atom& atom) const;
  [[nodiscard]] std::vector<Atom> getChainEnds(
    int crossLinkerType = 2,
    bool closePrimaryLoop = true) const;

  // computations
  [[nodiscard]] Eigen::Vector3d computeEndToEndVector() const;
  [[nodiscard]] double computeEndToEndDistance() const;
  [[nodiscard]] Eigen::Vector3d computeEndToEndVectorWithDerivedImageFlags()
    const;
  [[nodiscard]] double computeEndToEndDistanceWithDerivedImageFlags() const;
  [[nodiscard]] double computeRadiusOfGyration();
  [[nodiscard]] double computeRadiusOfGyrationWithDerivedImageFlags() const;
  [[nodiscard]] double computeTotalMass();
  [[nodiscard]] std::vector<double> computeBondLengths() const;
  [[nodiscard]] double computeTotalLength();

  /**
   * @brief Get the sum of all bond vectors, similar to
   * `computeEndToEndDistanceWithDerivedImageFlags`
   *
   * The offset is computed as if computing the vector of the first to the
   * last atom (coords of last minus coords of first).
   *
   * NOTE: even for primary loops, it is possible that this is not equal to
   * zero.
   * @param crossLinkerType
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d getOverallBondSum(const int crossLinkerType = 2,
                                    const bool closeLoop = true) const;

  /**
   * @brief Get the overall offset in terms of boxes (for PBC)
   *
   * The offset is computed as if computing the vector of the first to the
   * last atom (coords of last minus coords of first).
   *
   * NOTE: even for primary loops, it is possible that this is not equal to
   * zero.
   * @param atomIdFrom
   * @param atomIdTo
   * @param crossLinkerType
   * @param requireOrder whether to throw an error if atomIdTo is occurring
   * before atomIdFrom
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d getOverallBondSumFromTo(size_t atomIdFrom,
                                          size_t atomIdTo,
                                          const int crossLinkerType = 2,
                                          bool requireOrder = true) const;

  size_t getNrOfBondsFromTo(size_t atomIdFrom,
                            size_t atomIdTo,
                            const int crossLinkerType = 2,
                            bool requireOrder = true) const;

  // operators
  Atom operator[](const size_t index) const
  {
    return this->getAtomByVertexIdx(static_cast<igraph_integer_t>(index));
  }

private:
  Box parent;
  MoleculeType typeOfThisMolecule;
  int size = 0;
  std::string key;
  std::map<int, double> massPerType;
  std::unordered_map<long int, long int> atomIdToVertexIdx;

  void initializeFromGraph(const igraph_t* ingraph);
};
} // namespace pylimer_tools::entities
