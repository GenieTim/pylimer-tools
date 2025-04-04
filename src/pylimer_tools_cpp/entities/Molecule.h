#pragma once

extern "C"
{
#include <igraph.h>
}
#include "Atom.h"
#include "AtomGraphParent.h"
#include "Box.h"
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace pylimer_tools {
namespace entities {

  enum MoleculeType
  {
    UNDEFINED,
    NETWORK_STRAND,
    PRIMARY_LOOP,
    DANGLING_CHAIN,
    FREE_CHAIN
  };

  class Molecule : public AtomGraphParent
  {
  public:
    Molecule(const Box& parent,
             const igraph_t* graph,
             const MoleculeType type,
             const std::map<int, double>& massPerType);

    // rule of three:
    // 1. destructor (to destroy the graph)
    ~Molecule();
    // 2. copy constructor
    Molecule(const Molecule& src);
    // 3. copy assignment operator
    Molecule& operator=(Molecule src);

    // other operators
    bool operator==(const Molecule& ref) const;

    // getters
    int getLength() const;
    MoleculeType getType() const;
    std::vector<Atom> getAtoms() const;
    // std::map<std::string, std::vector<long int>> getBonds() const;
    std::vector<Atom> getAtomsLinedUp(int crossLinkerType = 2,
                                      bool assumedCoordinates = false,
                                      bool closeLoop = false) const;
    std::vector<long int> getVerticesLinedUp(int crossLinkerType = 2,
                                             bool closeLoop = false) const;
    int getNrOfAtoms() const;
    int getNrOfBonds() const;
    const Box& getBox() const;
    std::string getKey() const;
    std::vector<int> getAtomTypes()
    {
      return this->getPropertyValues<int>("type");
    }
    long int getAtomIdByIdx(const int vertexId) const override;
    long int getIdxByAtomId(const int atomId) const override;
    bool containsAtom(const Atom& atom) const;
    std::vector<Atom> getChainEnds(int crossLinkerType = 2,
                                   bool closePrimaryLoop = true) const;

    // computations
    Eigen::Vector3d computeEndToEndVector() const;
    double computeEndToEndDistance() const;
    Eigen::Vector3d computeEndToEndVectorWithDerivedImageFlags() const;
    double computeEndToEndDistanceWithDerivedImageFlags() const;
    double computeRadiusOfGyration();
    double computeRadiusOfGyrationWithDerivedImageFlags() const;
    double computeTotalMass();
    std::vector<double> computeBondLengths();
    double computeTotalLength();

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
    Atom operator[](size_t index) const
    {
      return this->getAtomByVertexIdx(index);
    }

  private:
    Box parent;
    MoleculeType typeOfThisMolecule;
    int size;
    std::string key;
    std::map<int, double> massPerType;
    std::unordered_map<long int, long int> atomIdToVertexIdx;

    void initializeFromGraph(const igraph_t* ingraph);
  };
} // namespace entities
} // namespace pylimer_tools
