#ifndef MOLECULE_H
#define MOLECULE_H

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
    Molecule(const Box* parent,
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
    // getters
    int getLength() const;
    MoleculeType getType() const;
    std::vector<Atom> getAtoms();
    // std::map<std::string, std::vector<long int>> getBonds() const;
    std::vector<Atom> getAtomsLinedUp(int crossLinkerType = 2,
                                      bool assumedCoordinates = false,
                                      bool closeLoop = false) const;
    std::vector<long int> getVerticesLinedUp(int crossLinkerType = 2,
                                             bool closeLoop = false) const;
    int getNrOfAtoms() const;
    const Box* getBox() const;
    std::string getKey() const;
    std::vector<int> getAtomTypes()
    {
      return this->getPropertyValues<int>("type");
    }
    long int getAtomIdByIdx(const int vertexId) const override;
    long int getIdxByAtomId(const int atomId) const override;

    // computations
    double computeEndToEndDistance();
    double computeEndToEndDistanceWithDerivedImageFlags() const;
    double computeRadiusOfGyration();
    double computeRadiusOfGyrationWithDerivedImageFlags() const;
    double computeTotalMass();
    std::vector<double> computeBondLengths()
    {
      return AtomGraphParent::computeBondLengths(this->parent);
    };
    Eigen::Vector3d getOverallBondSum(const int crosslinkerType = 2) const;
    Eigen::Vector3d getOverallBondSumFromTo(
      size_t atomIdFrom,
      size_t atomIdTo,
      const int crosslinkerType = 2) const;
    size_t getNrOfBondsFromTo(size_t atomIdFrom,
                              size_t atomIdTo,
                              const int crosslinkerType = 2) const;

    // operators
    Atom operator[](size_t index) const
    {
      return this->getAtomByVertexIdx(index);
    }

  private:
    Box _boxNoUse;
    const Box* parent;
    MoleculeType typeOfThisMolecule;
    int size;
    std::string key;
    std::map<int, double> massPerType;
    std::unordered_map<long int, long int> atomIdToVertexIdx;

    void initializeFromGraph(const igraph_t* ingraph);
  };
} // namespace entities
} // namespace pylimer_tools

#endif
