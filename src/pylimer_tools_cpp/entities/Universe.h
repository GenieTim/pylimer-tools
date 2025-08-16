#ifndef UNIVERSE_H
#define UNIVERSE_H

extern "C"
{
#include <igraph/igraph.h>
}
#include "../utils/CerealUtils.h"
#include "Atom.h"
#include "AtomGraphParent.h"
#include "Molecule.h"
#include <Eigen/Dense>

#ifdef CEREALIZABLE
#include <cereal/types/base_class.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/unordered_map.hpp>
#endif

#include <map>
#include <unordered_map>
#include <vector>

namespace pylimer_tools::entities {

struct LoopIntersectionInfo
{
  std::vector<Atom> involvedAtoms;
  long int edge1;
  long int edge2;
  Eigen::Vector3d intersectionPoint;
  double direction;
};

struct StrandAffiliationInfo
{
  std::vector<long int> strandIdOfVertex;
  std::vector<long int> indexOfVertexInStrand;
};

class Universe : public AtomGraphParent
{
public:
  explicit Universe(double Lx = 1., double Ly = 1., double Lz = 1.);
  explicit Universe(const Box& box);

  // rule of three:
  // 1. destructor (to destroy the graph)
  ~Universe() override;
  // 2. copy constructor
  Universe(const Universe& src);
  // 3. copy assignment operator
  Universe& operator=(Universe src);

  // equality operator
  bool operator==(const Universe& other) const;
  bool operator!=(const Universe& other) const;

  // initilaization/setters (and removers)
  void setBoxLengths(double Lx,
                     double Ly,
                     double Lz,
                     bool rescaleAtomCoordinates = false);
  // atoms
  void addAtoms(const std::vector<long int>& ids,
                const std::vector<int>& types,
                const std::vector<double>& x,
                const std::vector<double>& y,
                const std::vector<double>& z,
                const std::vector<int>& nx,
                const std::vector<int>& ny,
                const std::vector<int>& nz);
  void addAtoms(
    const std::vector<long int>& ids,
    const std::vector<int>& types,
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& z,
    const std::vector<int>& nx,
    const std::vector<int>& ny,
    const std::vector<int>& nz,
    const std::unordered_map<std::string, std::vector<double>>& additionalData);
  // void addAtoms(const std::vector<Atom>& atoms);
  void removeAtoms(const std::vector<long int>& ids);
  void replaceAtom(long int id, const Atom& replacement);
  void replaceAtomType(long int id, int newType);
  // bonds
  void addBonds(const std::vector<long int>& from,
                const std::vector<long int>& to);
  void addBonds(const std::vector<long int>& from,
                const std::vector<long int>& to,
                const std::vector<int>& types);
  void addBonds(size_t NNewBonds,
                const std::vector<long int>& from,
                const std::vector<long int>& to);
  void addBonds(size_t NNewBonds,
                const std::vector<long int>& from,
                const std::vector<long int>& to,
                const std::vector<int>& bondTypes,
                bool ignoreNonExistentAtoms = false,
                bool simplify = false);
  void removeBonds(const std::vector<long int>& atomIdsFrom,
                   const std::vector<long int>& atomIdsTo);
  void removeBondsOfType(int bondType);
  // others
  void addAngles(const std::vector<long int>& from,
                 const std::vector<long int>& via,
                 const std::vector<long int>& to,
                 const std::vector<int>& types);
  void addDihedralAngles(const std::vector<long int>& from,
                         const std::vector<long int>& via1,
                         const std::vector<long int>& via2,
                         const std::vector<long int>& to,
                         const std::vector<int>& types);
  void setMasses(const std::map<int, double>& atomMassPerType);
  void setMassForType(int atomType, double mass);
  void setBox(const Box& box, bool rescaleAtomCoordinates = false);
  void setTimestep(const long int newTimestep)
  {
    this->timestep = newTimestep;
  };
  void initializeFromGraph(const igraph_t* ingraph);
  void removeAllAngles();
  void removeAllDihedralAngles();

  // adjustments
  void resampleVelocities(double mean,
                          double variance,
                          std::string seed = "",
                          bool is2d = false);
  void inferCoordinates(int crosslinkerType);
  void simplify();

  // getters
  [[nodiscard]] bool containsAtom(const Atom& atom) const;
  [[nodiscard]] bool containsAtomWithId(int atomId) const;
  [[nodiscard]] Atom getAtom(int atomId) const;
  [[nodiscard]] std::vector<Atom> getAtoms() const;
  // std::map<std::st¨ring, std::vector<long int>> getBonds() const;
  [[nodiscard]] std::map<std::string, std::vector<long int>> getAngles() const;
  [[nodiscard]] std::vector<double> computeAngles() const;
  [[nodiscard]] std::map<std::string, std::vector<long int>> getDihedralAngles()
    const;
  [[nodiscard]] std::vector<Universe> getClusters() const;
  [[nodiscard]] StrandAffiliationInfo getStrandAffiliation(
    int crosslinkerType) const;
  [[nodiscard]] std::vector<Molecule> getMolecules(
    int atomTypeToOmit = -1) const;
  [[nodiscard]] std::vector<Molecule> getChainsWithCrosslinker(
    int crossLinkerType) const;
  [[nodiscard]] std::vector<pylimer_tools::entities::MoleculeType>
  identifyObviouslyDanglingAtoms(bool distinguishFree = false) const;
  [[nodiscard]] Universe getNetworkOfCrosslinker(int crossLinkerType) const;
  [[nodiscard]] Universe contractVerticesAlongBondType(int bondType) const;
  // TODO: find & implement a better return type, e.g. std::vector<Molecule>
  [[nodiscard]] std::vector<std::vector<igraph_integer_t>> findLoops(
    int crossLinkerType,
    int maxLength = -1,
    bool skipSelfLoops = false,
    std::vector<std::vector<igraph_integer_t>>* edges = nullptr) const;
  [[nodiscard]] std::unordered_map<int, int> countLoopLengths(
    int maxLength = -1) const;
  [[nodiscard]] std::map<int, std::vector<std::vector<Atom>>> findLoopsOfAtoms(
    int crossLinkerType,
    int maxLength = -1,
    bool skipSelfLoops = false) const;
  [[nodiscard]] std::vector<Atom> findMinimalOrderLoopFrom(
    long int loopStart,
    long int loopStep1,
    int maxLength = -1,
    bool skipSelfLoops = false) const;
  [[nodiscard]] bool hasInfiniteStrand(int crossLinkerType,
                                       int maxLength = -1) const;
  [[nodiscard]] std::vector<int> getAtomTypes() const
  {
    return this->getPropertyValues<int>("type");
  }
  [[nodiscard]] std::map<int, int> countAtomTypes() const;
  [[nodiscard]] std::vector<size_t> countAtomsInSkinDistance(
    const std::vector<double>& distances,
    bool unwrapped = false) const;
  template<typename IN>
  [[nodiscard]] long int findVertexIdForProperty(const char* propertyName,
                                                 IN propertyValue) const;
  [[nodiscard]] Box getBox() const;
  [[nodiscard]] double getVolume() const;
  [[nodiscard]] size_t getNrOfAtoms() const;
  [[nodiscard]] size_t getNrOfBonds() const;
  [[nodiscard]] size_t getNrOfAngles() const;
  [[nodiscard]] size_t getNrOfDihedralAngles() const;
  [[nodiscard]] std::map<int, double> getMasses();
  [[nodiscard]] long int getTimestep() const { return this->timestep; };
  [[nodiscard]] long int getAtomIdByIdx(
    igraph_integer_t vertexId) const override;
  [[nodiscard]] igraph_integer_t getIdxByAtomId(long int atomId) const override;

  // operators
  Atom operator[](const size_t index) const { return this->getAtom(index); }

  // computations
  [[nodiscard]] std::map<std::string, std::vector<long int>> detectAngles()
    const;
  [[nodiscard]] std::map<std::string, std::vector<long int>>
  detectDihedralAngles() const;
  [[nodiscard]] std::map<int, int> determineFunctionalityPerType() const;
  [[nodiscard]] std::map<int, double> determineEffectiveFunctionalityPerType()
    const;
  [[nodiscard]] std::map<int, double> computeWeightFractions() const;
  [[nodiscard]] double computeWeightFractionOfClustersAssociatedWith(
    std::vector<long int> atomIds) const;
  [[nodiscard]] std::vector<std::pair<size_t, size_t>> interpolateEdges(
    int crossLinkerType,
    double interpolationFactor) const;
  [[nodiscard]] std::vector<double> computeDxs(
    const std::vector<long int>& bondFrom,
    const std::vector<long int>& bondTo) const;
  [[nodiscard]] std::vector<double> computeDys(
    const std::vector<long int>& bondFrom,
    const std::vector<long int>& bondTo) const;
  [[nodiscard]] std::vector<double> computeDzs(
    const std::vector<long int>& bondFrom,
    const std::vector<long int>& bondTo) const;
  [[nodiscard]] std::vector<double> computeBondLengths() const
  {
    return AtomGraphParent::computeBondLengths(this->box);
  };
  [[nodiscard]] std::vector<Eigen::Vector3d> computeBondVectors() const
  {
    return AtomGraphParent::computeBondVectors(this->box);
  };
  [[nodiscard]] double computeMeanSquaredBondLength() const
  {
    return AtomGraphParent::computeMeanSquaredBondLength(this->box);
  };
  [[nodiscard]] double computeTemperature(int dimensions = 3,
                                          double kb = 1.) const;
  std::vector<LoopIntersectionInfo> findLoopEntanglements(
    const std::vector<igraph_integer_t>& vertexIndicesLoop1,
    const std::vector<igraph_integer_t>& vertexIndicesLoop2,
    const std::vector<igraph_integer_t>& edgeIndicesLoop1,
    const std::vector<igraph_integer_t>& edgeIndicesLoop2) const;
  [[nodiscard]] double getMeanStrandLength(int crossLinkerType) const;
  [[nodiscard]] std::vector<double> computeEndToEndDistances(
    int crossLinkerType,
    bool implyImageFlags = false) const;
  [[nodiscard]] double computeMeanEndToEndDistance(
    int crossLinkerType,
    bool implyImageFlags = false) const;
  [[nodiscard]] double computeMeanSquareEndToEndDistance(
    int crossLinkerType,
    bool onlyThoseWithTwoCrosslinkers = false,
    bool implyImageFlags = false) const;
  [[nodiscard]] double computeMeanBondLength() const;
  [[nodiscard]] double computeTotalMass() const;
  [[nodiscard]] double computeTotalMassWithMasses(
    std::map<int, double> massPerTypeToUse) const;
  [[nodiscard]] double computeWeightAverageMolecularWeight(
    int crossLinkerType) const;
  [[nodiscard]] double computeNumberAverageMolecularWeight(
    int crossLinkerType) const;
  [[nodiscard]] double computePolydispersityIndex(int crossLinkerType) const;
  bool validate() const;

#ifdef CEREALIZABLE
  template<class Archive>
  void serialize(Archive& archive)
  {
    archive(cereal::virtual_base_class<AtomGraphParent>(this),
            // properties
            timestep,
            NAtoms,
            NBonds,
            box,
            // connectivity
            atomIdToVertexIdx,
            // angles etc.
            angleFrom,
            angleTo,
            angleVia,
            angleType,
            // dihedral angles etc.
            dihedralAngleFrom,
            dihedralAngleTo,
            dihedralAngleVia1,
            dihedralAngleVia2,
            dihedralAngleType,
            // type's properties
            massPerType);
  }
#endif

protected:
  // properties of the universe
  long int timestep = 0;
  size_t NAtoms = 0;
  size_t NBonds = 0;
  Box box;
  // connectivity
  // igraph_t graph;
  std::unordered_map<long int, igraph_integer_t> atomIdToVertexIdx;
  // extra info
  // TODO: might want to move the angle business to the parent?!?
  // angles (NOTE: only atom-ids, not vertex-idxs are used!)
  std::vector<long int> angleFrom = {};
  std::vector<long int> angleTo = {};
  std::vector<long int> angleVia = {};
  std::vector<int> angleType = {};
  // dihedral angles (NOTE: only atom-ids, not vertex-idxs are used!)
  std::vector<long int> dihedralAngleFrom = {};
  std::vector<long int> dihedralAngleVia1 = {};
  std::vector<long int> dihedralAngleVia2 = {};
  std::vector<long int> dihedralAngleTo = {};
  std::vector<int> dihedralAngleType = {};

  // type's properties
  std::map<int, double>
    massPerType = {}; // a dictionary with key: type, and value: weight per atom
  // of this atom type.

  // internal functions
  igraph_vs_t getVerticesOfType(int type) const;
  igraph_vs_t getVerticesByIndices(std::vector<igraph_integer_t> indices) const;
  std::vector<double> computeDs(const std::vector<long int>& bondFrom,
                                const std::vector<long int>& bondTo,
                                const std::string& direction,
                                double boxLimit) const;
  void resetAtomIdMapping();
};
}

#ifdef CEREALIZABLE
CEREAL_REGISTER_TYPE(pylimer_tools::entities::Universe);
CEREAL_REGISTER_POLYMORPHIC_RELATION(pylimer_tools::entities::AtomGraphParent,
                                     pylimer_tools::entities::Universe);
#endif

#endif
