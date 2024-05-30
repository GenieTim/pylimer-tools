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
#include <cereal/types/base_class.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/unordered_map.hpp>
#include <map>
#include <unordered_map>
#include <vector>

namespace pylimer_tools {
namespace entities {

  struct LoopIntersectionInfo
  {
    std::vector<Atom> involvedAtoms;
    long int edge1;
    long int edge2;
    Eigen::Vector3d intersectionPoint;
    double direction;
  };

  class Universe : public AtomGraphParent
  {
  public:
    Universe(const double Lx = 1., const double Ly = 1., const double Lz = 1.);
    Universe(const Box& box);

    // rule of three:
    // 1. destructor (to destroy the graph)
    ~Universe();
    // 2. copy constructor
    Universe(const Universe& src);
    // 3. copy assignment operator
    Universe& operator=(Universe src);

    // initilaization/setters (and removers)
    void setBoxLengths(const double Lx,
                       const double Ly,
                       const double Lz,
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
    void addAtoms(const std::vector<long int>& ids,
                  const std::vector<int>& types,
                  const std::vector<double>& x,
                  const std::vector<double>& y,
                  const std::vector<double>& z,
                  const std::vector<int>& nx,
                  const std::vector<int>& ny,
                  const std::vector<int>& nz,
                  const std::unordered_map<std::string, std::vector<double>>&
                    additionalData);
    void removeAtoms(const std::vector<long int>& ids);
    void replaceAtom(const long int id, const Atom& replacement);
    // bonds
    void addBonds(const std::vector<long int>& from,
                  const std::vector<long int>& to);
    void addBonds(const std::vector<long int>& from,
                  const std::vector<long int>& to,
                  const std::vector<int>& types);
    void addBonds(const size_t NNewBonds,
                  const std::vector<long int>& from,
                  const std::vector<long int>& to);
    void addBonds(const size_t NNewBonds,
                  const std::vector<long int>& from,
                  const std::vector<long int>& to,
                  const std::vector<int>& bondTypes,
                  const bool ignoreNonExistentAtoms = false,
                  const bool simplify = false);
    void removeBonds(const std::vector<long int>& atomIdsFrom,
                     const std::vector<long int>& atomIdsTo);
    void removeBondsOfType(const int bondType);
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
    void setMasses(const std::map<int, double>& massPerType);
    void setBox(const Box& box, bool rescaleAtomCoordinates = false);
    void setTimestep(long int timestep) { this->timestep = timestep; };
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
    bool containsAtomWithId(const int atomId) const;
    Atom getAtom(const int atomId) const;
    std::vector<Atom> getAtoms() const;
    // std::map<std::st¨ring, std::vector<long int>> getBonds() const;
    std::map<std::string, std::vector<long int>> getAngles() const;
    std::map<std::string, std::vector<long int>> getDihedralAngles() const;
    std::vector<Universe> getClusters() const;
    std::vector<Molecule> getMolecules(const int atomTypeToOmit = -1) const;
    std::vector<Molecule> getChainsWithCrosslinker(
      const int crossLinkerType) const;
    Universe getNetworkOfCrosslinker(const int crossLinkerType) const;
    // TODO: find & implement a better return type, e.g. std::vector<Molecule>
    std::vector<std::vector<long int>> findLoops(
      const int crossLinkerType,
      const int maxLength = -1,
      bool skipSelfLoops = false,
      std::vector<std::vector<long int>>* edges = nullptr) const;
    std::map<int, std::vector<std::vector<Atom>>> findLoopsOfAtoms(
      const int crossLinkerType,
      const int maxLength = -1,
      bool skipSelfLoops = false) const;
    std::vector<Atom> findMinimalOrderLoopFrom(
      const long int loopStart,
      const long int loopStep1,
      const int maxLength = -1,
      bool skipSelfLoops = false) const;
    bool hasInfiniteStrand(const int crossLinkerType,
                           const int maxLength = -1) const;
    std::vector<int> getAtomTypes() const
    {
      return this->getPropertyValues<int>("type");
    }
    std::map<int, int> countAtomTypes() const;
    std::vector<size_t> countAtomsInSkinDistance(std::vector<double> distances,
                                                 bool unwrapped = false) const;
    template<typename IN>
    long int findVertexIdForProperty(const char* propertyName,
                                     IN propertyValue) const;
    Box getBox() const;
    double getVolume() const;
    size_t getNrOfAtoms() const;
    size_t getNrOfBonds() const;
    size_t getNrOfAngles() const;
    size_t getNrOfDihedralAngles() const;
    std::map<int, double> getMasses();
    long int getTimestep() { return this->timestep; };
    long int getAtomIdByIdx(const int vertexId) const override;
    long int getIdxByAtomId(const int atomId) const override;

    // operators
    Atom operator[](size_t index) const { return this->getAtom(index); }

    // computations
    std::map<std::string, std::vector<long int>> detectAngles() const;
    std::map<std::string, std::vector<long int>> detectDihedralAngles() const;
    std::map<int, int> determineFunctionalityPerType() const;
    std::map<int, double> determineEffectiveFunctionalityPerType() const;
    std::map<int, double> computeWeightFractions() const;
    double computeWeightFractionOfClustersAssociatedWith(
      std::vector<long int> atomIds) const;
    std::vector<double> computeDxs(const std::vector<long int>& bondFrom,
                                   const std::vector<long int>& bondTo);
    std::vector<double> computeDys(const std::vector<long int>& bondFrom,
                                   const std::vector<long int>& bondTo);
    std::vector<double> computeDzs(const std::vector<long int>& bondFrom,
                                   const std::vector<long int>& bondTo);
    std::vector<double> computeBondLengths()
    {
      return AtomGraphParent::computeBondLengths(this->box);
    };
    double computeTemperature(const int dimensions = 3,
                              const double kb = 1.) const;
    Eigen::Vector3d getPositionVectorForVertex(const int vertexId) const;
    Eigen::Vector3d getUnwrappedPositionVectorForVertex(
      const int vertexId) const;
    std::vector<LoopIntersectionInfo> findLoopEntanglements(
      const std::vector<long int>& vertexIndicesLoop1,
      const std::vector<long int>& vertexIndicesLoop2,
      const std::vector<long int>& edgeIndicesLoop1,
      const std::vector<long int>& edgeIndicesLoop2) const;
    double getMeanStrandLength(int crossLinkerType);
    std::vector<double> computeEndToEndDistances(int crossLinkerType);
    double computeMeanEndToEndDistance(int crossLinkerType);
    double computeMeanSquareEndToEndDistance(
      int crossLinkerType,
      bool onlyThoseWithTwoCrosslinkers = false);
    double computeMeanBondLength();
    double computeTotalMass() const;
    double computeTotalMassWithMasses(
      std::map<int, double> massPerTypeToUse) const;
    double computeWeightAverageMolecularWeight(int crossLinkerType) const;
    double computeNumberAverageMolecularWeight(int crossLinkerType) const;
    double computePolydispersityIndex(int crossLinkerType) const;
    bool validate();

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

  protected:
    // properties of the universe
    long int timestep;
    size_t NAtoms = 0;
    size_t NBonds = 0;
    Box box;
    // connectivity
    // igraph_t graph;
    std::unordered_map<int, int> atomIdToVertexIdx;
    // extra info
    // TODO: might want to move the angle business to the parent?!?
    // angles (NOTE: only atom-ids, not vertex-idxs are used!)
    std::vector<long int> angleFrom;
    std::vector<long int> angleTo;
    std::vector<long int> angleVia;
    std::vector<int> angleType;
    // dihedral angles (NOTE: only atom-ids, not vertex-idxs are used!)
    std::vector<long int> dihedralAngleFrom;
    std::vector<long int> dihedralAngleVia1;
    std::vector<long int> dihedralAngleVia2;
    std::vector<long int> dihedralAngleTo;
    std::vector<int> dihedralAngleType;

    // type's properties
    std::map<int, double>
      massPerType; // a dictionary with key: type, and value: weight per atom
                   // of this atom type.

    // internal functions
    igraph_vs_t getVerticesOfType(const int type) const;
    std::vector<long int> getIndicesOfType(const int type) const;
    igraph_vs_t getVerticesByIndices(std::vector<long int> indices) const;
    std::vector<double> computeDs(const std::vector<long int>& bondFrom,
                                  const std::vector<long int>& bondTo,
                                  const std::string& direction,
                                  const double boxLimit) const;
  };
} // namespace entities
} // namespace pylimer_tools

CEREAL_REGISTER_TYPE(pylimer_tools::entities::Universe);
CEREAL_REGISTER_POLYMORPHIC_RELATION(pylimer_tools::entities::AtomGraphParent,
                                     pylimer_tools::entities::Universe);

#endif
