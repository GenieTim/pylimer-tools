#ifndef MEHP_FORCE_RELAX_H
#define MEHP_FORCE_RELAX_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <map>
#include <string>
#include <tuple>

namespace pylimer_tools {
namespace calc {
namespace mehp {
#define EPSILON 0.077 /* relaxation parameter, adjusted by trial and errors */
#define RED2 1.e-8    /* second residual norm reduction */
#define EPS2 1.e-16   /* distance tolerance squared for the loop count */

typedef struct _Spring {
  long int a;    /* first node */
  long int b;    /* second node */
  long int len;  /* length */
  bool isActive; /* 1/0; active or not */
} Spring;

typedef struct _Node {
  double x; /* x-coordinate */
  double y; /* y-coordinate */
  double z; /* z-coordinate */
  long int
      nrOfActiveConnected; /* number of active springs connected to the node */
} Node;

typedef struct _Network {
  double L[3];                /* box sizes */
  double vol;                 /* box volume */
  long int nrOfNodes;         /* number of nodes */
  long int nrOfSprings;       /* number of springs */
  Node *nodes;                /* nodes */
  Spring *springs;            /* spings */
  double *nodalDisplacements; /* nodal displacements */
  double *nodalForces;        /* nodal forces */
  double averageSpringLength; /* average spring length */
  long int nrOfLoops;         /* loops */
} Network;

// heavily inspired by Prof. Dr. Andrei Gusev's Code
class MEHPForceRelaxation {
public:
  MEHPForceRelaxation(const pylimer_tools::entities::Universe u,
                      int crosslinkerType = 2)
      : universe(u) {
    // interpret network already to be able to give early results
    Network net;
    ConvertNetwork(&net, crosslinkerType);
    free(net.nodes);
    free(net.springs);
    this->finalConfig = net;
  };

  void configDoOutputSteps(const std::string outputFile, const int outputFreq) {
    this->stepOutputFile = outputFile;
    this->stepOutputFrequency = outputFreq;
  }

  void configDoOutputFinalCoordinates(const std::string outputFile) {
    this->outputEndNodes = true;
    this->endNodesFile = outputFile;
  }

  double getVolume() { return this->finalConfig.vol; }

  int getNrOfNodes() { return this->finalConfig.nrOfNodes; }

  int getNrOfSprings() { return this->finalConfig.nrOfSprings; }

  int getNrOfActiveNodes() { return this->nrOfActiveNodes; }

  int getNrOfActiveSprings() { return this->nrOfActiveSprings; }

  double getAverageSpringLength() { return sqrt(this->R2Mean); }

  double getGammaEq() { return this->gammaEq; }

  double getGammaX() { return this->gammaX; }

  double getGammaY() { return this->gammaY; }

  double getGammaZ() { return this->gammaZ; }

  double getNb2() { return this->Nb2; }

  void runForceRelaxation(int crosslinkerType, long int maxNrOfSteps = 250000,
                          double tol = 1e-8, double Nb2 = -1.0,
                          bool is2d = false, double eps = 0.025,
                          double kappa = 1.0) {
    this->is2d = is2d;
    int dimensions = 3;
    if (is2d) {
      dimensions = 2;
    }
    Network net;
    if (!ConvertNetwork(&net, crosslinkerType)) {
      return;
    }
    const int M = this->universe.getMolecules(crosslinkerType).size();
    const int N = this->universe.getMeanStrandLength(crosslinkerType) + 1;
    const double bM = this->universe.computeMeanBondLength();
    if (Nb2 == -1.0) {
      this->Nb2 = N * bM * bM;
    } else {
      this->Nb2 = Nb2;
    }
    const int f =
        this->universe.determineFunctionalityPerType()[crosslinkerType];

    for (size_t i = 0; i < net.nrOfSprings; i++) {

      if (std::fabs(net.springs[i].len) < 1.0e-10) {
        std::cout << "WARNING: Spring " << i << " has negligible length."
                  << std::endl;
      }
      // assert(net.springs[i].len != 0.0);
    }

    /* array allocation */
    double *force;
    force = (double *)calloc(3 * net.nrOfNodes, sizeof(double));
    for (size_t i = 0; i < 3 * net.nrOfNodes; ++i) {
      force[i] = 0.0;
    }

    // calculate initial absolute force
    for (size_t i = 0; i < net.nrOfSprings; ++i) {
      _Spring spring = net.springs[i];
      double distances[3];
      actualSpringDistance(net.nodes[spring.a], net.nodes[spring.b], distances,
                           net.L);
      for (int j = 0; j < dimensions; ++j) {
        force[spring.a * 3 + j] -= kappa * distances[j];
        force[spring.b * 3 + j] += kappa * distances[j];
      }
    }

    double absoluteForce = 0.0;
    for (size_t i = 0; i < 3 * net.nrOfNodes; ++i) {
      absoluteForce += force[i] * force[i];
    }
    double initialForce = absoluteForce;

    FILE *stepOutputFp;
    if (this->stepOutputFrequency > 0) {
      stepOutputFp = fopen(this->stepOutputFile.c_str(), "w");
      fprintf(stepOutputFp, "Step Fabs GammaEq\n");
      fprintf(stepOutputFp, "%d %.16f %.16f\n", 0, initialForce, 0.0);
    }

    // start of force relaxation
    long int stepsDone = 0;
    double Gamma_eq = 0.0;
    while (absoluteForce / initialForce > tol && stepsDone < maxNrOfSteps) {
      stepsDone += 1;
      // update coordinates
      for (size_t i = 0; i < net.nrOfNodes; ++i) {
        double coordinates[3];
        coordinates[0] = net.nodes[i].x + eps * force[i * 3 + 0];
        coordinates[1] = net.nodes[i].y + eps * force[i * 3 + 1];
        coordinates[2] = net.nodes[i].z + eps * force[i * 3 + 2];
        ImposePBC(coordinates, net.L);
        net.nodes[i].x = coordinates[0];
        net.nodes[i].y = coordinates[1];
        net.nodes[i].z = coordinates[2];
      }

      // calculate new force
      Gamma_eq = 0.0;
      for (size_t i = 0; i < net.nrOfNodes; ++i) {
        force[i] = 0.0;
      }
      for (size_t i = 0; i < net.nrOfSprings; ++i) {
        _Spring spring = net.springs[i];
        double distances[3];
        actualSpringDistance(net.nodes[spring.a], net.nodes[spring.b],
                             distances, net.L);
        for (int j = 0; j < dimensions; ++j) {
          force[spring.a * 3 + j] -= kappa * distances[j];
          force[spring.b * 3 + j] += kappa * distances[j];
          Gamma_eq += distances[j] * distances[j];
        }
      }
      absoluteForce = 0.0;
      for (int i = 0; i < 3 * net.nrOfNodes; ++i) {
        absoluteForce += force[i] * force[i];
      }

      // potentially output
      if (this->stepOutputFrequency > 0 &&
          (stepsDone <= 5 || this->stepOutputFrequency *
                                     (stepsDone / this->stepOutputFrequency) ==
                                 stepsDone)) {

        fprintf(stepOutputFp, "%ld %.16f %.16f\n", stepsDone, absoluteForce,
                Gamma_eq);
      }
    }

    // calculate equilibrium values
    Gamma_eq = 0.0;
    double Rx2_sum = 0.0;
    double Ry2_sum = 0.0;
    double Rz2_sum = 0.0;
    int Nact = 0;

    for (size_t i = 0; i < net.nrOfSprings; ++i) {
      _Spring spring = net.springs[i];
      double distances[3];
      actualSpringDistance(net.nodes[spring.a], net.nodes[spring.b], distances,
                           net.L);
      double springLen = 0.0;
      for (int j = 0; j < dimensions; ++j) {
        springLen += distances[j] * distances[j];
      }
      Gamma_eq += springLen;
      Rx2_sum += distances[0] * distances[0];
      Ry2_sum += distances[1] * distances[1];
      Rz2_sum += distances[2] * distances[2];
      if (springLen / this->Nb2 > 0.001) {
        Nact += 1;
        net.springs[i].isActive = true;
      }
    }

    this->gammaEq = Gamma_eq / ((double)net.nrOfSprings * this->Nb2);
    this->R2Mean = Gamma_eq / (double)net.nrOfSprings;
    this->gammaX = 3 * Rx2_sum / ((double)net.nrOfSprings * this->Nb2);
    this->gammaY = 3 * Ry2_sum / ((double)net.nrOfSprings * this->Nb2);
    this->gammaZ = 3 * Rz2_sum / ((double)net.nrOfSprings * this->Nb2);

    /* output */
    if (this->outputEndNodes) {
      FILE *fp = fopen(this->endNodesFile.c_str(), "w");
      fprintf(fp, "%.10f %.10f %.10f\n", net.L[0], net.L[1], net.L[2]);
      fprintf(fp, "%ld #nodes\n", net.nrOfNodes);
      fprintf(fp, "%ld #springs\n", net.nrOfSprings);
      for (size_t i = 0; i < net.nrOfNodes; i++) {
        fprintf(fp, "%.16f %.16f %.16f\n", net.nodes[i].x, net.nodes[i].y,
                net.nodes[i].z);
      }
      fclose(fp);
    }

    /* count active nodes */
    int Mact = 0;
    for (size_t i = 0; i < net.nrOfNodes; i++) {
      net.nodes[i].nrOfActiveConnected = 0; /* initial */
    }
    for (size_t i = 0; i < net.nrOfSprings; i++) {
      if (net.springs[i].isActive == true) /* active spring */
      {
        size_t a = net.springs[i].a;
        size_t b = net.springs[i].b;
        ++(net.nodes[a].nrOfActiveConnected);
        ++(net.nodes[b].nrOfActiveConnected);
      }
    }
    for (size_t i = 0; i < net.nrOfNodes; i++) {
      if (net.nodes[i].nrOfActiveConnected >= 2) {
        ++Mact;
      }
    }

    /* save results */
    this->finalConfig = net;
    this->nrOfActiveNodes = Mact;
    this->nrOfActiveSprings = Nact;

    /** array deallocation */
    free(force);
    free(net.nodes);
    free(net.springs);
  };

protected:
  /**
   * @brief Convert the universe to a network
   *
   * @param net the target network
   * @param crosslinkerType the atom type of the crosslinker
   * @return true
   * @return false
   */
  bool ConvertNetwork(Network *net, int crosslinkerType = 2) {
    pylimer_tools::entities::Universe crosslinkerUniverse =
        this->universe.getNetworkOfCrosslinker(crosslinkerType);
    // crosslinkerUniverse.simplify();
    pylimer_tools::entities::Box box = crosslinkerUniverse.getBox();
    net->L[0] = box.getLx();
    net->L[1] = box.getLy();
    net->L[2] = box.getLz();
    net->nrOfNodes = crosslinkerUniverse.getNrOfAtoms();
    net->nrOfSprings = crosslinkerUniverse.getNrOfBonds();

    int usualChainLen =
        this->universe.getMolecules(crosslinkerType)[0].getNrOfAtoms();

    net->nodes = (Node *)calloc(net->nrOfNodes, sizeof(Node));
    net->springs = (Spring *)calloc(net->nrOfSprings, sizeof(Spring));

    // convert beads
    std::vector<pylimer_tools::entities::Atom> allAtoms =
        crosslinkerUniverse.getAtoms();
    std::map<int, int> atomIdToNode;
    for (int i = 0; i < allAtoms.size(); ++i) {
      pylimer_tools::entities::Atom atom = allAtoms[i];
      net->nodes[i].x = atom.getX();
      net->nodes[i].y = atom.getY();
      net->nodes[i].z = atom.getZ();
      atomIdToNode[atom.getId()] = i;
    }

    // convert springs
    std::map<std::string, std::vector<long int>> allBonds =
        crosslinkerUniverse.getBonds();
    net->averageSpringLength = 0;
    for (int i = 0; i < net->nrOfSprings; ++i) {
      int atomIdFrom = allBonds["bond_from"][i];
      int atomIdTo = allBonds["bond_to"][i];
      net->springs[i].a = atomIdToNode.at(atomIdFrom);
      net->springs[i].b = atomIdToNode.at(atomIdTo);
      net->springs[i].len = usualChainLen;
      net->averageSpringLength += usualChainLen;
    }

    if (crosslinkerUniverse.getNrOfBonds() != net->nrOfSprings) {
      return false;
    };

    net->averageSpringLength /= net->nrOfSprings;

    /* box volume */
    net->vol = net->L[0] * net->L[1] * net->L[2];
    return true;
  };

  /**
   * @brief Adjust a vector of distances to lie within half the box
   *
   * @param s the distances
   * @param box the box lengths
   */
  void ImposePBC(double s[3], double box[3]) {
    long int i;
    double half;

    for (i = 0; i < 3; i++) {
      half = 0.5 * box[i];
      while (s[i] > half) {
        s[i] -= box[i];
      }
      while (s[i] < -half) {
        s[i] += box[i];
      }
    }

    return;
  }

  /**
   * @brief Compute the distance between two nodes
   *
   * @param a the first node
   * @param b the second node
   * @param coords the array to write the distances in
   * @param boxL the box sizes
   */
  void actualSpringDistance(_Node a, _Node b, double (&coords)[3],
                            double *boxL) {
    /* initial spring vector */
    coords[0] = a.x - b.x;
    coords[1] = a.y - b.y;
    if (this->is2d) {
      coords[2] = 0.0;
    } else {
      coords[2] = a.z - b.z;
    }

    /* periodic boundary conditions */
    ImposePBC(coords, boxL);
  }

private:
  pylimer_tools::entities::Universe universe;
  bool is2d = false;
  int stepOutputFrequency = 0;
  std::string stepOutputFile;
  bool outputEndNodes = false;
  std::string endNodesFile;
  _Network finalConfig;
  int nrOfActiveSprings = 0;
  int nrOfActiveNodes = 0;
  double Nb2 = 0.0;
  double gammaEq = 0.0;
  double gammaX = 0.0;
  double gammaY = 0.0;
  double gammaZ = 0.0;
  double R2Mean = 0.0;
};
} // namespace mehp
} // namespace calc
} // namespace pylimer_tools
#endif
