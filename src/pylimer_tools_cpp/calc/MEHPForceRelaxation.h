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
#define MAX 250000    /* maximum number of relaxation steps */
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

  double getAverageSpringLength() {
    return this->finalConfig.averageSpringLength;
  }

  int getFinalNrOfLoops() { return this->finalConfig.nrOfLoops; }

  double getInitialShearModulus() { return this->initialG; }

  double getFinalShearModulus() { return this->finalG; }

  double getInitialSquareDistance() { return this->initialS2; }

  double getFinalSquareDistance() { return this->finalS2; }

  double getInitialSquareRelativeDistance() { return this->initialS2len; }

  double getFinalSquareRelativeDistance() { return this->finalS2len; }

  void runForceRelaxation(int crosslinkerType) {

    long int i, step, Nact, Mact, a, b;
    Network net;
    if (!ConvertNetwork(&net, crosslinkerType)) {
      return;
    }
    double *u, *r, Fdef;
    double stress[3][3], r20, G0, r2, G, s20, s2, s20len, s2len;
    FILE *fp;
    const int M = this->universe.getMolecules(crosslinkerType).size();
    const int N = this->universe.getMeanStrandLength(crosslinkerType) + 1;
    const double bM = this->universe.computeMeanBondLength();
    const double Nb2 = N * bM * bM;
    const int f =
        this->universe.determineFunctionalityPerType()[crosslinkerType];

    for (int i = 0; i < net.nrOfSprings; i++) {

      if (std::fabs(net.springs[i].len) < 1.0e-10) {
        std::cout << "WARNING: " << i << std::endl;
      }
      // assert(net.springs[i].len != 0.0);
    }

    /* array allocation */

    u = (double *)calloc(3 * net.nrOfNodes, sizeof(double));
    r = (double *)calloc(3 * net.nrOfNodes, sizeof(double));

    /* initial */

    for (i = 0; i < 3 * net.nrOfNodes; i++) {
      u[i] = 0.0;
    }
    std::tie(s20, s20len) = computeStressAndSquareDistances(&net, u, stress);
    Fdef = Residual(&net, u, r);
    for (r20 = 0., i = 0; i < 3 * net.nrOfNodes; i++) {
      r20 += r[i] * r[i];
    }
    G0 = (stress[0][0] + stress[1][1] + stress[2][2]) / 3.;

    FILE *stepOutputFp;

    if (this->stepOutputFrequency > 0) {
      stepOutputFp = fopen(this->stepOutputFile.c_str(), "w");
      fprintf(stepOutputFp, "Step Fdef r2 G\n");
      fprintf(stepOutputFp, "%d %.10e %.16f %.16f\n", 0, Fdef, r20, G0);
    }

    /* force relaxation */
    double Gamma_eq = 0.0;
    double Rx2_sum = 0.0;
    double Ry2_sum = 0.0;
    double Rz2_sum = 0.0;
    for (step = 1; step <= MAX; step++) {
      Gamma_eq = 0.0;
      Rx2_sum = 0.0;
      Ry2_sum = 0.0;
      Rz2_sum = 0.0;
      for (i = 0; i < 3 * net.nrOfNodes; i++) {
        u[i] -= EPSILON * r[i];
      }
      Fdef = Residual(&net, u, r);
      for (r2 = 0., i = 0; i < 3 * net.nrOfNodes; i++) {
        r2 += r[i] * r[i];
        if (i % 3 == 0) {
          Rx2_sum += r[i] * r[i];
        } else if ((i + 1) % 3 == 0) {
          Ry2_sum += r[i] * r[i];
        } else {
          assert((i + 2) % 3 == 0);
          Rz2_sum += r[i] * r[i];
        }
        // Gamma_eq +=
      }
      Gamma_eq = r2 / ((double)net.nrOfSprings * Nb2);
      if (this->stepOutputFrequency > 0 &&
          (step <= 5 ||
           this->stepOutputFrequency * (step / this->stepOutputFrequency) ==
               step)) {
        computeStressAndSquareDistances(&net, u, stress);
        G = (stress[0][0] + stress[1][1] + stress[2][2]) / 3.;

        fprintf(stepOutputFp, "%ld %.10e %.16f %.16f\n", step, Fdef, r2,
                G);
      }

      if (r2 / r20 < RED2) {
        break;
      }
    }

    if (this->stepOutputFrequency > 0) {
      fclose(stepOutputFp);
    }

    // TODO: this is incorrect.
    double R2_mean = Gamma_eq / (double)net.nrOfSprings;
    double Gamma_x = 3 * Rx2_sum / ((double)net.nrOfSprings * Nb2);
    double Gamma_y = 3 * Ry2_sum / ((double)net.nrOfSprings * Nb2);
    double Gamma_z = 3 * Rz2_sum / ((double)net.nrOfSprings * Nb2);

    /* acquire equilibrium properties */
    std::tie(s2, s2len) = computeStressAndSquareDistances(&net, u, stress);
    G = (stress[0][0] + stress[1][1] + stress[2][2]) / 3.;

    /* output */
    if (this->outputEndNodes) {
      fp = fopen(this->endNodesFile.c_str(), "w");
      fprintf(fp, "%.10f %.10f %.10f\n", net.L[0], net.L[1], net.L[2]);
      fprintf(fp, "%ld #nodes\n", net.nrOfNodes);
      fprintf(fp, "%ld #springs\n", net.nrOfSprings);
      for (i = 0; i < net.nrOfNodes; i++) {
        fprintf(fp, "%.16f %.16f %.16f\n", net.nodes[i].x + u[3 * i],
                net.nodes[i].y + u[3 * i + 1], net.nodes[i].z + u[3 * i + 2]);
      }
      fclose(fp);
    }

    /* active springs */
    for (Nact = 0, i = 0; i < net.nrOfSprings; i++) {
      if (net.springs[i].isActive == 1) {
        ++Nact;
      }
    }

    /* active nodes */
    for (i = 0; i < net.nrOfNodes; i++)
      net.nodes[i].nrOfActiveConnected = 0; /* initial */
    for (i = 0; i < net.nrOfSprings; i++) {
      if (net.springs[i].isActive == 1) /* active spring */
      {
        a = net.springs[i].a;
        b = net.springs[i].b;
        ++(net.nodes[a].nrOfActiveConnected);
        ++(net.nodes[b].nrOfActiveConnected);
      }
    }
    for (Mact = 0, i = 0; i < net.nrOfNodes; i++) {
      if (net.nodes[i].nrOfActiveConnected >= 2) {
        ++Mact;
      }
    }

    /* save results */
    this->initialG = G0;
    this->initialS2 = s20;
    this->finalS2 = s2;
    this->finalS2len = s2len;
    this->initialS2len = s20len;
    this->finalG = G;
    this->finalConfig = net;
    this->nrOfActiveNodes = Mact;
    this->nrOfActiveSprings = Nact;

    /** array deallocation */
    free(u);
    free(r);
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
   * @brief the residual is minus the assembled force vector
   *
   * @param net the network to compute the residual for
   * @param u the displacement vector
   * @param r the residual
   * @return the definitive force
   */
  double Residual(Network *net, double *u, double *r) {
    long int i, a, b, len, j, glob[6];
    double kappa, s[3], ua[3], ub[3], s2, fele[6];

    /* initial */
    double Fdef = 0.;
    for (i = 0; i < 3 * net->nrOfNodes; i++) {
      r[i] = 0.;
    }

    /* assembly */
    for (i = 0; i < net->nrOfSprings; i++) {
      a = net->springs[i].a;
      b = net->springs[i].b;
      len = net->springs[i].len;
      kappa = 3. / len;

      actualSpringDistance(net->nodes[b], b, net->nodes[a], a, u, s, net->L);

      /* energy update */
      for (s2 = 0., j = 0; j < 3; j++) {
        s2 += s[j] * s[j];
      }
      Fdef += 0.5 * kappa * s2;

      /* element force vector */
      for (j = 0; j < 3; j++) {
        fele[j] = kappa * s[j];
        fele[j + 3] = -kappa * s[j];
      }

      /* global numbering scheme */
      for (j = 0; j < 3; j++) {
        glob[j] = 3 * a + j;
        glob[j + 3] = 3 * b + j;
      }

      /* residual update */
      for (j = 0; j < 6; j++) {
        r[glob[j]] -= fele[j];
      }
    };

    return Fdef;
  }

  std::pair<double, double>
  computeStressAndSquareDistances(Network *net, double *u,
                                  double stress[3][3]) {
    long int i, a, b, len, j, k;
    double kappa, s[3], ua[3], ub[3], s2 = 0, s2len = 0;

    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        stress[j][k] = 0.;
      }
    }

    /* assembly */
    for (i = 0; i < net->nrOfSprings; i++) {
      a = net->springs[i].a;
      b = net->springs[i].b;
      len = net->springs[i].len;
      if (isnan(a) || isnan(b) || isnan(len)) {
        std::cout << "ERROR: a or b or len is nan" << a << " " << b << " "
                  << len << std::endl;
      }
      // if (len == 0) {
      //   std::cout << "ERROR: len is 0 " << len << " between a or b " << a <<
      //   " "
      //             << b << " " << std::endl;
      //   kappa = 100;
      // }
      if (std::fabs(len) > 1e-10) {
        kappa = 3. / len;
      } else {
        kappa = 3e5;
      }

      actualSpringDistance(net->nodes[b], b, net->nodes[a], a, u, s, net->L);

      /* spring contribution to the overall stress tensor */
      for (j = 0; j < 3; j++) {
        for (k = 0; k < 3; k++) {
          stress[j][k] += kappa * s[j] * s[k];
        }
      }

      double s2local = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);

      /* update */
      s2 += s2local;
      s2len += s2local / len;

      /* loop count */
      net->springs[i].isActive = 1;
      if (s2 < EPS2) {
        net->nrOfLoops++;
        net->springs[i].isActive = false;
      }
    };

    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        stress[j][k] /= net->vol;
      }
    }

    return std::make_pair(s2 / (double)net->nrOfSprings,
                          s2len / (double)net->nrOfSprings);
  }

  /**
   * @brief Get the actual distance between two nodes
   *
   * @param a node a
   * @param indexA the index of the first node (needed for indexing into u)
   * @param b
   * @param indexB the index of the second node (needed for indexing into u)
   * @param u the displacement vector
   * @param coords the vector to write the coordinates into
   * @param boxL the box side lengths
   */
  void actualSpringDistance(_Node a, int indexA, _Node b, int indexB, double *u,
                            double *coords, double *boxL) {
    double s[3], ua[3], ub[3];

    /* initial spring vector */
    coords[0] = a.x - b.x;
    coords[1] = a.y - b.y;
    coords[2] = a.z - b.z;

    /* spring displacement vectors */
    for (size_t j = 0; j < 3; j++) {
      ua[j] = u[3 * indexA + j];
      ub[j] = u[3 * indexB + j];
    }

    /* current spring vector */
    for (size_t j = 0; j < 3; j++) {
      coords[j] += ub[j] - ua[j];
    }

    /* periodic boundary conditions */
    ImposePBC(coords, boxL);
  }

private:
  pylimer_tools::entities::Universe universe;
  int stepOutputFrequency = 0;
  std::string stepOutputFile;
  bool outputEndNodes = false;
  std::string endNodesFile;
  _Network finalConfig;
  double finalG = 0.0;
  double initialG = 0.0;
  int nrOfActiveSprings = 0;
  int nrOfActiveNodes = 0;
  double initialS2 = 0.0;
  double finalS2 = 0.0;
  double finalS2len = 0.0;
  double initialS2len = 0.0;
};
} // namespace mehp
} // namespace calc
} // namespace pylimer_tools
#endif
