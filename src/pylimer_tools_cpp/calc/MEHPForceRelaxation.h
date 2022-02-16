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

namespace pylimer_tools {
namespace calc {
namespace mehp {
#define OUTFREQ 250   /* output frequency */
#define MAX 250000    /* maximum number of relaxation steps */
#define EPSILON 0.077 /* relaxation parameter, adjusted by trial and errors */
#define RED2 1.e-8    /* second residual norm reduction */
#define OUTNODES 1    /* write output nodes_eq.dat file */
#define EPS2 1.e-16   /* distance tolerance squared for the loop count */

typedef int Integer; /* large enough integer */

typedef struct _Spring {
  Integer a;   /* first node */
  Integer b;   /* second node */
  Integer len; /* length */
  Integer isActive; /* 1/0; active or not */
} Spring;

typedef struct _Node {
  double x;    /* x-coordinate */
  double y;    /* y-coordinate */
  double z;    /* z-coordinate */
  Integer nrOfActiveConnected; /* number of active springs connected to the node */
} Node;

typedef struct _Network {
  double L[3];     /* box sizes */
  double vol;      /* box volume */
  Integer nrOfNodes;   /* number of nodes */
  Integer nrOfSprings; /* number of springs */
  Node *nodes;       /* nodes */
  Spring *springs;     /* spings */
  double *nodalDisplacements;       /* nodal displacements */
  double *nodalForces;   /* nodal forces */
  double averageSpringLength;    /* average spring length */
  Integer nrOfLoops;   /* loops */
} Network;

// heavily inspired by Prof. Dr. Andrei Gusev"s Code
class MEHPForceRelaxation {
public:
  MEHPForceRelaxation(const pylimer_tools::entities::Universe u)
      : universe(u){};
  void runForceRelaxation(int crosslinkerType = 2) {

    Integer i, step, Nact, Mact, a, b;
    Network net;
    double *u, *r, Fdef;
    double stress[3][3], r20, G0, r2, G, s20, s2, s20len, s2len;
    FILE *fp;
    const int M = this->universe.getMolecules(crosslinkerType).size();
    const int N = this->universe.getMeanStrandLength(crosslinkerType) + 1;
    const double bM = this->universe.computeMeanBondLength();
    const double Nb2 = N * bM * bM;
    const int f =
        this->universe.determineFunctionalityPerType()[crosslinkerType];

    /* read network */
    if (!ConvertNetwork(&net, crosslinkerType)) {
      return;
    };

    for (int i = 0; i < net.nrOfSprings; i++) {

      if (std::fabs(net.springs[i].len) < 1.0e-10) {
        std::cout << "WARNING: " << i << std::endl;
      }
      assert(net.springs[i].len != 0.0);
    }

    /* array allocation */

    u = (double *)calloc(3 * net.nrOfNodes, sizeof(double));
    r = (double *)calloc(3 * net.nrOfNodes, sizeof(double));

    /* initial */

    for (i = 0; i < 3 * net.nrOfNodes; i++) {
      u[i] = 0.0;
    }
    Stress(&net, u, stress);
    printf("\nInitial stress tensor:\n");
    for (i = 0; i < 3; i++) {
      printf("%.10f %.10f %.10f\n", stress[i][0], stress[i][1], stress[i][2]);
    }
    Residual(&net, u, r, &Fdef);
    for (r20 = 0., i = 0; i < 3 * net.nrOfNodes; i++) {
      r20 += r[i] * r[i];
    }
    G0 = (stress[0][0] + stress[1][1] + stress[2][2]) / 3.;
    printf("\nForce relaxation:\n");
    printf("%d %.10e %.16f %.16f\n", 0, Fdef, r20, G0);

    fp = fopen("control.dat", "w");
    fprintf(fp, "%d %.10e %.16f %.16f\n", 0, Fdef, r20, G0);
    fclose(fp);

    s20 = S2(&net, u);
    s20len = S2len(&net, u);

    printf("%f\n", Fdef);

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
      Residual(&net, u, r, &Fdef);
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

      if (r2 / r20 < RED2) {
        break;
      }
      Gamma_eq = r2 / ((double)net.nrOfSprings * Nb2);
      if (step <= 5 || OUTFREQ * (step / OUTFREQ) == step) {
        Stress(&net, u, stress);
        G = (stress[0][0] + stress[1][1] + stress[2][2]) / 3.;

        printf("%d %.10e %.16f %.16f\n", step, Fdef, r2 / r20, G);
        fp = fopen("control.dat", "a");
        fprintf(fp, "%d %.10e %.16f %.16f\n", step, Fdef, r2 / r20, G);
        fclose(fp);
      }
    }
    // TODO: this is incorrect.
    double R2_mean = Gamma_eq / (double)net.nrOfSprings;
    double Gamma_x = 3 * Rx2_sum / ((double)net.nrOfSprings * Nb2);
    double Gamma_y = 3 * Ry2_sum / ((double)net.nrOfSprings * Nb2);
    double Gamma_z = 3 * Rz2_sum / ((double)net.nrOfSprings * Nb2);

    /* acquire equilibrium properties */
    Stress(&net, u, stress);
    G = (stress[0][0] + stress[1][1] + stress[2][2]) / 3.;
    printf("\n%d %.10e %.16f %.16f\n", --step, Fdef, r2, G);
    fp = fopen("control.dat", "a");
    fprintf(fp, "\n%d %.10e %.16f %.16f\n", step, Fdef, r2, G);
    fclose(fp);
    printf("\nStress tensor:\n");
    for (i = 0; i < 3; i++) {
      printf("%.10f %.10f %.10f\n", stress[i][0], stress[i][1], stress[i][2]);
    }
    /* output */
    if (OUTNODES) {
      fp = fopen("nodes_eq.dat", "w");
      fprintf(fp, "%.10f %.10f %.10f\n", net.L[0], net.L[1], net.L[2]);
      fprintf(fp, "%d #nodes\n", net.nrOfNodes);
      fprintf(fp, "%d #springs\n", net.nrOfSprings);
      for (i = 0; i < net.nrOfNodes; i++) {
        fprintf(fp, "%.16f %.16f %.16f\n", net.nodes[i].x + u[3 * i],
                net.nodes[i].y + u[3 * i + 1], net.nodes[i].z + u[3 * i + 2]);
      }
      fclose(fp);
    }

    s2 = S2(&net, u);
    s2len = S2len(&net, u);

    fp = fopen("out.dat", "w");
    fprintf(fp, "%.10f %.10f %.10f %.10f   %.10f %.10f\n",
            net.nrOfSprings / net.vol, net.averageSpringLength, G0, G, s20, s2);
    fclose(fp);

    fp = fopen("s2len.dat", "w");
    fprintf(fp, "%.10f  %.10f %.10f  %d\n", net.nrOfSprings / net.vol, s20len,
            s2len, net.nrOfLoops);
    fclose(fp);

    /* active springs */
    for (Nact = 0, i = 0; i < net.nrOfSprings; i++) {
      if (net.springs[i].isActive == 1) {
        ++Nact;
      }
    }
    printf("\nactive springs: %d\n", Nact);

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
    printf("active nodes: %d\n", Mact);

    /* OUTPUT FOR THE MAIN FIGURES */
    /* data as in Fabian"s MATLAB, approx */
    std::string filename = "force_relax_info_M" + std::to_string(M) + "_N" +
                           std::to_string(N) + "_f" + std::to_string(f) +
                           ".dat";
    FILE *fileID = fopen(filename.c_str(), "w");
    fprintf(fileID, "# Gamma x | Gamma y | Gamma z | Gamma_eq | R2_mean \n");
    fprintf(fileID, "%10.8f %10.8f %10.8f %10.8f %10.8f\n", Gamma_x, Gamma_y,
            Gamma_z, Gamma_eq, R2_mean);
    fclose(fileID);

    /* exact and ANT shear moduli [MPa] as functions of the total strand number
     * density [nm^-3] */
    fp = fopen("ANT.dat", "w");
    fprintf(fp, "%.10f   %.10f %.10f\n", net.nrOfSprings / net.vol, 4.14195 * G,
            4.14195 * s2len * net.nrOfSprings / net.vol);
    fclose(fp);

    /* exact and ANM shear moduli [MPa] using total and active strand densities
     */
    fp = fopen("ANM.dat", "w");
    fprintf(fp, "%.10f   %.10f %.10f %.10f\n", net.nrOfSprings / net.vol,
            4.14195 * G, 4.14195 * net.nrOfSprings / net.vol,
            4.14195 * Nact / net.vol);
    fclose(fp);

    /* exact and PNM shear moduli [MPa] using total and active strand and node
     * densities */
    fp = fopen("PNM.dat", "w");
    fprintf(fp, "%.10f   %.10f %.10f %.10f\n", net.nrOfSprings / net.vol,
            4.14195 * G, 4.14195 * (net.nrOfSprings - net.nrOfNodes) / net.vol,
            4.14195 * (Nact - Mact) / net.vol);
    fclose(fp);

    /* total and active node and strand number densities */
    fp = fopen("Densities.dat", "w");
    fprintf(fp, "%.10f %.10f     %.10f %.10f\n", net.nrOfSprings / net.vol,
            Nact / net.vol, net.nrOfNodes / net.vol, Mact / net.vol);
    fclose(fp);
  };

protected:
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

    assert(crosslinkerUniverse.getNrOfBonds() == net->nrOfSprings);

    net->averageSpringLength /= net->nrOfSprings;

    /* box volume */
    net->vol = net->L[0] * net->L[1] * net->L[2];
    return true;
  };

  void ImposePBC(double s[3], double box[3]) {
    Integer i;
    double half;

    for (i = 0; i < 3; i++) {
      half = 0.5 * box[i];
      if (s[i] > half)
        s[i] -= box[i];
      if (s[i] < -half)
        s[i] += box[i];
    }

    return;
  }

  /* the residual is minus the assembled force vector */
  void Residual(Network *net, double *u, double *r, double *Fdef) {
    Integer i, a, b, len, j, glob[6];
    double kappa, s[3], ua[3], ub[3], s2, fele[6];

    /* initial */
    *Fdef = 0.;
    for (i = 0; i < 3 * net->nrOfNodes; i++)
      r[i] = 0.;

    /* assembly */
    for (i = 0; i < net->nrOfSprings; i++) {
      a = net->springs[i].a;
      b = net->springs[i].b;
      len = net->springs[i].len;
      kappa = 3. / len;

      /* initial spring vector */
      s[0] = net->nodes[b].x - net->nodes[a].x;
      s[1] = net->nodes[b].y - net->nodes[a].y;
      s[2] = net->nodes[b].z - net->nodes[a].z;

      /* spring displacement vectors */
      for (j = 0; j < 3; j++) {
        ua[j] = u[3 * a + j];
        ub[j] = u[3 * b + j];
      }

      /* current spring vector */
      for (j = 0; j < 3; j++)
        s[j] += ub[j] - ua[j];

      /* periodic boundary conditions */
      ImposePBC(s, net->L);

      /* energy update */
      for (s2 = 0., j = 0; j < 3; j++)
        s2 += s[j] * s[j];
      *Fdef += 0.5 * kappa * s2;

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
      for (j = 0; j < 6; j++)
        r[glob[j]] -= fele[j];
    };

    return;
  }

  void Stress(Network *net, double *u, double stress[3][3]) {
    Integer i, a, b, len, j, k;
    double kappa, s[3], ua[3], ub[3];

    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
        stress[j][k] = 0.;

    /* assembly */
    for (i = 0; i < net->nrOfSprings; i++) {
      a = net->springs[i].a;
      b = net->springs[i].b;
      len = net->springs[i].len;
      if (isnan(a) || isnan(b) || isnan(len)) {
        std::cout << "ERROR: a or b or len is nan" << a << " " << b << " "
                  << len << std::endl;
      }
      kappa = 3. / len;
      if (len == 0) {
        std::cout << "ERROR: len is 0 " << len << " between a or b " << a << " "
                  << b << " " << std::endl;
        kappa = 100;
      }

      /* initial spring vector */
      s[0] = net->nodes[b].x - net->nodes[a].x;
      s[1] = net->nodes[b].y - net->nodes[a].y;
      s[2] = net->nodes[b].z - net->nodes[a].z;

      /* spring displacement vectors */
      for (j = 0; j < 3; j++) {
        ua[j] = u[3 * a + j];
        ub[j] = u[3 * b + j];
      }

      /* current spring vector */
      for (j = 0; j < 3; j++)
        s[j] += ub[j] - ua[j];

      /* periodic boundary conditions */
      ImposePBC(s, net->L);

      /* spring contribution to the overall stress tensor */
      for (j = 0; j < 3; j++)
        for (k = 0; k < 3; k++)
          stress[j][k] += kappa * s[j] * s[k];
    };

    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
        stress[j][k] /= net->vol;

    return;
  }

  /* average ratio <s2> */
  double S2(Network *net, double *u) {
    Integer i, a, b, len, j;
    double kappa, s[3], ua[3], ub[3], s2;

    /* assembly */
    for (s2 = 0., i = 0; i < net->nrOfSprings; i++) {
      a = net->springs[i].a;
      b = net->springs[i].b;
      len = net->springs[i].len;
      kappa = 3. / len;

      /* initial spring vector */
      s[0] = net->nodes[b].x - net->nodes[a].x;
      s[1] = net->nodes[b].y - net->nodes[a].y;
      s[2] = net->nodes[b].z - net->nodes[a].z;

      /* spring displacement vectors */
      for (j = 0; j < 3; j++) {
        ua[j] = u[3 * a + j];
        ub[j] = u[3 * b + j];
      }

      /* current spring vector */
      for (j = 0; j < 3; j++)
        s[j] += ub[j] - ua[j];

      /* periodic boundary conditions */
      ImposePBC(s, net->L);

      /* update */
      s2 += (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
    };

    return (s2 / net->nrOfSprings);
  }

  /* average ratio <s2/len> */
  double S2len(Network *net, double *u) {
    Integer i, a, b, len, j;
    double kappa, s[3], ua[3], ub[3], s2, s2len;

    /* assembly */
    for (s2len = 0., net->nrOfLoops = 0, i = 0; i < net->nrOfSprings; i++) {
      a = net->springs[i].a;
      b = net->springs[i].b;
      len = net->springs[i].len;
      kappa = 3. / len;

      /* initial spring vector */
      s[0] = net->nodes[b].x - net->nodes[a].x;
      s[1] = net->nodes[b].y - net->nodes[a].y;
      s[2] = net->nodes[b].z - net->nodes[a].z;

      /* spring displacement vectors */
      for (j = 0; j < 3; j++) {
        ua[j] = u[3 * a + j];
        ub[j] = u[3 * b + j];
      }

      /* current spring vector */
      for (j = 0; j < 3; j++)
        s[j] += ub[j] - ua[j];

      /* periodic boundary conditions */
      ImposePBC(s, net->L);

      /* update */
      s2 = s[0] * s[0] + s[1] * s[1] + s[2] * s[2];
      s2len += s2 / len;

      /* loop count */
      net->springs[i].isActive = 1;
      if (s2 < EPS2) {
        net->nrOfLoops++;
        net->springs[i].isActive = 0;
      }
    };
    return (s2len / net->nrOfSprings);
  }

private:
  pylimer_tools::entities::Universe universe;
};
} // namespace mehp
} // namespace calc
} // namespace pylimer_tools
#endif
