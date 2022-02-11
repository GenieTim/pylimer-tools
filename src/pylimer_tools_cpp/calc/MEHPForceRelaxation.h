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
#define N0 100000     /* total number of strands as defined in pdms.m */

typedef int Integer; /* large enough integer */

typedef struct _Spring {
  Integer a;   /* first node */
  Integer b;   /* second node */
  Integer len; /* length */
  Integer act; /* 1/0; active or not */
} Spring;

typedef struct _Node {
  double x;    /* x-coordinate */
  double y;    /* y-coordinate */
  double z;    /* z-coordinate */
  Integer act; /* number of active springs connected to the node */
} Node;

typedef struct _Network {
  double L[3];     /* box sizes */
  double vol;      /* box volume */
  Integer nodes;   /* number of nodes */
  Integer springs; /* number of springs */
  Node *nod;       /* nodes */
  Spring *spr;     /* spings */
  double *u;       /* nodal displacements */
  double *force;   /* nodal forces */
  double avlen;    /* average spring length */
  Integer loops;   /* loops */
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

    for (int i = 0; i < net.springs; i++) {

      if (std::fabs(net.spr[i].len) < 1.0e-10) {
        std::cout << "WARNING: " << i << std::endl;
      }
      assert(net.spr[i].len != 0.0);
    }

    /* array allocation */

    u = (double *)calloc(3 * net.nodes, sizeof(double));
    r = (double *)calloc(3 * net.nodes, sizeof(double));

    /* initial */

    for (i = 0; i < 3 * net.nodes; i++) {
      u[i] = 0.0;
    }
    Stress(&net, u, stress);
    printf("\nInitial stress tensor:\n");
    for (i = 0; i < 3; i++) {
      printf("%.10f %.10f %.10f\n", stress[i][0], stress[i][1], stress[i][2]);
    }
    Residual(&net, u, r, &Fdef);
    for (r20 = 0., i = 0; i < 3 * net.nodes; i++) {
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
      for (i = 0; i < 3 * net.nodes; i++) {
        u[i] -= EPSILON * r[i];
      }
      Residual(&net, u, r, &Fdef);
      for (r2 = 0., i = 0; i < 3 * net.nodes; i++) {
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
      Gamma_eq = r2 / ((double)net.springs * Nb2);
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
    double R2_mean = Gamma_eq / (double)net.springs;
    double Gamma_x = 3 * Rx2_sum / ((double)net.springs * Nb2);
    double Gamma_y = 3 * Ry2_sum / ((double)net.springs * Nb2);
    double Gamma_z = 3 * Rz2_sum / ((double)net.springs * Nb2);

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
      fprintf(fp, "%d #nodes\n", net.nodes);
      fprintf(fp, "%d #springs\n", net.springs);
      for (i = 0; i < net.nodes; i++) {
        fprintf(fp, "%.16f %.16f %.16f\n", net.nod[i].x + u[3 * i],
                net.nod[i].y + u[3 * i + 1], net.nod[i].z + u[3 * i + 2]);
      }
      fclose(fp);
    }

    s2 = S2(&net, u);
    s2len = S2len(&net, u);

    fp = fopen("out.dat", "w");
    fprintf(fp, "%.10f %.10f %.10f %.10f   %.10f %.10f\n", N0 / net.vol,
            net.avlen, G0, G, s20, s2);
    fclose(fp);

    fp = fopen("s2len.dat", "w");
    fprintf(fp, "%.10f  %.10f %.10f  %d\n", N0 / net.vol, s20len, s2len,
            net.loops);
    fclose(fp);

    /* active springs */
    for (Nact = 0, i = 0; i < net.springs; i++) {
      if (net.spr[i].act == 1) {
        ++Nact;
      }
    }
    printf("\nactive springs: %d\n", Nact);

    /* active nodes */
    for (i = 0; i < net.nodes; i++)
      net.nod[i].act = 0; /* initial */
    for (i = 0; i < net.springs; i++) {
      if (net.spr[i].act == 1) /* active spring */
      {
        a = net.spr[i].a;
        b = net.spr[i].b;
        ++(net.nod[a].act);
        ++(net.nod[b].act);
      }
    }
    for (Mact = 0, i = 0; i < net.nodes; i++) {
      if (net.nod[i].act >= 2) {
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
    fprintf(fp, "%.10f   %.10f %.10f\n", N0 / net.vol, 4.14195 * G,
            4.14195 * s2len * N0 / net.vol);
    fclose(fp);

    /* exact and ANM shear moduli [MPa] using total and active strand densities
     */
    fp = fopen("ANM.dat", "w");
    fprintf(fp, "%.10f   %.10f %.10f %.10f\n", N0 / net.vol, 4.14195 * G,
            4.14195 * N0 / net.vol, 4.14195 * Nact / net.vol);
    fclose(fp);

    /* exact and PNM shear moduli [MPa] using total and active strand and node
     * densities */
    fp = fopen("PNM.dat", "w");
    fprintf(fp, "%.10f   %.10f %.10f %.10f\n", N0 / net.vol, 4.14195 * G,
            4.14195 * (N0 - net.nodes) / net.vol,
            4.14195 * (Nact - Mact) / net.vol);
    fclose(fp);

    /* total and active node and strand number densities */
    fp = fopen("Densities.dat", "w");
    fprintf(fp, "%.10f %.10f     %.10f %.10f\n", N0 / net.vol, Nact / net.vol,
            net.nodes / net.vol, Mact / net.vol);
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
    net->nodes = crosslinkerUniverse.getNrOfAtoms();
    net->springs = crosslinkerUniverse.getNrOfBonds();

    int usualChainLen =
        this->universe.getMolecules(crosslinkerType)[0].getNrOfAtoms();

    net->nod = (Node *)calloc(net->nodes, sizeof(Node));
    net->spr = (Spring *)calloc(net->springs, sizeof(Spring));

    // convert beads
    std::vector<pylimer_tools::entities::Atom> allAtoms =
        crosslinkerUniverse.getAtoms();
    std::map<int, int> atomIdToNode;
    for (int i = 0; i < allAtoms.size(); ++i) {
      pylimer_tools::entities::Atom atom = allAtoms[i];
      net->nod[i].x = atom.getX();
      net->nod[i].y = atom.getY();
      net->nod[i].z = atom.getZ();
      atomIdToNode[atom.getId()] = i;
    }

    // convert springs
    std::map<std::string, std::vector<long int>> allBonds =
        crosslinkerUniverse.getBonds();
    net->avlen = 0;
    for (int i = 0; i < net->springs; ++i) {
      int atomIdFrom = allBonds["bond_from"][i];
      int atomIdTo = allBonds["bond_to"][i];
      net->spr[i].a = atomIdToNode.at(atomIdFrom);
      net->spr[i].b = atomIdToNode.at(atomIdTo);
      net->spr[i].len = usualChainLen;
      net->avlen += usualChainLen;
    }

    assert(crosslinkerUniverse.getNrOfBonds() == net->springs);

    net->avlen /= net->springs;

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
    for (i = 0; i < 3 * net->nodes; i++)
      r[i] = 0.;

    /* assembly */
    for (i = 0; i < net->springs; i++) {
      a = net->spr[i].a;
      b = net->spr[i].b;
      len = net->spr[i].len;
      kappa = 3. / len;

      /* initial spring vector */
      s[0] = net->nod[b].x - net->nod[a].x;
      s[1] = net->nod[b].y - net->nod[a].y;
      s[2] = net->nod[b].z - net->nod[a].z;

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
    for (i = 0; i < net->springs; i++) {
      a = net->spr[i].a;
      b = net->spr[i].b;
      len = net->spr[i].len;
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
      s[0] = net->nod[b].x - net->nod[a].x;
      s[1] = net->nod[b].y - net->nod[a].y;
      s[2] = net->nod[b].z - net->nod[a].z;

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
    for (s2 = 0., i = 0; i < net->springs; i++) {
      a = net->spr[i].a;
      b = net->spr[i].b;
      len = net->spr[i].len;
      kappa = 3. / len;

      /* initial spring vector */
      s[0] = net->nod[b].x - net->nod[a].x;
      s[1] = net->nod[b].y - net->nod[a].y;
      s[2] = net->nod[b].z - net->nod[a].z;

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

    return (s2 / net->springs);
  }

  /* average ratio <s2/len> */
  double S2len(Network *net, double *u) {
    Integer i, a, b, len, j;
    double kappa, s[3], ua[3], ub[3], s2, s2len;

    /* assembly */
    for (s2len = 0., net->loops = 0, i = 0; i < net->springs; i++) {
      a = net->spr[i].a;
      b = net->spr[i].b;
      len = net->spr[i].len;
      kappa = 3. / len;

      /* initial spring vector */
      s[0] = net->nod[b].x - net->nod[a].x;
      s[1] = net->nod[b].y - net->nod[a].y;
      s[2] = net->nod[b].z - net->nod[a].z;

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
      net->spr[i].act = 1;
      if (s2 < EPS2) {
        net->loops++;
        net->spr[i].act = 0;
      }
    };
    /* For numerical efficiency, all dangling and free chains were discarded in
    pdms.m. They all have s2 = 0 and hence, their presence can now simply be
    exactly taken into account by normalizing by N0 instead of net->springs */
    return (s2len / N0);
  }

private:
  pylimer_tools::entities::Universe universe;
};
} // namespace mehp
} // namespace calc
} // namespace pylimer_tools
#endif
