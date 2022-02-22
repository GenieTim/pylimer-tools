#ifndef ATOM_H
#define ATOM_H

#include "Box.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <vector>

namespace pylimer_tools {
namespace entities {

class Atom {
public:
  Atom(const long int id, const int type, const double x, const double y,
       const double z, const int nx = 0, const int ny = 0, const int nz = 0) {
    this->id = id;
    this->type = type;
    this->x = x;
    this->y = y;
    this->z = z;
    this->nx = nx;
    this->ny = ny;
    this->nz = nz;
  };

  double _getDeltaDistanceUnwrapped(double c1, double c2, int n1, int n2,
                                    double boxL) const {
    double delta = std::fabs(c1 - c2);
    if (n1 != n2) {
      delta -= ((double)(n1 - n2)) * boxL;
    }
    return delta;
  }

  double _getDeltaDistance(double c1, double c2, int n1, int n2,
                           double boxL) const {
    double delta = c1 - c2;
    assert(!std::isnan(boxL) && !std::isnan(delta) && !std::isinf(delta) && !std::isinf(boxL));
    while (delta > 0.5 * boxL) {
      delta -= boxL;
    }
    while (delta < -0.5 * boxL) {
      delta += boxL;
    }

    return delta;
  }

  void vectorTo(Atom b, const Box *box, double *result) const {
    result[0] = this->_getDeltaDistance(this->x, b.getX(), this->nx, b.getNX(),
                                        box->getLx());
    result[1] = this->_getDeltaDistance(this->y, b.getY(), this->ny, b.getNY(),
                                        box->getLy());
    result[2] = this->_getDeltaDistance(this->z, b.getZ(), this->nz, b.getNZ(),
                                        box->getLz());
  }

  std::vector<double> computeVectorTo(Atom b, const Box box) const {
    double result[3];
    vectorTo(b, &box, result);
    std::vector<double> resultV;
    std::copy(std::begin(result), std::end(result),
              std::back_inserter(resultV));
    return resultV;
  }

  double distanceTo(Atom b, const Box *box) const {
    double distanceVec[3];
    vectorTo(b, box, distanceVec);
    // norm
    return sqrt(distanceVec[0] * distanceVec[0] +
                distanceVec[1] * distanceVec[1] +
                distanceVec[2] * distanceVec[2]);
  }

  void vectorToUnwrapped(Atom b, const Box *box, double *result) const {
    result[0] = this->_getDeltaDistanceUnwrapped(this->x, b.getX(), this->nx,
                                                 b.getNX(), box->getLx());
    result[1] = this->_getDeltaDistanceUnwrapped(this->y, b.getY(), this->ny,
                                                 b.getNY(), box->getLy());
    result[2] = this->_getDeltaDistanceUnwrapped(this->z, b.getZ(), this->nz,
                                                 b.getNZ(), box->getLz());
  }

  double distanceToUnwrapped(Atom b, const Box *box) const {
    double distanceVec[3];
    vectorToUnwrapped(b, box, distanceVec);
    // norm
    return sqrt(distanceVec[0] * distanceVec[0] +
                distanceVec[1] * distanceVec[1] +
                distanceVec[2] * distanceVec[2]);
  }

  inline long int getId() const { return this->id; }
  inline int getType() const { return this->type; }
  inline double getX() const { return this->x; }
  inline double getY() const { return this->y; }
  inline double getZ() const { return this->z; }
  inline double getUnwrappedX(const Box *box) const {
    return this->x * this->nx * box->getLx();
  }
  inline double getUnwrappedY(const Box *box) const {
    return this->y * this->ny * box->getLy();
  }
  inline double getUnwrappedZ(const Box *box) const {
    return this->z * this->nz * box->getLz();
  }
  int getNX() const { return this->nx; }
  int getNY() const { return this->ny; }
  int getNZ() const { return this->nz; }

private:
  long int id;
  int type;
  double x, y, z;
  int nx, ny, nz;
};
} // namespace entities
} // namespace pylimer_tools
#endif
