#ifndef ATOM_H
#define ATOM_H

#include "Box.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
// #include <iostream>
#include <iterator>
#include <vector>

namespace pylimer_tools {
namespace entities {

  class Atom
  {
  public:
    Atom(const long int id,
         const int type,
         const double x,
         const double y,
         const double z,
         const int nx = 0,
         const int ny = 0,
         const int nz = 0)
    {
      this->id = id;
      this->type = type;
      this->x = x;
      this->y = y;
      this->z = z;
      this->nx = nx;
      this->ny = ny;
      this->nz = nz;
    };

    bool operator==(const Atom& ref) const
    {
      return this->id == ref.id && this->type == ref.type && this->x == ref.x &&
             this->y == ref.y && this->z == ref.z && this->nx == ref.nx &&
             this->ny == ref.ny && this->nz == ref.nz;
    }

    Eigen::Vector3d vectorTo(const Atom& b, const Box& box) const
    {
      Eigen::Vector3d dist = b.getCoordinates() - this->getCoordinates();
      box.handlePBC(dist);
      return dist;
    }

    Eigen::Vector3d vectorToUnwrapped(const Atom& b, const Box& box) const
    {
      return b.getUnwrappedCoordinates(box) -
             this->getUnwrappedCoordinates(box);
    }

    Eigen::Vector3d meanPositionWith(const Atom& b, const Box& box) const
    {
      Eigen::Vector3d result =
        this->getCoordinates() + 0.5 * this->vectorTo(b, box);
      box.handlePBC(result); // move into box
      return result;
    }

    Eigen::Vector3d meanPositionWithUnwrapped(const Atom& b,
                                              const Box& box) const
    {
      Eigen::Vector3d result =
        this->getCoordinates() + 0.5 * this->vectorToUnwrapped(b, box);
      box.handlePBC(result); // move into box
      return result;
    }

    double distanceTo(const Atom& b, const Box& box) const
    {
      return this->vectorTo(b, box).norm();
    }

    double distanceToUnwrapped(const Atom& b, const Box& box) const
    {
      return this->vectorToUnwrapped(b, box).norm();
    }

    long int getId() const { return this->id; }
    int getType() const { return this->type; }
    double getX() const { return this->x; }
    double getY() const { return this->y; }
    double getZ() const { return this->z; }
    double getUnwrappedX(const Box& box) const
    {
      return this->x + (this->nx * box.getLx());
    }
    double getUnwrappedY(const Box& box) const
    {
      return this->y + (this->ny * box.getLy());
    }
    double getUnwrappedZ(const Box& box) const
    {
      return this->z + (this->nz * box.getLz());
    }
    int getNX() const { return this->nx; }
    int getNY() const { return this->ny; }
    int getNZ() const { return this->nz; }

    template<typename VectorType>
    void getCoordinates(VectorType& vec) const
    {
      INVALIDARG_EXP_IFN(vec.size() == 3,
                         "Expect coordinates to be in order x, y, z, i.e., a "
                         "vector or array of size 3.");
      vec[0] = this->getX();
      vec[1] = this->getY();
      vec[2] = this->getZ();
    }
    Eigen::Vector3d getCoordinates() const
    {
      Eigen::Vector3d coords = Eigen::Vector3d::Zero();
      this->getCoordinates<Eigen::Vector3d>(coords);
      return coords;
    }
    template<typename VectorType>
    void getUnwrappedCoordinates(VectorType& vec, const Box& box) const
    {
      INVALIDARG_EXP_IFN(vec.size() == 3,
                         "Expect coordinates to be in order x, y, z, i.e., a "
                         "vector or array of size 3.");
      vec[0] = this->getUnwrappedX(box);
      vec[1] = this->getUnwrappedY(box);
      vec[2] = this->getUnwrappedZ(box);
    }
    Eigen::Vector3d getUnwrappedCoordinates(const Box& box) const
    {
      Eigen::Vector3d coords = Eigen::Vector3d::Zero();
      this->getUnwrappedCoordinates<Eigen::Vector3d>(coords, box);
      return coords;
    }

  private:
    long int id;
    int type;
    double x, y, z;
    int nx, ny, nz;
  };
} // namespace entities
} // namespace pylimer_tools
#endif
