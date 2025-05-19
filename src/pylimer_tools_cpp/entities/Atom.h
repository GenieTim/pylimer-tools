#pragma once

#include "../io/DataFileParser.h"
#include "../utils/VectorUtils.h"
#include "Box.h"

#include <Eigen/Dense>
#include <cmath>
#include <unordered_map>
// #include <iostream>

namespace pylimer_tools::entities {

class Atom
{
public:
  Atom(const long int idNow,
       const int typeNow,
       const double xNow,
       const double yNow,
       const double zNow,
       const int nxNow = 0,
       const int nyNow = 0,
       const int nzNow = 0)
    : id(idNow)
    , type(typeNow)
    , x(xNow)
    , y(yNow)
    , z(zNow)
    , nx(nxNow)
    , ny(nyNow)
    , nz(nzNow) {};

  explicit Atom(std::unordered_map<std::string, double>& properties)
    : Atom(static_cast<long int>(std::lround(properties["id"])),
           static_cast<int>(std::rint(properties["type"])),
           properties["x"],
           properties["y"],
           properties["z"],
           static_cast<int>(std::rint(properties["nx"])),
           static_cast<int>(std::rint(properties["ny"])),
           static_cast<int>(std::rint(properties["nz"])))
  {
    this->extraData = properties;
  };

  bool operator==(const Atom& ref) const
  {
    // compare properties we know that they exist
    const bool stdPropertiesEqual =
      (this->id == ref.id && this->type == ref.type &&
       APPROX_EQUAL(this->x, ref.x, 1e-12) &&
       APPROX_EQUAL(this->y, ref.y, 1e-12) &&
       APPROX_EQUAL(this->z, ref.z, 1e-12) && this->nx == ref.nx &&
       this->ny == ref.ny && this->nz == ref.nz);
    if (!stdPropertiesEqual) {
      return false;
    }
    // if we have a match so far, compare the extra data
    for (const auto& [key, value] : ref.extraData) {
      if (!pylimer_tools::utils::map_has_key(this->extraData, key)) {
        return false;
      }
      if (!APPROX_EQUAL(this->extraData.at(key), value, 1e-12)) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] Eigen::Vector3d vectorTo(const Atom& b, const Box& box) const
  {
    Eigen::Vector3d dist = b.getCoordinates() - this->getCoordinates();
    box.handlePBC(dist);
    return dist;
  }

  [[nodiscard]] Eigen::Vector3d vectorToUnwrapped(const Atom& b,
                                                  const Box& box) const
  {
    return b.getUnwrappedCoordinates(box) - this->getUnwrappedCoordinates(box);
  }

  [[nodiscard]] Eigen::Vector3d meanPositionWith(const Atom& b,
                                                 const Box& box) const
  {
    Eigen::Vector3d result =
      this->getCoordinates() + 0.5 * this->vectorTo(b, box);
    box.handlePBC(result); // move into box
    return result;
  }

  [[nodiscard]] Eigen::Vector3d meanPositionWithUnwrapped(const Atom& b,
                                                          const Box& box) const
  {
    Eigen::Vector3d result =
      this->getCoordinates() + 0.5 * this->vectorToUnwrapped(b, box);
    box.handlePBC(result); // move into box
    return result;
  }

  [[nodiscard]] double distanceTo(const Atom& b, const Box& box) const
  {
    return this->vectorTo(b, box).norm();
  }

  [[nodiscard]] double distanceToUnwrapped(const Atom& b, const Box& box) const
  {
    return this->vectorToUnwrapped(b, box).norm();
  }

  [[nodiscard]] long int getId() const { return this->id; }
  [[nodiscard]] int getType() const { return this->type; }
  [[nodiscard]] double getX() const { return this->x; }
  [[nodiscard]] double getY() const { return this->y; }
  [[nodiscard]] double getZ() const { return this->z; }
  [[nodiscard]] double getUnwrappedX(const Box& box) const
  {
    return this->x + (this->nx * box.getLx());
  }
  [[nodiscard]] double getUnwrappedY(const Box& box) const
  {
    return this->y + (this->ny * box.getLy());
  }
  [[nodiscard]] double getUnwrappedZ(const Box& box) const
  {
    return this->z + (this->nz * box.getLz());
  }
  [[nodiscard]] int getNX() const { return this->nx; }
  [[nodiscard]] int getNY() const { return this->ny; }
  [[nodiscard]] int getNZ() const { return this->nz; }

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
  [[nodiscard]] Eigen::Vector3d getCoordinates() const
  {
    return Eigen::Vector3d(this->x, this->y, this->z);
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
  [[nodiscard]] Eigen::Vector3d getUnwrappedCoordinates(const Box& box) const
  {
    Eigen::Vector3d coords = Eigen::Vector3d::Zero();
    this->getUnwrappedCoordinates<Eigen::Vector3d>(coords, box);
    return coords;
  }

  [[nodiscard]] std::unordered_map<std::string, double> getExtraData() const
  {
    return this->extraData;
  }

  [[nodiscard]] double getProperty(const std::string property) const
  {
    return this->extraData.at(property);
  }

private:
  long int id;
  int type;
  double x, y, z;
  int nx, ny, nz;
  std::unordered_map<std::string, double> extraData;
};
}