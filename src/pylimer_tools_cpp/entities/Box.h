#pragma once

#include "../utils/utilityMacros.h"
#include <Eigen/Dense>
#include <array>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
// #include <iostream>

namespace pylimer_tools::entities {
// TODO: currently, this way, no tilt or more complicated boxes etc. is
// supported
class Box
{
private:
  Eigen::Array3d L = Eigen::Array3d::Zero();
  double volume = 0.;
  Eigen::Array3d boxHalfs = Eigen::Array3d::Zero();
  Eigen::Array3d oneOverL = Eigen::Array3d::Zero();
  Eigen::Array3d loCoords = Eigen::Array3d::Zero();
  Eigen::Array3d hiCoords = Eigen::Array3d::Zero();
  double simpleShearMagnitude = 0.0;
  int shearDirection = -1;

protected:
  [[nodiscard]] inline double minImageDistance(const double dcoord,
                                               const int coord) const
  {
    return dcoord -
           (this->L[coord] * std::nearbyint(dcoord * this->oneOverL[coord]));
  }

  // Eigen::Vector3d minImageDistance(Eigen::Vector3d coords) const
  // {
  //   return coords -
  //          (this->L * (coords.array() * this->oneOverL).round()).matrix();
  // }
  template<typename Derived>
  inline void minImageDistance(Eigen::MatrixBase<Derived>& coords) const
  {
    coords -=
      (this->L.replicate(coords.size() / 3, 1) *
       (coords.array() * this->oneOverL.replicate(coords.size() / 3, 1)).rint())
        .matrix();
  }

  template<typename VectorType>
  inline VectorType minImageDistancesStdContainer(VectorType& coords) const
  {
    INVALIDARG_EXP_IFN(
      coords.size() % 3 == 0,
      "Expect coordinates to be in order x, y, z, repeatedly.");

    for (size_t i = 0; i < coords.size(); ++i) {
      coords[i] = this->minImageDistance(coords[i], i % 3);
    }

    return coords;
  }

  static double iterateForPlacementIn(double coord,
                                      const double min,
                                      const double max)
  {
    int min_iterations = 0;
    assert(!std::isinf(coord) && !std::isnan(coord));
    double coord0 = coord;
    // assert(max > min);
    while (coord > max) {
      coord -= (max - min) / 2;
      min_iterations++;
      if (min_iterations > 50) {
        throw std::runtime_error(
          "Too many iterations in PBC, currently at " + std::to_string(coord) +
          " from " + std::to_string(coord0) + " in box with min/max " +
          std::to_string(min) + "/" + std::to_string(max) + " after " +
          std::to_string(min_iterations) + " iterations");
      }
    }
    int max_iterations = 0;
    while (coord < min) {
      coord += (max - min) / 2;
      max_iterations++;
      if (max_iterations > 50) {
        throw std::runtime_error(
          "Too many iterations in PBC, currently at " + std::to_string(coord) +
          " from " + std::to_string(coord0) + " in box with min/max " +
          std::to_string(min) + "/" + std::to_string(max) + " after " +
          std::to_string(max_iterations) + " iterations (and " +
          std::to_string(min_iterations) + " before that)");
      }
    }
    return coord;
  }

  void recomputeBoxProperties()
  {
    this->L = this->hiCoords - this->loCoords;
    this->boxHalfs = 0.5 * this->L;
    this->oneOverL = 1.0 / this->L;
    this->volume = this->L.prod();
  }

public:
  explicit Box(const double Lx = 1.0,
               const double Ly = 1.0,
               const double Lz = 1.0)
  {
    INVALIDARG_EXP_IFN(Lx > 0 && Ly > 0 && Lz > 0,
                       "The box must be instantiated with positive lengths.");
    this->loCoords[0] = 0.0;
    this->hiCoords[0] = Lx;
    this->loCoords[1] = 0.0;
    this->hiCoords[1] = Ly;
    this->loCoords[2] = 0.0;
    this->hiCoords[2] = Lz;
    this->recomputeBoxProperties();
  }

  Box(const double xLo,
      const double xHi,
      const double yLo,
      const double yHi,
      const double zLo,
      const double zHi)
  {
    INVALIDARG_EXP_IFN(xHi > xLo && yHi > yLo && zHi > zLo,
                       "The box must be instantiated with positive lengths.");
    this->loCoords[0] = xLo;
    this->hiCoords[0] = xHi;
    this->loCoords[1] = yLo;
    this->hiCoords[1] = yHi;
    this->loCoords[2] = zLo;
    this->hiCoords[2] = zHi;
    this->recomputeBoxProperties();
  }

  Box(const Box& other)
  {
    this->loCoords = other.loCoords;
    this->hiCoords = other.hiCoords;
    this->shearDirection = other.shearDirection;
    this->simpleShearMagnitude = other.simpleShearMagnitude;
    this->recomputeBoxProperties();
  }

  Box& operator=(const Box& other)
  {
    if (this != &other) {
      this->loCoords = other.loCoords;
      this->hiCoords = other.hiCoords;
      this->shearDirection = other.shearDirection;
      this->simpleShearMagnitude = other.simpleShearMagnitude;
      this->recomputeBoxProperties();
    }
    return *this;
  }

  inline bool operator==(const Box& rhs) const
  {
    const Box& lhs = *this;
    return lhs.getLowX() == rhs.getLowX() && lhs.getLowY() == rhs.getLowY() &&
           lhs.getLowZ() == rhs.getLowZ() && lhs.getHighX() == rhs.getHighX() &&
           lhs.getHighY() == rhs.getHighY() &&
           lhs.getHighZ() == rhs.getHighZ() &&
           lhs.getShearDirection() == rhs.getShearDirection() &&
           lhs.getShearMagnitude() == rhs.getShearMagnitude();
  }

  void applySimpleShear(const double shearMagnitude,
                        const int newShearDirection = 0)
  {
    this->simpleShearMagnitude = shearMagnitude;
    this->shearDirection = newShearDirection;
  }

  [[nodiscard]] double getVolume() const { return this->volume; }

  [[nodiscard]] double getL(const int i) const { return this->L[i % 3]; }
  [[nodiscard]] Eigen::Array3d getL() const { return this->L; }
  [[nodiscard]] double getLowL(const int i) const
  {
    return this->loCoords[i % 3];
  }
  [[nodiscard]] Eigen::Array3d getLowL() const { return this->loCoords; }
  [[nodiscard]] double getHighL(const int i) const
  {
    return this->hiCoords[i % 3];
  }
  [[nodiscard]] Eigen::Array3d getHighL() const { return this->hiCoords; }

  [[nodiscard]] double getLx() const { return this->L[0]; }
  [[nodiscard]] double getLy() const { return this->L[1]; }
  [[nodiscard]] double getLz() const { return this->L[2]; }

  [[nodiscard]] double getLowX() const { return this->loCoords[0]; }
  [[nodiscard]] double getLowY() const { return this->loCoords[1]; }
  [[nodiscard]] double getLowZ() const { return this->loCoords[2]; }

  [[nodiscard]] double getHighX() const { return this->hiCoords[0]; }
  [[nodiscard]] double getHighY() const { return this->hiCoords[1]; }
  [[nodiscard]] double getHighZ() const { return this->hiCoords[2]; }

  [[nodiscard]] double getShearMagnitude() const
  {
    return this->simpleShearMagnitude;
  }
  [[nodiscard]] int getShearDirection() const { return this->shearDirection; }
  [[nodiscard]] bool isSheared() const
  {
    return this->getShearMagnitude() != 0 &&
           (this->getShearDirection() >= 0 && this->getShearDirection() <= 2);
    ;
  }

  template<typename T>
  inline std::vector<T> minImageDistances(std::vector<T>& coords) const
  {
    return this->minImageDistancesStdContainer(coords);
  }

  template<typename T, size_t N>
  inline std::array<T, N> minImageDistances(std::array<T, N>& coords) const
  {
    return this->minImageDistancesStdContainer(coords);
  }

  template<typename Derived>
  inline void minImageDistances(Eigen::MatrixBase<Derived>& coords) const
  {
    INVALIDARG_EXP_IFN(
      coords.size() % 3 == 0,
      "Expect coordinates to be in order x, y, z, repeatedly.");

    this->minImageDistance(coords);
  }

  /**
   * @brief Compute the offset required to compensate for PBC
   *
   * @param distanceVec
   * @return Eigen::VectorXd
   */
  [[nodiscard]] Eigen::VectorXd getOffset(
    const Eigen::VectorXd& distanceVec) const
  {
    return -(this->L.replicate(distanceVec.size() / 3, 1) *
             (distanceVec.array() *
              this->oneOverL.replicate(distanceVec.size() / 3, 1))
               .rint())
              .matrix();
  }

  [[nodiscard]] bool isValidOffset(const Eigen::VectorXd& offset,
                                   const double precision = 1e-5) const
  {
    if (!(offset.size() % 3 == 0)) {
      return false;
    }
    for (int i = 0; i < offset.size(); ++i) {
      if (offset[i] > precision) {
        double multiple = offset[i] / this->L[i % 3];
        bool rowIsValid = std::abs(multiple - round(multiple)) < precision;
        if (!rowIsValid) {
          return false;
        }
      }
    }
    return true;
  }

  template<typename VectorType>
  void adjustCoordinatesTo(VectorType& coords, const Box& newBox) const
  {
    INVALIDARG_EXP_IFN(
      coords.size() % 3 == 0,
      "Expect coordinates to be in order x, y, z, repeatedly.");

    double scalingFactorX = newBox.getLx() / this->L[0];
    double scalingFactorY = newBox.getLy() / this->L[1];
    double scalingFactorZ = newBox.getLz() / this->L[2];
    RUNTIME_EXP_IFN(scalingFactorX > 0.,
                    "Requiring scaling factor to be > 0, got in x-direction " +
                      std::to_string(scalingFactorX) + ".");
    RUNTIME_EXP_IFN(scalingFactorY > 0.,
                    "Requiring scaling factor to be > 0, got in y-direction " +
                      std::to_string(scalingFactorY) + ".");
    RUNTIME_EXP_IFN(scalingFactorZ > 0.,
                    "Requiring scaling factor to be > 0, got in z-direction " +
                      std::to_string(scalingFactorZ) + ".");

    // first, scale back to non-sheared.
    if (this->shearDirection >= 0 && this->shearDirection <= 3) {
      for (int i = 0; i < coords.size() / 3; ++i) {
        if (this->getShearDirection() == 0) {
          coords[3 * i] -=
            this->getShearMagnitude() * coords[3 * i + 1]; // x' = x + ɣ*y
        }
        if (this->getShearDirection() == 1) {
          coords[3 * i + 1] -=
            this->getShearMagnitude() * coords[3 * i + 2]; // y' = y + ɣ*z
        }
        if (this->getShearDirection() == 2) {
          coords[3 * i + 2] -=
            this->getShearMagnitude() * coords[3 * i]; // z' = z + ɣ*x
        }
      }
    }

    // actually do the deformation as appropriate
    for (int i = 0; i < coords.size() / 3; ++i) {
      coords[3 * i] *= scalingFactorX;
      coords[3 * i + 1] *= scalingFactorY;
      coords[3 * i + 2] *= scalingFactorZ;
      if (newBox.getShearDirection() == 0) {
        coords[3 * i] +=
          newBox.getShearMagnitude() * coords[3 * i + 1]; // x' = x + ɣ*y
      }
      if (newBox.getShearDirection() == 1) {
        coords[3 * i + 1] +=
          newBox.getShearMagnitude() * coords[3 * i + 2]; // y' = y + ɣ*z
      }
      if (newBox.getShearDirection() == 2) {
        coords[3 * i + 2] +=
          newBox.getShearMagnitude() * coords[3 * i]; // z' = z + ɣ*x
      }
    }
  };

  template<typename VectorType>
  inline void handlePBC(VectorType& distances) const
  {
    const bool sheared = this->isSheared();
    if (sheared) {
      INVALIDARG_EXP_IFN(distances.size() % 3 == 0,
                         "Require distances to be a multiple of 3 to handle "
                         "PBC for sheared box.");
      // scaled coordinates in the initial cubic box
      for (int j = 0; j < distances.size() / 3; ++j) {
        if (this->getShearDirection() == 0) {
          distances[3 * j] -= this->getShearMagnitude() * distances[3 * j + 1];
        }
        if (this->getShearDirection() == 1) {
          distances[3 * j + 1] -=
            this->getShearMagnitude() * distances[3 * j + 2];
        }
        if (this->getShearDirection() == 2) {
          distances[3 * j + 2] -= this->getShearMagnitude() * distances[3 * j];
        }
      }
    }
    // actually do PBC
    this->minImageDistances(distances);
    // back to the physical space
    if (sheared) {
      // scaled coordinates in the initial cubic box
      for (int j = 0; j < distances.size() / 3; ++j) {
        if (this->getShearDirection() == 0) {
          distances[3 * j] += this->getShearMagnitude() * distances[3 * j + 1];
        }
        if (this->getShearDirection() == 1) {
          distances[3 * j + 1] +=
            this->getShearMagnitude() * distances[3 * j + 2];
        }
        if (this->getShearDirection() == 2) {
          distances[3 * j + 2] += this->getShearMagnitude() * distances[3 * j];
        }
      }
    }
  }

  /**
   * @brief Find a ("linear" in terms of volume) intermediate between two
   * different boxes
   *
   * NOTE: The implementation entails that the deformation might not be
   * symmetric.
   *
   * @param other
   * @param f interpolation factor
   * @return Box the new box
   */
  [[nodiscard]] Box interpolate(const Box& other, const double f) const
  {
    INVALIDARG_EXP_IFN(f >= 0. && f <= 1., "Cannot extrapolate box.");
    INVALIDARG_EXP_IFN(
      other.getShearDirection() == this->getShearDirection() ||
        (other.getShearDirection() == -1 || this->getShearDirection() == -1),
      "Cannot interpolate more than one shear direction");
    double newLx = (1 - f) * this->getLx() + f * other.getLx();
    double newLy = (this->getLx() * this->getLy() * (1 - f) +
                    other.getLx() * other.getLy() * f) /
                   (newLx);
    double newLz =
      (this->getVolume() * (1 - f) + other.getVolume() * f) / (newLy * newLx);

    Box newBox = Box(this->getLowX() * (newLx / (this->getLx())),
                     this->getHighX() * (newLx / (this->getLx())),
                     this->getLowY() * (newLy / (this->getLy())),
                     this->getHighY() * (newLy / (this->getLy())),
                     this->getLowZ() * (newLz / (this->getLz())),
                     this->getHighZ() * (newLz / (this->getLz())));
    newBox.applySimpleShear(
      (1 - f) * this->getShearMagnitude() + f * other.getShearMagnitude(),
      this->getShearDirection() == -1 ? other.getShearDirection()
                                      : this->getShearDirection());

    return newBox;
  }

  [[nodiscard]] Box getBoundingBox() const
  {
    if (this->isSheared()) {
      Eigen::Vector3d lowerCorner = this->getLowL();
      // this->handlePBC(lowerCorner);
      Eigen::Vector3d upperCorner = this->getHighL();
      upperCorner[this->getShearDirection()] *=
        (1. + this->getShearMagnitude());

      return { lowerCorner[0], upperCorner[0], lowerCorner[1],
               upperCorner[1], lowerCorner[2], upperCorner[2] };
    }
    return { *this };
  }

#ifdef CEREALIZABLE
  /**
   * @brief Serialize this box using Cereal
   *
   * @tparam Archive
   * @param ar
   */
  template<class Archive>
  void serialize(Archive& ar)
  {
    ar(L,
       volume,
       boxHalfs,
       oneOverL,
       loCoords,
       hiCoords,
       simpleShearMagnitude,
       shearDirection);
  }
#endif // CEREALIZABLE
};
}
