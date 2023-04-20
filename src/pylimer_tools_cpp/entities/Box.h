#ifndef BOX_H
#define BOX_H

#include "../utils/utilityMacros.h"
#include <array>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>

namespace pylimer_tools {
namespace entities {
  // TODO: currently, this way, no tilt or more complicated boxes etc. is
  // supported
  class Box
  {
  private:
    double L[3];
    double boxHalfs[3];
    double oneOverL[3];
    double loCoords[3];
    double hiCoords[3];
    double simpleShearMagnitude = 0.0;
    int shearDirection = -1;

  protected:
    double minImageDistance(double dcoord, const int coord) const
    {
      return dcoord -
             (this->L[coord] * std::nearbyint(dcoord * this->oneOverL[coord]));
    }

    double placeIn(double dcoord, const int coord) const
    {
      return this->minImageDistance(dcoord, coord) + this->loCoords[coord];
    }

    double iterateForPlacementIn(double coord, double min, double max) const
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
            "Too many iterations in PBC, currently at " +
            std::to_string(coord) + " from " + std::to_string(coord0) +
            " in box with min/max " + std::to_string(min) + "/" +
            std::to_string(max) + " after " + std::to_string(min_iterations) +
            " iterations");
        }
      }
      int max_iterations = 0;
      while (coord < min) {
        coord += (max - min) / 2;
        max_iterations++;
        if (max_iterations > 50) {
          throw std::runtime_error(
            "Too many iterations in PBC, currently at " +
            std::to_string(coord) + " from " + std::to_string(coord0) +
            " in box with min/max " + std::to_string(min) + "/" +
            std::to_string(max) + " after " + std::to_string(max_iterations) +
            " iterations (and " + std::to_string(min_iterations) +
            " before that)");
        }
      }
      return coord;
    }

    void recomputeBoxProperties()
    {
      for (unsigned int i = 0; i < 3; ++i) {
        this->L[i] = this->hiCoords[i] - this->loCoords[i];
        this->boxHalfs[i] = L[i] * 0.5;
        this->oneOverL[i] = 1.0 / L[i];
      }
    }

  public:
    Box(const double Lx = 0.0, const double Ly = 0.0, const double Lz = 0.0)
    {
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
      this->loCoords[0] = xLo;
      this->hiCoords[0] = xHi;
      this->loCoords[1] = yLo;
      this->hiCoords[1] = yHi;
      this->loCoords[2] = zLo;
      this->hiCoords[2] = zHi;
      this->recomputeBoxProperties();
    }

    void applySimpleShear(const double shearMagnitude,
                          int newShearDirection = 0)
    {
      this->simpleShearMagnitude = shearMagnitude;
      this->shearDirection = newShearDirection;
    }

    double getVolume() const
    {
      return this->getLx() * this->getLy() * this->getLz();
    }

    double getLx() const { return this->L[0]; }
    double getLy() const { return this->L[1]; }
    double getLz() const { return this->L[2]; }

    double getLowX() const { return this->loCoords[0]; }
    double getLowY() const { return this->loCoords[1]; }
    double getLowZ() const { return this->loCoords[2]; }

    double getHighX() const { return this->hiCoords[0]; }
    double getHighY() const { return this->hiCoords[1]; }
    double getHighZ() const { return this->hiCoords[2]; }

    double getShearMagnitude() const { return this->simpleShearMagnitude; }
    int getShearDirection() const { return this->shearDirection; }

    template<typename VectorType>
    VectorType placeInBox(VectorType& coords) const
    {
      INVALIDARG_EXP_IFN(
        coords.size() % 3 == 0,
        "Expect coordinates to be in order x, y, z, repeatedly.");

      for (size_t i = 0; i < coords.size(); ++i) {
        coords[i] = this->placeIn(coords[i], i % 3);
      }

      return coords;
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
      RUNTIME_EXP_IFN(
        scalingFactorX > 0.,
        "Requiring scaling factor to be > 0, got in x-direction " +
          std::to_string(scalingFactorX) + ".");
      RUNTIME_EXP_IFN(
        scalingFactorY > 0.,
        "Requiring scaling factor to be > 0, got in y-direction " +
          std::to_string(scalingFactorY) + ".");
      RUNTIME_EXP_IFN(
        scalingFactorZ > 0.,
        "Requiring scaling factor to be > 0, got in z-direction " +
          std::to_string(scalingFactorZ) + ".");

      // first, scale back to non-sheared.
      if (this->shearDirection >= 0 && this->shearDirection <= 3) {
        for (size_t i = 0; i < coords.size() / 3; ++i) {
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
      for (size_t i = 0; i < coords.size() / 3; ++i) {
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
    }

    template<typename VectorType>
    void handlePBC(VectorType& distances) const
    {
      const bool isSheared =
        (this->getShearDirection() >= 0 && this->getShearDirection() <= 2);
      if (isSheared) {
        INVALIDARG_EXP_IFN(distances.size() % 3 == 0,
                           "Require distances to be a multiple of 3 to handle "
                           "PBC for sheared box.");
        // scaled coordinates in the initial cubic box
        for (size_t j = 0; j < distances.size() / 3; ++j) {
          if (this->getShearDirection() == 0) {
            distances[3 * j] -=
              this->getShearMagnitude() * distances[3 * j + 1];
          }
          if (this->getShearDirection() == 1) {
            distances[3 * j + 1] -=
              this->getShearMagnitude() * distances[3 * j + 2];
          }
          if (this->getShearDirection() == 2) {
            distances[3 * j + 2] -=
              this->getShearMagnitude() * distances[3 * j];
          }
        }
      }
      // actually do PBC
      this->placeInBox(distances);
      // back to the physical space
      if (isSheared) {
        // scaled coordinates in the initial cubic box
        for (size_t j = 0; j < distances.size() / 3; ++j) {
          if (this->getShearDirection() == 0) {
            distances[3 * j] +=
              this->getShearMagnitude() * distances[3 * j + 1];
          }
          if (this->getShearDirection() == 1) {
            distances[3 * j + 1] +=
              this->getShearMagnitude() * distances[3 * j + 2];
          }
          if (this->getShearDirection() == 2) {
            distances[3 * j + 2] +=
              this->getShearMagnitude() * distances[3 * j];
          }
        }
      }
    }

    inline bool operator==(const Box& rhs)
    {
      const Box& lhs = *this;
      return lhs.getLowX() == rhs.getLowX() && lhs.getLowY() == rhs.getLowY() &&
             lhs.getLowZ() == rhs.getLowZ() &&
             lhs.getHighX() == rhs.getHighX() &&
             lhs.getHighY() == rhs.getHighY() &&
             lhs.getHighZ() == rhs.getHighZ() &&
             lhs.getShearDirection() == rhs.getShearDirection() &&
             lhs.getShearMagnitude() == rhs.getShearMagnitude();
    }
  };
} // namespace entities
} // namespace pylimer_tools

#endif
