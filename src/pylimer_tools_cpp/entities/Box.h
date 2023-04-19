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
    double xLo, xHi, yLo, yHi, zLo, zHi;
    double simpleShearMagnitude = 0.0;
    int shearDirection = -1;

  protected:
    double iterateForPlacementIn(double coord, double min, double max) const
    {
      // assert(max > min);
      while (coord > max) {
        coord -= (max - min) / 2;
      }
      while (coord < min) {
        coord += (max - min) / 2;
      }
      return coord;
    }

    void recomputeBoxHalfs()
    {
      this->boxHalfs[0] = L[0] * 0.5;
      this->boxHalfs[1] = L[1] * 0.5;
      this->boxHalfs[2] = L[2] * 0.5;
    }

  public:
    Box(const double Lx = 0.0, const double Ly = 0.0, const double Lz = 0.0)
    {
      this->L[0] = Lx;
      this->L[1] = Ly;
      this->L[2] = Lz;
      this->recomputeBoxHalfs();
      this->xLo = 0.0;
      this->xHi = Lx;
      this->yLo = 0.0;
      this->yHi = Ly;
      this->zLo = 0.0;
      this->zHi = Lz;
    }

    Box(const double xLo,
        const double xHi,
        const double yLo,
        const double yHi,
        const double zLo,
        const double zHi)
    {
      this->L[0] = xHi - xLo;
      this->L[1] = yHi - yLo;
      this->L[2] = zHi - zLo;
      this->recomputeBoxHalfs();
      this->xLo = xLo;
      this->xHi = xHi;
      this->yLo = yLo;
      this->yHi = yHi;
      this->zLo = zLo;
      this->zHi = zHi;
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

    double getLowX() const { return this->xLo; }
    double getLowY() const { return this->yLo; }
    double getLowZ() const { return this->zLo; }

    double getHighX() const { return this->xHi; }
    double getHighY() const { return this->yHi; }
    double getHighZ() const { return this->zHi; }

    double getShearMagnitude() const { return this->simpleShearMagnitude; }
    int getShearDirection() const { return this->shearDirection; }

    template<typename VectorType>
    VectorType placeInBox(VectorType coords) const
    {
      INVALIDARG_EXP_IFN(
        coords.size() % 3 == 0,
        "Expect coordinates to be in order x, y, z, repeatedly.");

      for (size_t i = 0; i < coords.size(); i += 3) {
        coords[i + 0] = this->iterateForPlacementIn(
          coords[i + 0], this->getLowX(), this->getHighX());
        coords[i + 1] = this->iterateForPlacementIn(
          coords[i + 1], this->getLowY(), this->getHighY());
        coords[i + 2] = this->iterateForPlacementIn(
          coords[i + 2], this->getLowZ(), this->getHighZ());
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
      // possibly improveable PBC
      for (size_t j = 0; j < distances.size(); ++j) {
        int min_iterations = 0;
        int j_mod_3 = j % 3;
        assert(!std::isinf(distances[j]) && !std::isnan(distances[j]));
        const double distance0 = distances[j];
        while (distances[j] > this->boxHalfs[j_mod_3]) {
          distances[j] -= this->L[j_mod_3];
          min_iterations++;
          if (min_iterations > 50) {
            throw std::runtime_error(
              "Too many iterations in PBC at distance index " +
              std::to_string(j) + ", currently at " +
              std::to_string(distances[j]) + " from " +
              std::to_string(distance0) + " in box with halfs " +
              std::to_string(this->boxHalfs[j_mod_3]) + " after " +
              std::to_string(min_iterations) + " iterations");
          }
        }
        int max_iterations = 0;
        while (distances[j] < -this->boxHalfs[j_mod_3]) {
          distances[j] += this->L[j_mod_3];
          max_iterations++;
          if (max_iterations > 50) {
            throw std::runtime_error(
              "Too many iterations in PBC at distance index " +
              std::to_string(j) + ", currently at " +
              std::to_string(distances[j]) + " from " +
              std::to_string(distance0) + " in box with halfs " +
              std::to_string(this->boxHalfs[j_mod_3]) + " after " +
              std::to_string(max_iterations) + " iterations (and " +
              std::to_string(min_iterations) + " before that)");
          }
        }
      }
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
  };
} // namespace entities
} // namespace pylimer_tools

#endif
