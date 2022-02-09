#ifndef BOX_H
#define BOX_H

namespace pylimer_tools {
namespace entities {
// TODO: currently, this way, no tilt or more complicated boxes etc. is
// supported
class Box {
private:
  double Lx, Ly, Lz;
  double xLo, xHi, yLo, yHi, zLo, zHi;

public:
  Box(const double Lx = 0.0, const double Ly = 0.0, const double Lz = 0.0) {
    this->Lx = Lx;
    this->Ly = Ly;
    this->Lz = Lz;
    this->xLo = 0.0;
    this->xHi = Lx;
    this->yLo = 0.0;
    this->yHi = Ly;
    this->zLo = 0.0;
    this->zHi = Lz;
  }

  Box(const double xLo, const double xHi, const double yLo, const double yHi,
      const double zLo, const double zHi) {
    this->Lx = xHi - xLo;
    this->Ly = yHi - yLo;
    this->Lz = zHi - zLo;
    this->xLo = xLo;
    this->xHi = xHi;
    this->yLo = yLo;
    this->yHi = yHi;
    this->zLo = zLo;
    this->zHi = zHi;
  }

  double getVolume() const { return this->Lx * this->Ly * this->Lz; }

  double getLx() const { return this->Lx; }
  double getLy() const { return this->Ly; }
  double getLz() const { return this->Lz; }

  double getLowX() const { return this->xLo; }
  double getLowY() const { return this->yLo; }
  double getLowZ() const { return this->zLo; }

  double getHighX() const { return this->xHi; }
  double getHighY() const { return this->yHi; }
  double getHighZ() const { return this->zHi; }
};
} // namespace entities
} // namespace pylimer_tools

#endif
