#ifndef BOX_H
#define BOX_H

namespace pylimer_tools
{
  namespace entities
  {
    class Box
    {
    public:
      Box(const double Lx = 0.0, const double Ly = 0.0, const double Lz = 0.0)
      {
        this->Lx = Lx;
        this->Ly = Ly;
        this->Lz = Lz;
      }

      double getVolume()
      {
        return this->Lx * this->Ly * this->Lz;
      }

      inline double getLx() { return this->Lx; }
      inline double getLy() { return this->Ly; }
      inline double getLz() { return this->Lz; }

      double Lx, Ly, Lz;
    };
  }
}

#endif
