#ifndef BOX_H
#define BOX_H

namespace pylimer_tools
{
  namespace entities
  {
    // TODO: currently, this way, no tilt or more complicated boxes etc. is supported
    class Box
    {
    private:
      double Lx, Ly, Lz;

    public:
      Box(const double Lx = 0.0, const double Ly = 0.0, const double Lz = 0.0)
      {
        this->Lx = Lx;
        this->Ly = Ly;
        this->Lz = Lz;
      }

      const double getVolume()
      {
        return this->Lx * this->Ly * this->Lz;
      }

      const double getLx() const { return this->Lx; }
      const double getLy() const { return this->Ly; }
      const double getLz() const { return this->Lz; }
    };
  }
}

#endif
