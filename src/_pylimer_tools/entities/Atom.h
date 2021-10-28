#ifndef ATOM_H
#define ATOM_H

#include "Box.h"
#include "math.h"

namespace pylimer_tools
{
  namespace entities
  {

    class Atom
    {
    public:
      Atom(int id, int type, double x, double y, double z, int nx, int ny, int nz)
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

      double _getDeltaDistance(double c1, double c2, int n1, int n2, double boxL)
      {
        double delta = abs(c1 - c2);
        if (n1 != n2)
        {
          delta -= ((double)(n1 - n2)) * boxL;
        }
        return delta;
      }

      void vectorTo(Atom b, Box *box, double *result)
      {
        result[0] = this->_getDeltaDistance(this->x, b.getX(), this->nx, b.getNX(), box->getLx());
        result[1] = this->_getDeltaDistance(this->y, b.getY(), this->ny, b.getNY(), box->getLy());
        result[2] = this->_getDeltaDistance(this->z, b.getZ(), this->nz, b.getNZ(), box->getLz());
      }

      double distanceTo(Atom b, Box *box)
      {
        double distanceVec[3];
        vectorTo(b, box, distanceVec);
        // norm
        return sqrt(distanceVec[0] * distanceVec[0] + distanceVec[1] * distanceVec[1] + distanceVec[2] * distanceVec[2]);
      };

      inline double getX() { return this->x; }
      inline double getY() { return this->y; }
      inline double getZ() { return this->z; }
      int getNX() { return this->nx; }
      int getNY() { return this->ny; }
      int getNZ() { return this->nz; }

    private:
      int id, type;
      double x, y, z;
      int nx, ny, nz;
    };
  }
}
#endif
