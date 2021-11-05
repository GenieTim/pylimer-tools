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
      Atom(const long int id, const int type, const double x, const double y, const double z, const int nx = 0, const int ny = 0, const int nz = 0)
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

      inline long int getId() { return this->id; }
      inline int getType() { return this->type; }
      inline double getX() { return this->x; }
      inline double getY() { return this->y; }
      inline double getZ() { return this->z; }
      int getNX() { return this->nx; }
      int getNY() { return this->ny; }
      int getNZ() { return this->nz; }

    private:
      long int id;
      int type;
      double x, y, z;
      int nx, ny, nz;
    };
  }
}
#endif
