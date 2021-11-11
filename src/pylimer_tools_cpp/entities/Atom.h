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

      double _getDeltaDistanceUnwrapped(double c1, double c2, int n1, int n2, double boxL) const
      {
        double delta = abs(c1 - c2);
        if (n1 != n2)
        {
          delta -= ((double)(n1 - n2)) * boxL;
        }
        return delta;
      }

      double _getDeltaDistance(double c1, double c2, int n1, int n2, double boxL) const
      {
        double delta = abs(c1 - c2);
        if (delta >= boxL * 0.5)
        {
          delta -= floor(delta / boxL) * boxL;
          // delta = abs(delta);
        }
        return delta;
      }

      const void vectorTo(Atom b, Box *box, double *result) const
      {
        result[0] = this->_getDeltaDistance(this->x, b.getX(), this->nx, b.getNX(), box->getLx());
        result[1] = this->_getDeltaDistance(this->y, b.getY(), this->ny, b.getNY(), box->getLy());
        result[2] = this->_getDeltaDistance(this->z, b.getZ(), this->nz, b.getNZ(), box->getLz());
      }

      const double distanceTo(Atom b, Box *box) const
      {
        double distanceVec[3];
        vectorTo(b, box, distanceVec);
        // norm
        return sqrt(distanceVec[0] * distanceVec[0] + distanceVec[1] * distanceVec[1] + distanceVec[2] * distanceVec[2]);
      }

      const void vectorToUnwrapped(Atom b, Box *box, double *result) const
      {
        result[0] = this->_getDeltaDistanceUnwrapped(this->x, b.getX(), this->nx, b.getNX(), box->getLx());
        result[1] = this->_getDeltaDistanceUnwrapped(this->y, b.getY(), this->ny, b.getNY(), box->getLy());
        result[2] = this->_getDeltaDistanceUnwrapped(this->z, b.getZ(), this->nz, b.getNZ(), box->getLz());
      }

      const double distanceToUnwrapped(Atom b, Box *box) const
      {
        double distanceVec[3];
        vectorToUnwrapped(b, box, distanceVec);
        // norm
        return sqrt(distanceVec[0] * distanceVec[0] + distanceVec[1] * distanceVec[1] + distanceVec[2] * distanceVec[2]);
      }

      const inline long int getId() const { return this->id; }
      const inline int getType() const { return this->type; }
      const inline double getX() const { return this->x; }
      const inline double getY() const { return this->y; }
      const inline double getZ() const { return this->z; }
      const inline double getUnwrappedX(Box *box) const { return this->x * this->nx * box->getLx(); }
      const inline double getUnwrappedY(Box *box) const { return this->y * this->ny * box->getLy(); }
      const inline double getUnwrappedZ(Box *box) const { return this->z * this->nz * box->getLz(); }
      const int getNX() const { return this->nx; }
      const int getNY() const { return this->ny; }
      const int getNZ() const { return this->nz; }

    private:
      long int id;
      int type;
      double x, y, z;
      int nx, ny, nz;
    };
  }
}
#endif
