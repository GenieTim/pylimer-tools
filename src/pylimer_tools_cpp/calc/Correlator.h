/**
Copyright (c) 2010 Jorge Ramirez

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */
#ifndef CORRELATOR_H
#define CORRELATOR_H

#include "../utils/CerealUtils.h"
#include <Eigen/Dense>
#include <string>

namespace pylimer_tools::calc {
typedef Eigen::Matrix<unsigned long int, Eigen::Dynamic, Eigen::Dynamic>
  MatrixXuli;
typedef Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> VectorXui;

////////////////////////////////////////////////////
/// Standard Scalar Correlator f(tau)=<A(t)A(t+tau)>
class Correlator
{

protected:
  /** Where the coming values are stored */
  Eigen::MatrixXd shift;
  /** Array containing the actual calculated correlation function */
  Eigen::MatrixXd correlation;
  /** Number of values accumulated in cor */
  MatrixXuli ncorrelation;

  /** Accumulator in each correlator */
  Eigen::VectorXd accumulator;
  /** Index that controls accumulation in each correlator */
  VectorXui naccumulator;
  /** Index pointing at the position at which the current value is inserted */
  VectorXui insertindex;

  /** Number of Correlators */
  unsigned int numcorrelators;

  /** Minimum distance between points for correlators k>0; dmin = p/m */
  unsigned int dmin;

  /*  SCHEMATIC VIEW OF EACH CORRELATOR
                                                  p=N
          <----------------------------------------------->
          _________________________________________________
          |0|1|2|3|.|.|.| | | | | | | | | | | | | | | |N-1|
          -------------------------------------------------
          */

  /** Length of result arrays */
  unsigned int length;
  /** Maximum correlator attained during simulation */
  unsigned int kmax;

public:
  /** Points per correlator */
  unsigned int p;
  /** Number of points over which to average; RECOMMENDED: p mod m = 0 */
  unsigned int m;
  Eigen::VectorXd t, f;
  unsigned int npcorr;

  /** Accumulated result of incoming values **/
  double accval;

  /** Constructor */
  Correlator(unsigned int numcorrin = 32,
             unsigned int pin = 16,
             unsigned int min = 2);

  /** Set size of correlator */
  void setsize(unsigned int numcorrin = 32,
               unsigned int pin = 16,
               unsigned int min = 2);

  /** Reset accumulated state while keeping current (numcorrelators, p, m) */
  void reset();

  /** Add a scalar to the correlator number k */
  void add(double w, unsigned int k = 0);

  /** Add a batch of scalar observations at level k = 0 */
  void addBatch(const Eigen::Ref<const Eigen::VectorXd>& values);

  /** Evaluate the current state of the correlator */
  void evaluate(bool norm = false);

  /** Return evaluated lag times (prefix of length npcorr) */
  Eigen::VectorXd getTimes() const;

  /** Return evaluated correlation values (prefix of length npcorr) */
  Eigen::VectorXd getCorrelation() const;

  /** Number of correlator levels configured */
  unsigned int getNumCorrelators() const { return numcorrelators; }

  /** Total number of samples added at level 0 */
  unsigned long int getNumSamples() const;

  /** Running sample mean at level 0 */
  double getMean() const;

#ifdef CEREALIZABLE
  static Correlator fromString(const std::string& in)
  {
    Correlator c;
    pylimer_tools::utils::deserializeFromString(c, in);
    return c;
  }
#endif

  /** serialize this object */
  template<class Archive>
  void serialize(Archive& ar)
  {
    ar(
      // private...
      shift,
      correlation,
      ncorrelation,
      accumulator,
      naccumulator,
      insertindex,
      numcorrelators,
      dmin,
      length,
      kmax,
      // ...and public members
      p,
      m,
      t,
      f,
      npcorr,
      accval);
  }
};

} // namespace pylimer_tools::calc
#endif
