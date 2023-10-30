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
#ifndef __correlator_h
#define __correlator_h

#include <stdio.h>

namespace pylimer_tools {
namespace calc {
  ////////////////////////////////////////////////////
  /// Standard Scalar Correlator f(tau)=<A(t)A(t+tau)>
  class Correlator
  {

  protected:
    /** Where the coming values are stored */
    double** shift;
    /** Array containing the actual calculated correlation function */
    double** correlation;
    /** Number of values accumulated in cor */
    unsigned long int** ncorrelation;

    /** Accumulator in each correlator */
    double* accumulator;
    /** Index that controls accumulation in each correlator */
    unsigned int* naccumulator;
    /** Index pointing at the position at which the current value is inserted */
    unsigned int* insertindex;

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

    /** Lenght of result arrays */
    unsigned int length;
    /** Maximum correlator attained during simulation */
    unsigned int kmax;

  public:
    /** Points per correlator */
    unsigned int p;
    /** Number of points over which to average; RECOMMENDED: p mod m = 0 */
    unsigned int m;
    double *t, *f;
    unsigned int npcorr;

    /** Accumulated result of incoming values **/
    double accval;

    /** Constructor */
    // Correlator() { numcorrelators = 0; };
    Correlator(const unsigned int numcorrin = 32,
               const unsigned int pin = 16,
               const unsigned int min = 2);
    ~Correlator();

    /** Set size of correlator */
    void setsize(const unsigned int numcorrin = 32,
                 const unsigned int pin = 16,
                 const unsigned int min = 2);

    /** Add a scalar to the correlator number k */
    void add(const double w, const unsigned int k = 0);

    /** Evaluate the current state of the correlator */
    void evaluate(const bool norm = false);

    /** Initialize all values (current and average) to zero */
    void initialize();
  };

}
}
#endif
