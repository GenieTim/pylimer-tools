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

#include "Correlator.h"
#include <algorithm>
#include <math.h>
#include <cstring> // for memcpy

namespace pylimer_tools {
namespace calc {
  /////////////////////////////////////////
  // Correlator class
  /////////////////////////////////////////
  Correlator::Correlator(const unsigned int numcorrin,
                         const unsigned int pin,
                         const unsigned int min)
  {
    setsize(numcorrin, pin, min);
  }

  Correlator::~Correlator()
  {

    if (numcorrelators == 0)
      return;

    delete[] shift;
    delete[] correlation;
    delete[] ncorrelation;
    delete[] accumulator;
    delete[] naccumulator;
    delete[] insertindex;

    delete[] t;
    delete[] f;
  }

  // copy-constructor
  Correlator::Correlator(const Correlator& src)
    : Correlator(src.numcorrelators, src.p, src.m)
  {
    // copy values/etc.
    std::memcpy(this->shift, src.shift, sizeof(this->shift));
    std::memcpy(this->correlation, src.correlation, sizeof(this->correlation));
    std::memcpy(
      this->ncorrelation, src.ncorrelation, sizeof(this->ncorrelation));
    std::memcpy(this->accumulator, src.accumulator, sizeof(this->accumulator));
    std::memcpy(
      this->naccumulator, src.naccumulator, sizeof(this->naccumulator));
    std::memcpy(this->insertindex, src.insertindex, sizeof(this->insertindex));
    std::memcpy(this->t, src.t, sizeof(this->t));
    std::memcpy(this->f, src.f, sizeof(this->f));
  }

  // copy-assignment operator
  Correlator& Correlator::operator=(Correlator src)
  {
    if (this == &src) {
      return *this;
    }

    std::swap(this->shift, src.shift);
    std::swap(this->correlation, src.correlation);
    std::swap(this->ncorrelation, src.ncorrelation);
    std::swap(this->accumulator, src.accumulator);
    std::swap(this->naccumulator, src.naccumulator);
    std::swap(this->insertindex, src.insertindex);
    std::swap(this->numcorrelators, src.numcorrelators);
    std::swap(this->dmin, src.dmin);
    std::swap(this->length, src.length);
    std::swap(this->kmax, src.kmax);
    std::swap(this->p, src.p);
    std::swap(this->m, src.m);
    std::swap(this->t, src.t);
    std::swap(this->f, src.f);
    std::swap(this->npcorr, src.npcorr);
    std::swap(this->accval, src.accval);

    return *this;
  }

  void Correlator::setsize(const unsigned int numcorrin,
                           const unsigned int pin,
                           const unsigned int min)
  {
    numcorrelators = numcorrin;
    p = pin;
    m = min;
    dmin = p / m;

    length = numcorrelators * p;

    shift = new double*[numcorrelators];
    correlation = new double*[numcorrelators];
    ncorrelation = new unsigned long int*[numcorrelators];
    accumulator = new double[numcorrelators];
    naccumulator = new unsigned int[numcorrelators];
    insertindex = new unsigned int[numcorrelators];

    for (unsigned int j = 0; j < numcorrelators; ++j) {
      shift[j] = new double[p];

      /* It can be optimized: Apart from correlator 0, correlation and
       * ncorrelation arrays only use p/2 values */
      correlation[j] = new double[p];
      ncorrelation[j] = new unsigned long int[p];
    }

    t = new double[length];
    f = new double[length];

    this->initialize();
  }

  void Correlator::initialize()
  {
    for (unsigned int j = 0; j < numcorrelators; ++j) {
      for (unsigned int i = 0; i < p; ++i) {
        shift[j][i] = -2E10;
        correlation[j][i] = 0;
        ncorrelation[j][i] = 0;
      }
      accumulator[j] = 0.0;
      naccumulator[j] = 0;
      insertindex[j] = 0;
    }

    for (unsigned int i = 0; i < length; ++i) {
      t[i] = 0;
      f[i] = 0;
    }

    npcorr = 0;
    kmax = 0;
    accval = 0;
  }

  void Correlator::add(const double w, const unsigned int k)
  {

    /// If we exceed the correlator side, the value is discarded
    if (k == numcorrelators)
      return;
    if (k > kmax)
      kmax = k;

    /// Insert new value in shift array
    shift[k][insertindex[k]] = w;

    /// Add to average value
    if (k == 0)
      accval += w;

    /// Add to accumulator and, if needed, add to next correlator
    accumulator[k] += w;
    ++naccumulator[k];
    if (naccumulator[k] == m) {
      add(accumulator[k] / m, k + 1);
      accumulator[k] = 0;
      naccumulator[k] = 0;
    }

    /// Calculate correlation function
    unsigned int ind1 = insertindex[k];
    if (k == 0) { /// First correlator is different
      int ind2 = ind1;
      for (unsigned int j = 0; j < p; ++j) {
        if (shift[k][ind2] > -1e10) {
          correlation[k][j] += shift[k][ind1] * shift[k][ind2];
          ++ncorrelation[k][j];
        }
        --ind2;
        if (ind2 < 0)
          ind2 += p;
      }
    } else {
      int ind2 = ind1 - dmin;
      for (unsigned int j = dmin; j < p; ++j) {
        if (ind2 < 0)
          ind2 += p;
        if (shift[k][ind2] > -1e10) {
          correlation[k][j] += shift[k][ind1] * shift[k][ind2];
          ++ncorrelation[k][j];
        }
        --ind2;
      }
    }

    ++insertindex[k];
    if (insertindex[k] == p)
      insertindex[k] = 0;
  }

  void Correlator::evaluate(const bool norm)
  {
    unsigned int im = 0;

    double aux = 0;
    if (norm)
      aux = (accval / ncorrelation[0][0]) * (accval / ncorrelation[0][0]);

    // First correlator
    for (unsigned int i = 0; i < p; ++i) {
      if (ncorrelation[0][i] > 0) {
        t[im] = i;
        f[im] = correlation[0][i] / ncorrelation[0][i] - aux;
        ++im;
      }
    }

    // Subsequent correlators
    for (int k = 1; k < kmax; ++k) {
      for (int i = dmin; i < p; ++i) {
        if (ncorrelation[k][i] > 0) {
          t[im] = i * pow((double)m, k);
          f[im] = correlation[k][i] / ncorrelation[k][i] - aux;
          ++im;
        }
      }
    }

    npcorr = im;
  }

}
}
