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
#include "../utils/utilityMacros.h"
#include <algorithm>
#include <cstring> // for memcpy
#include <limits.h>
#include <math.h>

namespace pylimer_tools::calc {
/////////////////////////////////////////
// Correlator class
/////////////////////////////////////////
Correlator::Correlator(const unsigned int numcorrin,
                       const unsigned int pin,
                       const unsigned int min)
{
  setsize(numcorrin, pin, min);
}

void
Correlator::setsize(const unsigned int numcorrin,
                    const unsigned int pin,
                    const unsigned int min)
{
  INVALIDARG_EXP_IFN(numcorrin < INT_MAX && pin < INT_MAX && min < INT_MAX,
                     "Arguments must be able to cast to integers.");
  numcorrelators = numcorrin;
  p = pin;
  m = min;
  dmin = p / m;

  length = numcorrelators * p;

  // initialize all values, set to 0
  shift = Eigen::MatrixXd::Constant(numcorrelators, p, -2e10);
  correlation = Eigen::MatrixXd::Zero(numcorrelators, p);
  ncorrelation = MatrixXuli::Zero(numcorrelators, p);
  accumulator = Eigen::VectorXd::Zero(numcorrelators);
  naccumulator = VectorXui::Zero(numcorrelators);
  insertindex = VectorXui::Zero(numcorrelators);

  t = Eigen::VectorXd::Zero(length);
  f = Eigen::VectorXd::Zero(length);

  npcorr = 0;
  kmax = 0;
  accval = 0;
}

void
Correlator::add(const double w, const unsigned int k)
{

  /// If we exceed the correlator side, the value is discarded
  if (k >= numcorrelators) {
    return;
  }
  if (k > kmax) {
    kmax = k;
  }

  /// Insert new value in shift array
  shift(k, insertindex[k]) = w;

  /// Add to average value
  if (k == 0) {
    accval += w;
  }

  /// Add to accumulator and, if needed, add to next correlator
  accumulator[k] += w;
  ++naccumulator[k];
  if (naccumulator[k] == m) {
    add(accumulator[k] / m, k + 1);
    accumulator[k] = 0;
    naccumulator[k] = 0;
  }

  /// Calculate correlation function
  const unsigned int ind1 = insertindex[k];
  // TODO: change to asserts or remove once we don't get any crashes here
  // anymore
  RUNTIME_EXP_IFN(k < numcorrelators,
                  "Cannot evaluate correlator outside its bounds.");
  RUNTIME_EXP_IFN(ind1 < p, "Cannot evaluate correlator outside its bounds.");
  if (k == 0) { /// First correlator is different
    int ind2 = static_cast<int>(ind1);
    for (unsigned int j = 0; j < p; ++j) {
      if (shift(k, ind2) > -1e10) {
        correlation(k, j) += shift(k, ind1) * shift(k, ind2);
        ncorrelation(k, j) += 1;
      }
      --ind2;
      if (ind2 < 0) {
        ind2 += p;
      }
    }
  } else {
    int ind2 = static_cast<int>(ind1) - static_cast<int>(dmin);
    for (unsigned int j = dmin; j < p; ++j) {
      if (ind2 < 0) {
        ind2 += p;
      }
      if (shift(k, ind2) > -1e10) {
        correlation(k, j) += shift(k, ind1) * shift(k, ind2);
        ncorrelation(k, j) += 1;
      }
      --ind2;
    }
  }

  ++insertindex[k];
  if (insertindex[k] == p) {
    insertindex[k] = 0;
  }
}

void
Correlator::evaluate(const bool norm)
{
  unsigned int im = 0;

  double aux = 0.;
  if (norm) {
    aux = (accval / ncorrelation(0, 0)) * (accval / ncorrelation(0, 0));
  }

  // First correlator
  for (unsigned int i = 0; i < p; ++i) {
    if (ncorrelation(0, i) > 0) {
      t[im] = i;
      f[im] = correlation(0, i) / ncorrelation(0, i) - aux;
      ++im;
    }
  }

  // Subsequent correlators
  for (int k = 1; k < kmax; ++k) {
    for (int i = dmin; i < p; ++i) {
      if (ncorrelation(k, i) > 0) {
        t[im] = i * pow(static_cast<double>(m), k);
        f[im] = correlation(k, i) / ncorrelation(k, i) - aux;
        ++im;
      }
    }
  }

  npcorr = im;
}

}
