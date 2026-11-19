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
#include <climits>
#include <cmath>
#include <cstdint>
#include <limits>

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
  INVALIDARG_EXP_IFN(numcorrin > 0 && pin > 0 && min > 0,
                     "numcorrelators, p, and m must all be positive.");
  INVALIDARG_EXP_IFN(numcorrin < INT_MAX && pin < INT_MAX && min < INT_MAX,
                     "Arguments must be able to cast to integers.");
  INVALIDARG_EXP_IFN(static_cast<uint64_t>(numcorrin) *
                         static_cast<uint64_t>(pin) <=
                       static_cast<uint64_t>(INT_MAX),
                     "numcorrelators * p overflows integer limit.");
  INVALIDARG_EXP_IFN(
    numcorrin == 1 || (min >= 2 && pin >= min),
    "For multi-level correlators (numcorrelators > 1), require m >= 2 and "
    "p >= m.");
  numcorrelators = numcorrin;
  p = pin;
  m = min;
  dmin = p / m;

  length = numcorrelators * p;

  const Eigen::Index nCorrIdx = static_cast<Eigen::Index>(numcorrelators);
  const Eigen::Index pIdx = static_cast<Eigen::Index>(p);
  const Eigen::Index lenIdx = static_cast<Eigen::Index>(length);

  // Initialize shift slots to NaN to distinguish empty slots from any valid
  // finite signal values (including large negative numbers)
  shift = Eigen::MatrixXd::Constant(
    nCorrIdx, pIdx, std::numeric_limits<double>::quiet_NaN());
  correlation = Eigen::MatrixXd::Zero(nCorrIdx, pIdx);
  ncorrelation = MatrixXuli::Zero(nCorrIdx, pIdx);
  accumulator = Eigen::VectorXd::Zero(nCorrIdx);
  naccumulator = VectorXui::Zero(nCorrIdx);
  insertindex = VectorXui::Zero(nCorrIdx);

  t = Eigen::VectorXd::Zero(lenIdx);
  f = Eigen::VectorXd::Zero(lenIdx);

  npcorr = 0;
  kmax = 0;
  accval = 0.0;
}

void
Correlator::reset()
{
  setsize(numcorrelators, p, m);
}

void
Correlator::add(const double w, const unsigned int k)
{
  INVALIDARG_EXP_IFN(std::isfinite(w),
                     "Correlator input value must be finite.");

  /// If we exceed the correlator size, the value is discarded
  if (k >= numcorrelators) {
    return;
  }
  if (k > kmax) {
    kmax = k;
  }

  const Eigen::Index kIdx = static_cast<Eigen::Index>(k);

  /// Insert new value in shift array
  shift(kIdx, static_cast<Eigen::Index>(insertindex[kIdx])) = w;

  /// Add to average value
  if (k == 0) {
    accval += w;
  }

  /// Add to accumulator and, if needed, add to next correlator
  accumulator[kIdx] += w;
  ++naccumulator[kIdx];
  if (naccumulator[kIdx] == m) {
    add(accumulator[kIdx] / static_cast<double>(m), k + 1);
    accumulator[kIdx] = 0.0;
    naccumulator[kIdx] = 0;
  }

  /// Calculate correlation function
  const unsigned int ind1 = insertindex[kIdx];
  const Eigen::Index ind1Idx = static_cast<Eigen::Index>(ind1);
  RUNTIME_EXP_IFN(k < numcorrelators,
                  "Cannot evaluate correlator outside its bounds.");
  RUNTIME_EXP_IFN(ind1 < p, "Cannot evaluate correlator outside its bounds.");
  if (k == 0) { /// First correlator is different
    int ind2 = static_cast<int>(ind1);
    for (unsigned int j = 0; j < p; ++j) {
      const Eigen::Index ind2Idx = static_cast<Eigen::Index>(ind2);
      if (!std::isnan(shift(kIdx, ind2Idx))) {
        const Eigen::Index jIdx = static_cast<Eigen::Index>(j);
        correlation(kIdx, jIdx) += shift(kIdx, ind1Idx) * shift(kIdx, ind2Idx);
        ncorrelation(kIdx, jIdx) += 1;
      }
      --ind2;
      if (ind2 < 0) {
        ind2 += static_cast<int>(p);
      }
    }
  } else {
    int ind2 = static_cast<int>(ind1) - static_cast<int>(dmin);
    for (unsigned int j = dmin; j < p; ++j) {
      if (ind2 < 0) {
        ind2 += static_cast<int>(p);
      }
      const Eigen::Index ind2Idx = static_cast<Eigen::Index>(ind2);
      if (!std::isnan(shift(kIdx, ind2Idx))) {
        const Eigen::Index jIdx = static_cast<Eigen::Index>(j);
        correlation(kIdx, jIdx) += shift(kIdx, ind1Idx) * shift(kIdx, ind2Idx);
        ncorrelation(kIdx, jIdx) += 1;
      }
      --ind2;
    }
  }

  ++insertindex[kIdx];
  if (insertindex[kIdx] == p) {
    insertindex[kIdx] = 0;
  }
}

void
Correlator::addBatch(const Eigen::Ref<const Eigen::VectorXd>& values)
{
  INVALIDARG_EXP_IFN(values.allFinite(),
                     "Correlator batch values must all be finite.");
  for (Eigen::Index i = 0; i < values.size(); ++i) {
    add(values(i));
  }
}

void
Correlator::evaluate(const bool norm)
{
  Eigen::Index im = 0;

  double aux = 0.0;
  if (norm && ncorrelation(0, 0) > 0) {
    const double n0 = static_cast<double>(ncorrelation(0, 0));
    aux = (accval / n0) * (accval / n0);
  }

  // First correlator
  for (unsigned int i = 0; i < p; ++i) {
    const Eigen::Index iIdx = static_cast<Eigen::Index>(i);
    if (ncorrelation(0, iIdx) > 0) {
      t[im] = static_cast<double>(i);
      f[im] =
        correlation(0, iIdx) / static_cast<double>(ncorrelation(0, iIdx)) - aux;
      ++im;
    }
  }

  // Subsequent correlators (including k == kmax)
  for (unsigned int k = 1; k <= kmax && k < numcorrelators; ++k) {
    const Eigen::Index kIdx = static_cast<Eigen::Index>(k);
    const double scale = std::pow(static_cast<double>(m), static_cast<int>(k));
    for (unsigned int i = dmin; i < p; ++i) {
      const Eigen::Index iIdx = static_cast<Eigen::Index>(i);
      if (ncorrelation(kIdx, iIdx) > 0) {
        t[im] = static_cast<double>(i) * scale;
        f[im] = correlation(kIdx, iIdx) /
                  static_cast<double>(ncorrelation(kIdx, iIdx)) -
                aux;
        ++im;
      }
    }
  }

  npcorr = static_cast<unsigned int>(im);
}

Eigen::VectorXd
Correlator::getTimes() const
{
  return t.head(static_cast<Eigen::Index>(npcorr));
}

Eigen::VectorXd
Correlator::getCorrelation() const
{
  return f.head(static_cast<Eigen::Index>(npcorr));
}

unsigned long int
Correlator::getNumSamples() const
{
  return (ncorrelation.rows() > 0 && ncorrelation.cols() > 0)
           ? ncorrelation(0, 0)
           : 0;
}

double
Correlator::getMean() const
{
  const unsigned long int n = getNumSamples();
  return n > 0 ? accval / static_cast<double>(n) : 0.0;
}

} // namespace pylimer_tools::calc
