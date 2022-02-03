#ifndef VECTOR_UTILS_H
#define VECTOR_UTILS_H

#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
extern "C"
{
#include <igraph/igraph.h>
}
#include <cassert>

namespace pylimer_tools
{
  namespace utils
  {
    template <typename IN>
    static inline std::vector<IN> interleave(std::vector<IN> in1,
                                             std::vector<IN> in2)
    {
      size_t size = in1.size();
      assert(size == in2.size());
      std::vector<IN> out;
      out.reserve(2 * size);
      // interleave until at least one container is done
      for (size_t i = 0; i < size; ++i)
      {
        out.push_back(in1[i]);
        out.push_back(in2[i]);
      }

      return out; // both done
    }

    template <typename IN>
    static inline void eraseIndices(std::vector<IN> from, std::vector<long int> indices)
    {
      for (auto index : indices)
      {
        from.erase(index);
      }
    }

    template <typename IN1>
    static inline void StdVectorToIgraphVectorT(IN1 &vectR, igraph_vector_t *v)
    {
      size_t n = vectR.size();

      /* Make sure that there is enough space for the items in v */
      igraph_vector_resize(v, n);

      /* Copy all the items */
      for (size_t i = 0; i < n; ++i)
      {
        igraph_vector_set(v, i, vectR[i]);
      }
    }

    template <typename IN>
    static inline void igraphVectorTToStdVector(igraph_vector_t *v, std::vector<IN> &vectL)
    {
      long n = igraph_vector_size(v);

      vectL.clear();
      vectL.reserve(n);

      for (long i = 0; i < n; ++i)
      {
        vectL.push_back(igraph_vector_e(v, i));
      }
    }

    template <typename IN>
    static inline std::vector<IN> initializeWithValue(size_t n, IN value)
    {
      std::vector<IN> result;
      result.reserve(n);
      for (size_t i = 0; i < n; ++i)
      {
        result.push_back(value);
      }
      return result;
    }
  }
}

#endif
