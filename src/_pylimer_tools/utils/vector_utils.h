#ifndef VECTOR_UTILS_H
#define VECTOR_UTILS_H

#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <igraph/igraph.h>

namespace pylimer_tools
{
  namespace utils
  {
    template <typename IN1,
              typename IN2,
              typename OUT>
    inline OUT interleave(IN1 it1,
                          IN1 end1,
                          IN2 it2,
                          IN2 end2,
                          OUT out)
    {
      // interleave until at least one container is done
      while (it1 != end1 && it2 != end2)
      {
        // insert from container 1
        *out = *it1;
        out++;
        it1++;
        // insert from container 2
        *out = *it2;
        out++;
        it2++;
      }
      if (it1 != end1) // check and finish container 1
      {
        return std::copy(it1, end1, out);
      }
      else if (it2 != end2) // check and finish container 2
      {
        return std::copy(it2, end2, out);
      }
      return out; // both done
    }

    template <typename IN1>
    inline void StdVectorToIgraphVectorT(IN1 &vectR, igraph_vector_t *v)
    {
      size_t n = vectR.size();

      /* Make sure that there is enough space for the items in v */
      igraph_vector_resize(v, n);

      /* Copy all the items */
      for (size_t i = 0; i < n; i++)
      {
        VECTOR(*v)
        [i] = vectR[i];
      }
    }
  }
}

#endif
