#pragma once

#include <string>
#include <vector>

namespace magneto {
   using IndexPairVector = std::vector<std::pair<int, int>>;

   template<class T>
   using LatticeTType = std::vector<std::vector<T>>;
   using LatticeType = LatticeTType<int>;
   using LatticeTemps = LatticeTType<double>;

   template<class T>
   struct ResultWithThreadId {
      operator T() {
         return m_result;
      }
      T m_result;
      std::string m_thread_id;
   };

   using UniformRandomBufferReturn = ResultWithThreadId<std::vector<double>>;
   using IndexPairVectorReturn = ResultWithThreadId<IndexPairVector>;
}
