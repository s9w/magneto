#pragma once

#include <string>
#include <vector>

namespace magneto {
   using IndexPairVector = std::vector<std::pair<int, int>>;

   template<class T>
   using LatticeTType = std::vector<std::vector<T>>;

   using LatticeType = LatticeTType<char>;
   using LatticeIType = LatticeTType<int>;
   using LatticeTemps = LatticeTType<double>;
}
