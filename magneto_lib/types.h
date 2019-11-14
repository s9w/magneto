#pragma once

#include <string>
#include <vector>

namespace magneto {
   using IndexPairVector = std::vector<std::pair<int, int>>;

   template<class T>
   using LatticeTType = std::vector<std::vector<T>>;

   using LatticeType = LatticeTType<char>;
   using LatticeIType = LatticeTType<int>;
   using LatticeDType = LatticeTType<double>;

   template<class T>
   std::pair<unsigned int, unsigned int> get_dimensions_of_lattice(const magneto::LatticeTType<T>& lattice) {
      const unsigned int Ly = static_cast<unsigned int>(lattice.size());
      const unsigned int Lx = static_cast<unsigned int>(lattice[0].size());
      return { Lx, Ly };
   }
}
