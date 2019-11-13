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
      const int Ly = lattice.size();
      const int Lx = lattice[0].size();
      return { Lx, Ly };
   }
}
