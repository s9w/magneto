#pragma once

#include "system.h"

namespace magneto {
   struct PhysicsResult {
      double temp;
      double energy;
      double cv;
      double magnetization;
      double chi;
   };

	PhysicsResult get_physical_results(const std::vector<magneto::PhysicalProperties>& properties, const unsigned int L, const double T);
}

