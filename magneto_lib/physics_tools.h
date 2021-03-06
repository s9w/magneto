#pragma once

#include "IsingSystem.h"

namespace magneto {
   struct PhysicsResult {
      double temp;
      double energy;
      double cv;
      double magnetization;
      double chi;
   };

	PhysicsResult get_physical_results(const PhysicalProperties& properties);
}

