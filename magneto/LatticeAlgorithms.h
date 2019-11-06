#pragma once

#include "types.h"
#include <future>


namespace magneto {   
   class Metropolis {
   public:
      Metropolis(const int J, const double T, const int L, const int max_rng_threads);
      ~Metropolis();
      void run(LatticeType& lattice);

   private:
      IndexPairVector m_lattice_index_buffer;
      std::vector<double> m_random_buffer;
      std::vector<double> m_cached_exp_values;
      std::vector<std::future<IndexPairVectorReturn>> m_lattice_index_futures;
      std::vector<std::future<UniformRandomBufferReturn>> m_future_random_buffers;
      const int m_J;
   };

   /// <summary>calculates all possible values of the exp-function
   /// <para>The exponential function exp(-beta*(H2-H1)) is used extensively during the core loop 
   /// of the metropolis algorithm. Luckily, because of the nature of the Ising model, there is 
   /// only a finite amount of possible values for the parameter. In a 2D Ising model, there can 
   /// only be the values [-8J,8J] for dH. Since the values will be stored in a vector, 
   /// the parameter will be shifted to start at 0.</para>
   /// </summary>
   std::vector<double> get_cached_exp_values(const int J, const double T);

}
