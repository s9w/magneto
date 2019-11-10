#include "physics_tools.h"
#include <numeric>

namespace {
   ///<summary>Computes the mean of anything that implements + and / operators</summary>
   template <class T>
   T get_mean(const std::vector<T>& values) {
      const T sum = std::accumulate(values.cbegin(), values.cend(), T());
      return sum / static_cast<int>(values.size());
   }


   ///<summary>Mean of squared elements</summary>
   template <class T>
   T get_squared_mean(const std::vector<T>& values) {
      const auto& square_op = [](const T init, const T elem) {
         return init + elem * elem;
      };
      const T sum = std::accumulate(values.cbegin(), values.cend(), T(), square_op);
      return sum / static_cast<int>(values.size());
   }


   ///<summary>Calculates variance with the old trick of Var(x) = |x^2|-|x|^2</summary>
   double get_variance(const std::vector<double>& values) {
      const double mean = get_mean(values);
      return get_squared_mean(values) - mean * mean;
   }


   std::vector<double> get_energies(const std::vector<magneto::PhysicalProperties>& properties) {
      std::vector<double> energies;
      energies.reserve(properties.size());
      for (const magneto::PhysicalProperties& property : properties)
         energies.emplace_back(property.energy);
      return energies;
   }


   std::vector<double> get_mags(const std::vector<magneto::PhysicalProperties>& properties) {
      std::vector<double> energies;
      energies.reserve(properties.size());
      for (const magneto::PhysicalProperties& property : properties)
         energies.emplace_back(property.magnetization);
      return energies;
   }


   double get_energy_variance(const std::vector<magneto::PhysicalProperties>& properties) {
      return get_variance(get_energies(properties));
   }

   double get_mag_variance(const std::vector<magneto::PhysicalProperties>& properties) {
      return get_variance(get_mags(properties));
   }

} // namespace {}


magneto::PhysicsResult magneto::get_physical_results(
	const std::vector<magneto::PhysicalProperties>& properties, const unsigned int L, const double T
){
   const double mean_energy = get_mean(get_energies(properties));
   const double mean_magnetization = get_mean(get_mags(properties));
	const double cv = get_energy_variance(properties) * L * L / (T * T);
	const double chi = get_mag_variance(properties) * L * L / T;
	return { T, mean_energy, cv, mean_magnetization, chi };
}
