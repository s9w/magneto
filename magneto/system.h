#pragma once

#include <filesystem>
#include <vector>
#include <random>

namespace magneto {
	class LatticeIndexRng;

	using GridType = std::vector<std::vector<int>>;

	struct PhysicsSettings {
		int J = 1;
		double T = 2.2;
	};

	GridType get_randomized_system(const int L);
	int get_dE(const GridType& grid, int i, int j);
	void metropolis_sweeps(GridType& grid, LatticeIndexRng& lattice_rng, const std::vector<double>& exp_values, const std::vector<double>& rng_buffer, const PhysicsSettings& physics);

	/// <summary>calculates all possible values of the exp-function
	/// <para>The exponential function exp(-beta*(H2-H1)) is used extensively during the core loop 
	/// of the metropolis algorithm. Luckily, because of the nature of the Ising model, there is 
	/// only a finite amount of possible values for the parameter. In a 2D Ising model, there can 
	/// only be the values [-8J,8J] for dH. Since the values will be stored in a vector, 
	/// the parameter will be shifted to start at 0.</para>
	/// </summary>
	std::vector<double> get_cached_exp_values(const PhysicsSettings& physics);
}
