#pragma once

#include <filesystem>
#include <vector>
#include <random>

using IndexPairVector = std::vector<std::pair<int, int>>;

namespace magneto {
	class LatticeIndexRng;

	using LatticeType = std::vector<std::vector<int>>;

	class IsingSystem {
	public:
		IsingSystem(const int j, const double T, const int L);
		IsingSystem(const int j, const double T, const int L, unsigned char* png_data, const int bpp);
		void metropolis_sweeps(const IndexPairVector& lattice_indices, const std::vector<double>& rng_buffer);
		[[nodiscard]] const LatticeType& get_lattice() const;
		size_t get_L() const;

	private:
		LatticeType m_lattice;
		int m_J = 1;
		double m_T = 2.2;
		std::vector<double> m_cached_exp_values;
	};

	LatticeType get_randomized_system(const int L);
	LatticeType get_empty_system(const int L);
	int get_dE(const LatticeType& grid, int i, int j);

	/// <summary>calculates all possible values of the exp-function
	/// <para>The exponential function exp(-beta*(H2-H1)) is used extensively during the core loop 
	/// of the metropolis algorithm. Luckily, because of the nature of the Ising model, there is 
	/// only a finite amount of possible values for the parameter. In a 2D Ising model, there can 
	/// only be the values [-8J,8J] for dH. Since the values will be stored in a vector, 
	/// the parameter will be shifted to start at 0.</para>
	/// </summary>
	std::vector<double> get_cached_exp_values(const int J, const double T);
}
