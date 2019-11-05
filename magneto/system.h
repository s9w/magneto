#pragma once

#include <filesystem>
#include <vector>
#include <random>
#include <variant>
#include <optional>

using IndexPairVector = std::vector<std::pair<int, int>>;

namespace magneto {
	class LatticeIndexRng;

	template<class T>
	using LatticeTType = std::vector<std::vector<T>>;
	using LatticeType = LatticeTType<int>;
	using LatticeTemps = LatticeTType<double>;

	//using LatticeTType = std::vector<std::vector<double>>;

	class IsingSystem {
	public:
		IsingSystem(const int j, const double T, const int L);
		IsingSystem(const int j, const double T, const std::filesystem::path& input_path);
		IsingSystem(const int j, const std::filesystem::path& lattice_png_path, const std::filesystem::path& temp_png_path);
		void metropolis_sweeps(const IndexPairVector& lattice_indices, const std::vector<double>& random_buffer);
		void wang_sweeps(const int n = 1);
		[[nodiscard]] const LatticeType& get_lattice() const;
		size_t get_L() const;
		std::optional<double> get_temp() const;

	private:
		void metropolis_sweeps_uniform_t(const IndexPairVector& lattice_indices, const std::vector<double>& random_buffer);
		void metropolis_sweeps_variable_t(const IndexPairVector& lattice_indices, const std::vector<double>& random_buffer);

		LatticeType m_lattice;
		int m_J = 1;
		std::variant<double, LatticeTemps> m_T;
		std::vector<double> m_cached_exp_values;
	};

   struct PropertySnapshot {
      double energy;
      double energy_squared;
   };

   PropertySnapshot get_properties(const IsingSystem& system);


	LatticeType get_randomized_system(const int L);
	LatticeType get_empty_system(const int L);
	int get_dE(const LatticeType& grid, int i, int j);

   /// <summary>Returns normalized energy (per size)</summary>
   double get_E(const LatticeType& grid);

	/// <summary>Returns normalized energy squared</summary>
   double get_E_squared(const LatticeType& grid);

	/// <summary>calculates all possible values of the exp-function
	/// <para>The exponential function exp(-beta*(H2-H1)) is used extensively during the core loop 
	/// of the metropolis algorithm. Luckily, because of the nature of the Ising model, there is 
	/// only a finite amount of possible values for the parameter. In a 2D Ising model, there can 
	/// only be the values [-8J,8J] for dH. Since the values will be stored in a vector, 
	/// the parameter will be shifted to start at 0.</para>
	/// </summary>
	std::vector<double> get_cached_exp_values(const int J, const double T);
}
