#pragma once

#include <filesystem>
#include <random>
#include <variant>
#include <optional>

#include "types.h"

namespace magneto {
	class IsingSystem {
	public:
		IsingSystem(const int j, const double T, const int L);
		IsingSystem(const int j, const double T, const std::filesystem::path& input_path);
		IsingSystem(const int j, const std::filesystem::path& lattice_png_path, const std::filesystem::path& temp_png_path);
		[[nodiscard]] const LatticeType& get_lattice() const;
		[[nodiscard]] LatticeType& get_lattice_nc();
		size_t get_L() const;
		std::optional<double> get_temp() const;

	private:
		LatticeType m_lattice;
		int m_J = 1;
		std::variant<double, LatticeTemps> m_T;	};

   struct PhysicalProperties {
      double energy = 0.0;
      double magnetization = 0.0;
   };

   PhysicalProperties operator+(const PhysicalProperties& a, const PhysicalProperties& b);
   PhysicalProperties operator/(const PhysicalProperties& a, const unsigned int d);

   PhysicalProperties get_properties(const IsingSystem& system);


	LatticeType get_randomized_system(const int L);
	LatticeType get_empty_system(const int L);
	int get_dE(const LatticeType& grid, int i, int j);

   /// <summary>Returns normalized energy (per size)</summary>
   double get_E(const LatticeType& grid);

   /// <summary>Returns normalized absolute magnetization</summary>
   double get_m_abs(const LatticeType& grid);

	
	
}
