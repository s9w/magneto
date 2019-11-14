#include "IsingSystem.h"
#include "logging.h"
#include "file_tools.h"

namespace {

} // namespace {}


magneto::PhysicalMeasurement magneto::operator+(const PhysicalMeasurement& a, const PhysicalMeasurement& b){
   PhysicalMeasurement sum_result(a);
   sum_result.energy += b.energy;
   sum_result.magnetization += b.magnetization;
   return sum_result;
}


magneto::PhysicalMeasurement magneto::operator/(const PhysicalMeasurement& a, const unsigned int d){
   PhysicalMeasurement div_result(a);
   div_result.energy /= d;
   div_result.magnetization /= d;
   return div_result;
}

magneto::PhysicalMeasurement magneto::get_properties(const IsingSystem& system){
   const double energy = get_E(system.get_lattice());
   const double m = get_m_abs(system.get_lattice());
   return { energy, m };
}


__declspec(noinline)
int magneto::get_dE(const LatticeType& grid, int i, int j){
   const auto [Lx, Ly] = get_dimensions_of_lattice(grid);
	return 2 * grid[i][j] * (
		grid[i][(j + 1) % Lx] +
		grid[(i + 1) % Ly][j] +
		grid[i][(j - 1 + Lx) % Lx] +
		grid[(i - 1 + Ly) % Ly][j]);
}


double magneto::get_E(const LatticeType& grid){
   const auto [Lx, Ly] = get_dimensions_of_lattice(grid);
   int E = 0;
   for (unsigned int i = 0; i < Ly; ++i) {
      for (unsigned int j = 0; j < Lx; ++j)
         E += -grid[i][j] * (grid[i][(j + 1) % Lx] + grid[(i + 1) % Ly][j]);
   }
   return E * 1.0 / (Lx * Ly);
}


double magneto::get_m_abs(const LatticeType& grid){
   const auto [Lx, Ly] = get_dimensions_of_lattice(grid);
   int m = 0;
   for (unsigned int i = 0; i < Ly; ++i) {
      for (unsigned int j = 0; j < Lx; ++j)
         m += grid[i][j];
   }
   return std::abs(m) * 1.0 / (Lx * Ly);
}


const magneto::LatticeType& magneto::IsingSystem::get_lattice() const{
	return m_lattice;
}


magneto::LatticeType& magneto::IsingSystem::get_lattice_nc() {
   return m_lattice;
}


size_t magneto::IsingSystem::get_L() const{
	return m_lattice.size();
}

int magneto::IsingSystem::get_J() const{
   return m_J;
}


magneto::IsingSystem::IsingSystem(const int j, const magneto::LatticeType& initial_state)
	: m_J(j)
	, m_lattice(initial_state)
{}


magneto::LatticeType magneto::get_randomized_system(const int Lx, const int Ly) {
   unsigned seed1 = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
   std::mt19937_64 generator(seed1);
   std::uniform_int_distribution<int> dist(0, 1);
   magneto::LatticeType grid(Ly, std::vector<char>(Lx));
   for (int i = 0; i < Ly; ++i) {
      for (int j = 0; j < Lx; ++j) {
         grid[i][j] = static_cast<char>(dist(generator) * 2 - 1);
      }
   }
   return grid;
}
