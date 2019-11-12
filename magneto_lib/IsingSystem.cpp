#include "IsingSystem.h"
#include "logging.h"
#include "file_tools.h"

namespace {

   magneto::LatticeType get_randomized_system(const int L) {
      unsigned seed1 = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
      std::mt19937_64 generator(seed1);
      std::uniform_int_distribution<int> dist(0, 1);
      magneto::LatticeType grid(L, std::vector<char>(L));
      for (int i = 0; i < L; ++i) {
         for (int j = 0; j < L; ++j) {
            grid[i][j] = static_cast<char>(dist(generator) * 2 - 1);
         }
      }
      return grid;
   }


   magneto::LatticeType get_empty_system(const int L) {
      return magneto::LatticeType(L, std::vector<char>(L, -1));
   }

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
	const int L = static_cast<int>(grid.size());
	return 2 * grid[i][j] * (
		grid[i][(j + 1) % L] +
		grid[(i + 1) % L][j] +
		grid[i][(j - 1 + L) % L] +
		grid[(i - 1 + L) % L][j]);
}


double magneto::get_E(const LatticeType& grid){
   const int L = static_cast<int>(grid.size());
   int E = 0;
   for (int i = 0; i < L; ++i) {
      for (int j = 0; j < L; ++j)
         E += -grid[i][j] * (grid[i][(j + 1) % L] + grid[(i + 1) % L][j]);
   }
   return E * 1.0 / (L * L);
}


double magneto::get_m_abs(const LatticeType& grid){
   const int L = static_cast<int>(grid.size());
   int m = 0;
   for (int i = 0; i < L; ++i) {
      for (int j = 0; j < L; ++j)
         m += grid[i][j];
   }
   return std::abs(m) * 1.0 / (L * L);
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

std::optional<double> magneto::IsingSystem::get_temp() const{
	if (!std::holds_alternative<double>(m_T))
		return std::nullopt;
	return std::get<double>(m_T);
}


magneto::IsingSystem::IsingSystem(const int j, const double T, const int L)
	: m_J(j)
	, m_T(T)
	, m_lattice(get_randomized_system(L))
{}


magneto::IsingSystem::IsingSystem(const int j, const LatticeDType& T, const int L)
   : m_J(j)
   , m_T(T)
   , m_lattice(get_randomized_system(L))
{}


magneto::IsingSystem::IsingSystem(const int j, const double T, const std::filesystem::path& input_path)
	: m_J(j)
	, m_T(T)
	, m_lattice(get_lattice_from_png_file(input_path))
{}

