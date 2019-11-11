#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>

#include "IsingSystem.h"
#include "logging.h"

namespace {

   /// <summary>Converts [0,255] to [-1,1]. </summary>
   template<class TIn, class TOut>
	TOut get_pm1_from_255_value(const TIn value) {
		return static_cast<TOut>(value / 255 * 2 - 1);
	}


   // Making sure this will not be used in an unsigned char
   template<>
   unsigned char get_pm1_from_255_value(const unsigned char value) = delete;


	magneto::LatticeType get_lattice_from_monochrome_bitmap_data(unsigned char* png_data, const int L) {
		magneto::LatticeType lattice(L, std::vector<char>(L));
		for (int i = 0; i < L; ++i) {
			for (int j = 0; j < L; ++j) {
				lattice[i][j] = get_pm1_from_255_value<unsigned char, char>(png_data[i * L + j]);
			}
		}
		return lattice;
	}


	magneto::LatticeType get_lattice_from_rgba_bitmap_data(unsigned char* png_data, const int L) {
		magneto::LatticeType lattice(L, std::vector<char>(L));
		for (int i = 0; i < L; ++i) {
			for (int j = 0; j < L; ++j) {
				unsigned int value = 0;
				for (int k = 0; k < 4; ++k) {
					value += png_data[(i * L + j) * 4 + k];
				}
				value /= 4;
				lattice[i][j] = get_pm1_from_255_value<unsigned int, char>(value);
			}
		}
		return lattice;
	}


	magneto::LatticeType get_lattice_from_png_data(unsigned char* png_data, const int bpp, const int L) {
		if (bpp == 1)
			return get_lattice_from_monochrome_bitmap_data(png_data, L);
		else if (bpp == 4)
			return get_lattice_from_rgba_bitmap_data(png_data, L);
	}


	magneto::LatticeType get_lattice_from_png_file(const std::filesystem::path& path) {
		int x, y, bpp;
		unsigned char* data = stbi_load(path.string().c_str(), &x, &y, &bpp, 0);
		magneto::LatticeType lattice = get_lattice_from_png_data(data, bpp, x);
		stbi_image_free(data);
		return lattice;
	}


	magneto::LatticeTemps get_lattice_temps_from_png_file(
		const std::filesystem::path& path, const double temp_min, const double temp_max
	) {
		int x, y, bpp;
		unsigned char* image_data = stbi_load(path.string().c_str(), &x, &y, &bpp, 0);
		magneto::LatticeTemps temps(x, std::vector<double>(x, 2.26));
      if (x != y) {
         magneto::get_logger()->error("Temperature input image is not square.");
         return temps;
      }
		const double temp_factor = temp_max - temp_min;
		for (int i = 0; i < x; ++i) {
			for (int j = 0; j < x; ++j) {
				int value = 0;
				for (int k = 0; k < bpp; ++k) {
					value += image_data[(i * x + j) * bpp + k];
				}
				value /= bpp;
				temps[i][j] = temp_min + value * 1.0 / 256 * temp_factor;
			}
		}
		stbi_image_free(image_data);
		return temps;
	}


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
}


magneto::PhysicalProperties magneto::operator+(const PhysicalProperties& a, const PhysicalProperties& b){
   PhysicalProperties sum_result(a);
   sum_result.energy += b.energy;
   sum_result.magnetization += b.magnetization;
   return sum_result;
}


magneto::PhysicalProperties magneto::operator/(const PhysicalProperties& a, const unsigned int d){
   PhysicalProperties div_result(a);
   div_result.energy /= d;
   div_result.magnetization /= d;
   return div_result;
}

magneto::PhysicalProperties magneto::get_properties(const IsingSystem& system){
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


magneto::IsingSystem::IsingSystem(const int j, const double T, const std::filesystem::path& input_path)
	: m_J(j)
	, m_T(T)
	, m_lattice(get_lattice_from_png_file(input_path))
{}


magneto::IsingSystem::IsingSystem(
	const int j, 
	const std::filesystem::path& lattice_png_path,
	const std::filesystem::path& temp_png_path
)
	: m_J(j)
	, m_T(get_lattice_temps_from_png_file(temp_png_path, 1.0, 2.3))
	, m_lattice(get_lattice_from_png_file(lattice_png_path))
{
}
