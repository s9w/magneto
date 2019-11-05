#include <chrono>
#include <random>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "system.h"

namespace {
	constexpr int get_exp_buffer_offset(int J) {
		// std::abs is not constexpr :(
		return J>0? 8 * J : -8 * J;
	}


	int get_pm1_from_255_value(const int value) {
		return value / 255 * 2 - 1;
	}


	magneto::LatticeType get_lattice_from_monochrome_bitmap_data(unsigned char* png_data, const int L) {
		magneto::LatticeType lattice(L, std::vector<int>(L));
		for (int i = 0; i < L; ++i) {
			for (int j = 0; j < L; ++j) {
				lattice[i][j] = get_pm1_from_255_value(png_data[i * L + j]);
			}
		}
		return lattice;
	}


	magneto::LatticeType get_lattice_from_rgba_bitmap_data(unsigned char* png_data, const int L) {
		magneto::LatticeType lattice(L, std::vector<int>(L));
		for (int i = 0; i < L; ++i) {
			for (int j = 0; j < L; ++j) {
				int value = 0;
				for (int k = 0; k < 4; ++k) {
					value += png_data[(i * L + j) * 4 + k];
				}
				value /= 4;
				lattice[i][j] = get_pm1_from_255_value(value);
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
		unsigned char* data = stbi_load(path.string().c_str(), &x, &y, &bpp, 0);
		magneto::LatticeTemps temps(x, std::vector<double>(x));
		const double temp_factor = temp_max - temp_min;
		for (int i = 0; i < x; ++i) {
			for (int j = 0; j < x; ++j) {
				int value = 0;
				for (int k = 0; k < 4; ++k) {
					value += data[(i * x + j) * 4 + k];
				}
				value /= 4;
				temps[i][j] = temp_min + value * 1.0 / 256 * temp_factor;
			}
		}
		stbi_image_free(data);
		return temps;
	}
}


magneto::LatticeType magneto::get_randomized_system(const int L){
	unsigned seed1 = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
	std::default_random_engine generator(seed1);
	std::uniform_int_distribution <int> dist(0, 1);
	LatticeType grid(L, std::vector<int>(L));
	for (int i = 0; i < L; ++i) {
		for (int j = 0; j < L; ++j) {
			grid[i][j] = dist(generator) * 2 - 1;
		}
	}
	return grid;
}


magneto::LatticeType magneto::get_empty_system(const int L){
	return LatticeType(L, std::vector<int>(L, -1));
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

void magneto::IsingSystem::metropolis_sweeps(
	const IndexPairVector& lattice_indices,
	const std::vector<double>& random_buffer
){
	if (std::holds_alternative<double>(m_T))
		metropolis_sweeps_uniform_t(lattice_indices, random_buffer);
	else
		metropolis_sweeps_variable_t(lattice_indices, random_buffer);
}


const magneto::LatticeType& magneto::IsingSystem::get_lattice() const{
	return m_lattice;
}


size_t magneto::IsingSystem::get_L() const{
	return m_lattice.size();
}

__declspec(noinline)
void magneto::IsingSystem::metropolis_sweeps_uniform_t(
	const IndexPairVector& lattice_indices, 
	const std::vector<double>& random_buffer
){
	const int L = static_cast<int>(m_lattice.size());
	std::uniform_int_distribution <int> dist_grid(0, L - 1);
	std::uniform_real_distribution <double > dist_one(0.0, 1.0);
	const int buffer_offset = get_exp_buffer_offset(m_J);
	int flip_i, flip_j;
	int dE;
	for (int i = 0; i < random_buffer.size(); ++i) {
		flip_i = lattice_indices[i].first;
		flip_j = lattice_indices[i].second;
		dE = m_J * get_dE(m_lattice, flip_i, flip_j);
		if (dE <= 0 || (random_buffer[i] < m_cached_exp_values[dE + buffer_offset]))
			m_lattice[flip_i][flip_j] *= -1;
	}
}

__declspec(noinline)
void magneto::IsingSystem::metropolis_sweeps_variable_t(
	const IndexPairVector& lattice_indices,
	const std::vector<double>& random_buffer
){
	const int L = static_cast<int>(m_lattice.size());
	std::uniform_int_distribution <int> dist_grid(0, L - 1);
	std::uniform_real_distribution <double > dist_one(0.0, 1.0);
	const int buffer_offset = get_exp_buffer_offset(m_J);
	int flip_i, flip_j;
	int dE;
	for (int i = 0; i < random_buffer.size(); ++i) {
		flip_i = lattice_indices[i].first;
		flip_j = lattice_indices[i].second;
		dE = m_J * get_dE(m_lattice, flip_i, flip_j);
		const double exp_value = exp(-dE / std::get<LatticeTemps>(m_T)[flip_i][flip_j]);
		if (dE <= 0 || (random_buffer[i] < exp_value))
			m_lattice[flip_i][flip_j] *= -1;
	}
}


std::vector<double> magneto::get_cached_exp_values(const int J, const double T){
	std::vector<double> exp_values;
	const int min_value = -8 * J;
	const int max_value = 8 * J;
	const int value_count = max_value - min_value + 1;
	const int buffer_offset = get_exp_buffer_offset(J);
	for (int dE = 0; dE < value_count; ++dE)
		exp_values.push_back(exp(-(dE - buffer_offset) / T));
	return exp_values;
}


magneto::IsingSystem::IsingSystem(const int j, const double T, const int L)
	: m_J(j)
	, m_T(T)
	, m_lattice(get_randomized_system(L))
	, m_cached_exp_values(get_cached_exp_values(m_J, std::get<double>(m_T)))
{}


magneto::IsingSystem::IsingSystem(const int j, const double T, const std::filesystem::path& input_path)
	: m_J(j)
	, m_T(T)
	, m_lattice(get_lattice_from_png_file(input_path))
	, m_cached_exp_values(get_cached_exp_values(m_J, std::get<double>(m_T)))
{}


magneto::IsingSystem::IsingSystem(
	const int j, 
	const std::filesystem::path& lattice_png_path,
	const std::filesystem::path& temp_png_path
)
	: m_J(j)
	, m_T(get_lattice_temps_from_png_file(temp_png_path, 0.0, 3.0))
	, m_lattice(get_lattice_from_png_file(lattice_png_path))
{
}
