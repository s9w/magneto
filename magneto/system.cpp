#include <chrono>
#include <random>

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

int magneto::get_dE(const LatticeType& grid, int i, int j){
	const int L = static_cast<int>(grid.size());
	return 2 * grid[i][j] * (
		grid[i][(j + 1) % L] +
		grid[(i + 1) % L][j] +
		grid[i][(j - 1 + L) % L] +
		grid[(i - 1 + L) % L][j]);
}


void magneto::IsingSystem::metropolis_sweeps(
	const IndexPairVector& lattice_indices,
	const std::vector<double>& random_buffer
){
	const int L = static_cast<int>(m_lattice.size());
	std::uniform_int_distribution <int> dist_grid(0, L - 1);
	std::uniform_real_distribution <double > dist_one(0.0, 1.0);
	const int buffer_offset = get_exp_buffer_offset(m_J);
	int flip_i, flip_j;
	int dE;
	for(int i=0; i < random_buffer.size(); ++i){
		flip_i = lattice_indices[i].first;
		flip_j = lattice_indices[i].second;
		dE = m_J * get_dE(m_lattice, flip_i, flip_j);
		if (dE <= 0 || (random_buffer[i] < m_cached_exp_values[dE + buffer_offset]))
			m_lattice[flip_i][flip_j] *= -1;
	}
}


const magneto::LatticeType& magneto::IsingSystem::get_lattice() const{
	return m_lattice;
}


size_t magneto::IsingSystem::get_L() const{
	return m_lattice.size();
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
	: m_J(j), m_T(T)
	, m_lattice(get_randomized_system(L))
	, m_cached_exp_values(get_cached_exp_values(m_J, m_T))
{}


magneto::IsingSystem::IsingSystem(const int j, const double T, const int L, unsigned char* png_data, const int bpp)
	: m_J(j), m_T(T)
	, m_lattice(get_lattice_from_png_data(png_data, bpp, L))
	, m_cached_exp_values(get_cached_exp_values(m_J, m_T))
{}
