#include <chrono>
#include <random>

#include "system.h"
#include "LatticeIndexRng.h"

namespace {
	constexpr int get_exp_buffer_offset(int J) {
		// std::abs is not constexpr :(
		return J>0? 8 * J : -8 * J;
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


int magneto::get_dE(const LatticeType& grid, int i, int j){
	const int L = static_cast<int>(grid.size());
	return 2 * grid[i][j] * (
		grid[i][(j + 1) % L] +
		grid[(i + 1) % L][j] +
		grid[i][(j - 1 + L) % L] +
		grid[(i - 1 + L) % L][j]);
}


void magneto::metropolis_sweeps(
	LatticeType& grid,
	LatticeIndexRng& lattice_rng,
	const std::vector<double>& exp_values,
	const std::vector<double>& random_buffer,
	const PhysicsSettings& physics
){
	const int L = static_cast<int>(grid.size());
	std::uniform_int_distribution <int> dist_grid(0, L - 1);
	std::uniform_real_distribution <double > dist_one(0.0, 1.0);
	const int buffer_offset = get_exp_buffer_offset(physics.J);
	int flip_index;
	int flipi, flipj;
	int dE;
	for(const double random_value : random_buffer){
		flip_index = lattice_rng.get();
		flipi = flip_index / L;
		flipj = flip_index % L;
		dE = physics.J * get_dE(grid, flipi, flipj);
		if (dE <= 0 || (random_value < exp_values[dE + buffer_offset]))
			grid[flipi][flipj] *= -1;
	}
}


std::vector<double> magneto::get_cached_exp_values(const PhysicsSettings& physics){
	std::vector<double> exp_values;
	const int min_value = -8 * physics.J;
	const int max_value = 8 * physics.J;
	const int value_count = max_value - min_value + 1;
	const int buffer_offset = get_exp_buffer_offset(physics.J);
	for (int dE = 0; dE < value_count; ++dE)
		exp_values.push_back(exp(-(dE - buffer_offset) / physics.T));
	return exp_values;
}
