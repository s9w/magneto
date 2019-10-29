#include <chrono>
#include <random>

#include "system.h"

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


void magneto::metropolis_sweeps(
	LatticeType& grid,
	const IndexPairVector& lattice_indices,
	const std::vector<double>& exp_values,
	const std::vector<double>& random_buffer,
	const PhysicsSettings& physics
){
	const int L = static_cast<int>(grid.size());
	std::uniform_int_distribution <int> dist_grid(0, L - 1);
	std::uniform_real_distribution <double > dist_one(0.0, 1.0);
	const int buffer_offset = get_exp_buffer_offset(physics.J);
	int flip_i, flip_j;
	int dE;
	for(int i=0; i < random_buffer.size(); ++i){
		flip_i = lattice_indices[i].first;
		flip_j = lattice_indices[i].second;
		dE = physics.J * get_dE(grid, flip_i, flip_j);
		if (dE <= 0 || (random_buffer[i] < exp_values[dE + buffer_offset]))
			grid[flip_i][flip_j] *= -1;
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
