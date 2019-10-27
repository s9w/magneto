#include <chrono>
#include <random>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "system.h"

namespace {
	unsigned char get_png_value_from_pm_one(const int value) {
		return (value + 1) / 2 * 255;
	}

	constexpr int get_exp_buffer_offset(int J) {
		// std::abs is not constexpr :(
		return J>0? 8 * J : -8 * J;
	}
}


magneto::GridType magneto::get_randomized_system(const int L){
	unsigned seed1 = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
	std::default_random_engine generator(seed1);
	std::uniform_int_distribution <int> dist(0, 1);
	GridType grid(L, std::vector<int>(L));
	for (int i = 0; i < L; ++i) {
		for (int j = 0; j < L; ++j) {
			grid[i][j] = dist(generator) * 2 - 1;
		}
	}
	return grid;
}


void magneto::write_png(const GridType& grid, const std::filesystem::path& path){
	std::vector<unsigned char> grid_png;
	const int L = static_cast<int>(grid.size());
	grid_png.reserve(L * L);
	for (const auto& row : grid)
		for (const auto& elem : row)
			grid_png.emplace_back(get_png_value_from_pm_one(elem));

	stbi_write_png(path.string().c_str(), L, L, 1, grid_png.data(), L * 1);
}


int magneto::get_dE(const GridType& grid, int i, int j){
	const int L = static_cast<int>(grid.size());
	return 2 * grid[i][j] * (
		grid[i][(j + 1) % L] +
		grid[(i + 1) % L][j] +
		grid[i][(j - 1 + L) % L] +
		grid[(i - 1 + L) % L][j]);
}


void magneto::metropolis_sweeps(
	GridType& grid, 
	std::mt19937& rng,
	const std::vector<double>& exp_values, 
	const PhysicsSettings& physics, 
	const int n
){
	const int L = static_cast<int>(grid.size());
	std::uniform_int_distribution <int> dist_grid(0, L - 1);
	std::uniform_real_distribution <double > dist_one(0.0, 1.0);
	const int buffer_offset = get_exp_buffer_offset(physics.J);
	int flipi, flipj;
	int dE;
	for (int i = 0; i < L * L * n; ++i) {
		flipi = dist_grid(rng);
		flipj = dist_grid(rng);
		dE = physics.J * get_dE(grid, flipi, flipj);
		if (dE <= 0 || (dist_one(rng) < exp_values[dE + buffer_offset]))
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
