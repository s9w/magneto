#include <chrono>
#include <random>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "system.h"

namespace {
	unsigned char get_png_value_from_pm_one(const int value) {
		return (value + 1) / 2 * 255;
	}
}


GridType get_randomized_system(const int L){
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


void write_png(const GridType& grid, const std::filesystem::path& path){
	std::vector<unsigned char> grid_png;
	const int L = static_cast<int>(grid.size());
	grid_png.reserve(L * L);
	for (const auto& row : grid)
		for (const auto& elem : row)
			grid_png.emplace_back(get_png_value_from_pm_one(elem));

	stbi_write_png("test.png", L, L, 1, grid_png.data(), L * 1);
}
