#include "system.h"

#include <string>
#include <iostream>


int main() {
	const magneto::PhysicsSettings physics{ 1, 2.2 };
	const std::vector<double> exp_values = get_cached_exp_values(physics);
	long long int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 rng = std::mt19937(seed1);

	constexpr int L = 1000;
	magneto::GridType grid = magneto::get_randomized_system(L);

	auto t0 = std::chrono::high_resolution_clock::now();

	magneto::write_png(grid, "test_0.png");
	for (int i = 1; i < 100; ++i) {
		magneto::metropolis_sweeps(grid, rng, exp_values, physics);
		const std::string filename = std::string("test_") + std::to_string(i) + ".png";
		magneto::write_png(grid, filename);
	}
	auto t1 = std::chrono::high_resolution_clock::now();
	float secs = (std::chrono::duration_cast <std::chrono::milliseconds> (t1 - t0).count()) * 0.001f;
	std::cout << "runtime: " << secs << std::endl;

	return 0;
}
