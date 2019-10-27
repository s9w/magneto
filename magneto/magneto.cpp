#include "system.h"

#include <string>
#include <iostream>

int main() {
	std::vector<double> exp_values;
	const int J = 1;
	const double T = 2.2;
	long long int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 rng = std::mt19937(seed1);
	int buffer_offset = 8 * std::abs(J);
	for (int dE = 0; dE < buffer_offset * 2 + 1; ++dE)
		exp_values.push_back(exp(-(dE - buffer_offset) / T));

	constexpr int L = 50;
	GridType grid = get_randomized_system(L);

	write_png(grid, "test_0.png");
	for (int i = 1; i < 100; ++i) {
		metropolis_sweeps(grid, rng, exp_values);
		const std::string filename = std::string("test_") + std::to_string(i) + ".png";
		write_png(grid, filename);
	}

	

	return 0;
}
