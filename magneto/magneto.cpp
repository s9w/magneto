#include "system.h"
#include "Output.h"
#include "LatticeIndexRng.h"

#include <future>
#include <string>
#include <iostream>
#include <thread>


std::vector<double> get_random_buffer(
	const int buffer_size
) {
	long long int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 rng = std::mt19937(seed1);
	std::uniform_real_distribution <double > dist_one(0.0, 1.0);
	std::vector<double> normal_random_vector;
	normal_random_vector.reserve(buffer_size);
	for (int i = 0; i < buffer_size; ++i)
		normal_random_vector.emplace_back(dist_one(rng));
	return normal_random_vector;
}


void refill_exp_function_buffer(
	std::vector<std::future<std::vector<double>>>& futures,
	std::vector<double>& random_buffer,
	const int exp_function_buffer_size
) {
	// usually there should be a thread already finished
	for (auto& future : futures) {
		// future.is_ready() is ironically not yet in the standard, but it is in MSVC... and way too convenient to ignore
		if (future._Is_ready()) {
			random_buffer = std::move(future.get());
			future = std::async(std::launch::async, get_random_buffer, exp_function_buffer_size);
			return;
		}
	}
	// no thread finished yet.. should only happen at the start after the first sweep
	random_buffer = std::move(futures.front().get());
	futures.front() = std::async(std::launch::async, get_random_buffer, exp_function_buffer_size);
}


int main() {
	constexpr int L = 400;
	const magneto::PhysicsSettings physics{ 1, 2.2 };
	const int random_buffer_size = L * L;

	std::vector<double> random_buffer = get_random_buffer(random_buffer_size);
	const std::vector<double> exp_values = get_cached_exp_values(physics);

	magneto::Output writer(L, 8);
	const int max_rng_threads = 2; // two threads compting random values saturate one thread of metropolis sweeps
	std::vector<std::future<std::vector<double>>> future_random_buffers;
	magneto::LatticeIndexRng lattice_rng(L);
	magneto::GridType grid = magneto::get_randomized_system(L);

	auto t0 = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < max_rng_threads; ++i)
		future_random_buffers.emplace_back(std::async(std::launch::async, get_random_buffer, random_buffer_size));

	for (int i = 1; i < 1000; ++i) {
		writer.photograph(grid);
		magneto::metropolis_sweeps(grid, lattice_rng, exp_values, random_buffer, physics);
        refill_exp_function_buffer(future_random_buffers, random_buffer, random_buffer_size);
	}
	for (auto& future : future_random_buffers)
		future.get();
	writer.make_movie();
	auto t1 = std::chrono::high_resolution_clock::now();
	float secs = (std::chrono::duration_cast <std::chrono::milliseconds> (t1 - t0).count()) * 0.001f;
	std::cout << "runtime: " << secs << std::endl;

	return 0;
}
