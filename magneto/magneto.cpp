#include "system.h"
#include "Output.h"
#include "ProgressIndicator.h"

#include <future>
#include <string>
#include <iostream>
#include <thread>


std::vector<double> get_random_buffer(const size_t buffer_size) {
	unsigned int seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
	std::mt19937 rng = std::mt19937(seed);
	std::uniform_real_distribution <double > dist_one(0.0, 1.0);
	std::vector<double> normal_random_vector;
	normal_random_vector.reserve(buffer_size);
	for (int i = 0; i < buffer_size; ++i)
		normal_random_vector.emplace_back(dist_one(rng));
	return normal_random_vector;
}


IndexPairVector get_lattice_indices(const size_t buffer_size, const int lattice_size) {
	unsigned int seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
	std::mt19937 rng = std::mt19937(seed);
	std::uniform_int_distribution<> dist_lattice(0, lattice_size-1);
	IndexPairVector indices;
	indices.reserve(buffer_size);
	for (int i = 0; i < buffer_size; ++i)
		indices.emplace_back(dist_lattice(rng), dist_lattice(rng));
	return indices;
}


template <class T>
void assign_target_buffer_and_relaunch_future(
	T& target_buffer,
	std::future<T>& future,
	const std::function<T()>& thread_function
) {
	target_buffer = std::move(future.get());
	future = std::async(std::launch::async, thread_function);
}


/// <summary>Fills the target buffer with the result of any of the futures. If no future is finished,
/// it's being waited until one is. The future is restarted at the end.</summary>
template <class T>
void refill_target_buffer_from_futures(
	std::vector<std::future<T>>& futures,
	T& target_buffer,
	const std::function<T()>& thread_function
) {
	// usually there should be a thread already finished
	for (auto& future : futures) {
		// future.is_ready() is ironically not yet in the standard, but it is in MSVC... and way too convenient to ignore
		if (future._Is_ready()) {
			assign_target_buffer_and_relaunch_future(target_buffer, future, thread_function);
			return;
		}
	}
	// no thread finished yet.. should only happen at the start after the first sweep or when the system is too busy
	assign_target_buffer_and_relaunch_future(target_buffer, futures.front(), thread_function);
}


int main() {
	ProgressIndicator progress;
	magneto::IsingSystem system(1, 2.0, 500);
	const size_t random_buffer_size = system.get_L() * system.get_L();
	auto get_random_buffer_baked = [&] {return get_random_buffer(random_buffer_size); };
	std::vector<double> random_buffer = get_random_buffer_baked();
	auto get_lattice_indices_baked = [&] {return get_lattice_indices(random_buffer_size, static_cast<int>(system.get_L())); };
	IndexPairVector lattice_index_buffer = get_lattice_indices_baked();
	std::vector<std::future<IndexPairVector>> lattice_index_futures;

	magneto::Output writer(system.get_L(), 8);
	const int max_rng_threads = 2; // two threads computing random values saturate one thread of metropolis sweeps
	std::vector<std::future<std::vector<double>>> future_random_buffers;

	auto t0 = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < max_rng_threads; ++i) {
		future_random_buffers.emplace_back(std::async(std::launch::async, get_random_buffer_baked));
		lattice_index_futures.emplace_back(std::async(std::launch::async, get_lattice_indices_baked));
	}

	const int metropolis_sweets = 1000;
	for (int i = 1; i < metropolis_sweets; ++i) {
		writer.photograph(system.get_lattice());
		system.metropolis_sweeps(lattice_index_buffer, random_buffer);
		refill_target_buffer_from_futures<std::vector<double>>(future_random_buffers, random_buffer, get_random_buffer_baked);
		refill_target_buffer_from_futures<IndexPairVector>(lattice_index_futures, lattice_index_buffer, get_lattice_indices_baked);
		progress.set_progress(i, metropolis_sweets);
		progress.write_progress();
	}
	for (auto& future : future_random_buffers)
		future.get();
	writer.make_movie();

	auto t1 = std::chrono::high_resolution_clock::now();
	float secs = (std::chrono::duration_cast <std::chrono::milliseconds> (t1 - t0).count()) * 0.001f;
	std::cout << std::endl << "runtime: " << secs << std::endl;

	return 0;
}
