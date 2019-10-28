#include "system.h"
#include "Output.h"
#include "LatticeIndexRng.h"

#include <string>
#include <iostream>
#include <thread>
#include <mutex>
#include <map>

std::mutex rng_buffers_mutex;


std::vector<double> get_buffer_vector(const int buffer_size) {
   long long int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
   std::mt19937 rng = std::mt19937(seed1);
   std::uniform_real_distribution <double > dist_one(0.0, 1.0);
   std::vector<double> buffer(buffer_size);

   for (int i = 0; i < buffer_size; ++i)
      buffer[i] = dist_one(rng);
   return buffer;
}


void make_rng_numbers(
	std::map<std::thread::id, std::vector<double>>& rng_buffers, 
	const int buffer_size,
	std::map<std::thread::id, bool>& threads_done
) {
	std::pair<std::map<std::thread::id, std::vector<double>>::iterator, bool> emplace_pair;
   {
      std::scoped_lock guard(rng_buffers_mutex);
      emplace_pair = rng_buffers.emplace(std::this_thread::get_id(), std::vector<double>(buffer_size));
      threads_done.emplace(std::this_thread::get_id(), false);
   }

   long long int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
   std::mt19937 rng = std::mt19937(seed1);
	std::uniform_real_distribution <double > dist_one(0.0, 1.0);
   std::vector<double>& buffer = emplace_pair.first->second;
   for (int i = 0; i < buffer_size; ++i)
      buffer[i] = dist_one(rng);
   threads_done[std::this_thread::get_id()] = true;
}


void inner_loop(
   std::vector<std::thread>& rng_threads, 
   std::vector<double>& rng_front_buffer,
   std::map<std::thread::id, std::vector<double>>& rng_back_buffers,
   int& i_rng_buffer,
   const int n_rng_buffer,
   std::map<std::thread::id, bool>& threads_done
) {
   while (true) {
      for (std::thread& t : rng_threads) {
         const auto thread_id = t.get_id();
         if (!threads_done[thread_id])
            continue;
         
         t.join();
         {
            std::scoped_lock guard(rng_buffers_mutex);
            rng_front_buffer = rng_back_buffers.at(thread_id);
            rng_back_buffers.erase(thread_id);
            threads_done.erase(thread_id);
         }
         i_rng_buffer = 0;
         t = std::thread(make_rng_numbers, std::ref(rng_back_buffers), n_rng_buffer, std::ref(threads_done));
         return;
      }
   }
}


int main() {
	constexpr int L = 400;
	const int n_metropolis = 1;
	const magneto::PhysicsSettings physics{ 1, 2.2 };
	
	const int n_rng_buffer = L * L * 10;

	std::vector<double> rng_front_buffer = get_buffer_vector(n_rng_buffer);
   std::map<std::thread::id, std::vector<double>> rng_back_buffers;
	std::map<std::thread::id, bool> threads_done;

	const std::vector<double> exp_values = get_cached_exp_values(physics);

	const int max_rng_threads = 2;
   std::vector<std::thread> rng_threads;
	magneto::LatticeIndexRng lattice_rng(L);
	magneto::GridType grid = magneto::get_randomized_system(L);

	auto t0 = std::chrono::high_resolution_clock::now();

	magneto::Output writer(L, 8);

   // start two threads already
	for (int i = 0; i < 2; ++i)
		rng_threads.emplace_back(make_rng_numbers, std::ref(rng_back_buffers), n_rng_buffer, std::ref(threads_done));

	int i_rng_buffer = 0;
	for (int i = 1; i < 1000; ++i) {
		writer.photograph(grid);
		magneto::metropolis_sweeps(grid, lattice_rng, exp_values, rng_front_buffer, i_rng_buffer, physics, n_metropolis);
		i_rng_buffer++;
		if (i_rng_buffer >= rng_front_buffer.size() - 1 - L*L* n_metropolis) {
         inner_loop(rng_threads, rng_front_buffer, rng_back_buffers, i_rng_buffer, n_rng_buffer, threads_done);
		}
	}
   for (std::thread& t : rng_threads)
	   t.join();
	writer.make_movie();
	auto t1 = std::chrono::high_resolution_clock::now();
	float secs = (std::chrono::duration_cast <std::chrono::milliseconds> (t1 - t0).count()) * 0.001f;
	std::cout << "runtime: " << secs << std::endl;

	return 0;
}
