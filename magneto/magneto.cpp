#include "system.h"
#include "Output.h"
#include "ProgressIndicator.h"
#include "windows.h"

#include <future>
#include <string>
#include <iostream>
#include <sstream>
#include <thread>
#include <numeric>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/async.h"
#include <execution>

std::shared_ptr<spdlog::logger> logger;

std::string thread_id_to_string(const std::thread::id& id) {
   std::stringstream ss;
   ss << id;
   return ss.str();
}

template<class T>
struct ResultWithThreadId {
   operator T() {
      return m_result;
   }
   T m_result;
   std::string m_thread_id;
};

using UniformRandomBufferReturn = ResultWithThreadId<std::vector<double>>;
using IndexPairVectorReturn = ResultWithThreadId<IndexPairVector>;

UniformRandomBufferReturn get_random_buffer(const size_t buffer_size) {
   const std::string thread_id = thread_id_to_string(std::this_thread::get_id());
	unsigned int seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
	std::mt19937 rng = std::mt19937(seed);
	std::uniform_real_distribution <double > dist_one(0.0, 1.0);
	std::vector<double> normal_random_vector;
	normal_random_vector.reserve(buffer_size);
	for (int i = 0; i < buffer_size; ++i)
		normal_random_vector.emplace_back(dist_one(rng));
   logger->debug("get_random_buffer() done from thread {}", thread_id);
   return { normal_random_vector, thread_id };
}


IndexPairVectorReturn get_lattice_indices(const size_t buffer_size, const int lattice_size) {
   const std::string thread_id = thread_id_to_string(std::this_thread::get_id());
	unsigned int seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
	std::mt19937 rng = std::mt19937(seed);
	std::uniform_int_distribution<> dist_lattice(0, lattice_size-1);
	IndexPairVector indices;
	indices.reserve(buffer_size);
	for (int i = 0; i < buffer_size; ++i)
		indices.emplace_back(dist_lattice(rng), dist_lattice(rng));
   logger->debug("get_lattice_indices() done from thread {}", thread_id);
   return { indices, thread_id };
}


template <class T>
void assign_target_buffer_and_relaunch_future(
	T& target_buffer,
	std::future<ResultWithThreadId<T>>& future,
	const std::function<ResultWithThreadId<T>()>& thread_function
) {
   const auto fresult = future.get();
   logger->debug("fetched results and restarting. thread: {}", fresult.m_thread_id);
	target_buffer = std::move(fresult.m_result);
	future = std::async(std::launch::async, thread_function);
}


/// <summary>Fills the target buffer with the result of any of the futures. If no future is finished,
/// it's being waited until one is. The future is restarted at the end.</summary>
template <class T>
void refill_target_buffer_from_futures(
	std::vector<std::future<ResultWithThreadId<T>>>& futures,
	T& target_buffer,
	const std::function<ResultWithThreadId<T>()>& thread_function
) {
	// usually there should be a thread already finished
   while (true) {
      for (auto& future : futures) {
         // future.is_ready() is ironically not yet in the standard, but it is in MSVC... and way too convenient to ignore
         if (future._Is_ready()) {
            assign_target_buffer_and_relaunch_future(target_buffer, future, thread_function);
            return;
         }
      }
   }
}


/// <summary>Self-explanatory, but doesn't seem to work on powershell</summary>
void set_console_cursor_visibility(bool visibility) {
	HANDLE out = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_CURSOR_INFO     cursorInfo;
	GetConsoleCursorInfo(out, &cursorInfo);
	cursorInfo.bVisible = visibility; // set the cursor visibility
	cursorInfo.dwSize = 100;
	SetConsoleCursorInfo(out, &cursorInfo);
}


struct PhysicsResult {
   double energy;
   double cv;
   double magnetization;
   double chi;
};


template <class T>
T get_mean(const std::vector<T>& values) {
   const T sum = std::accumulate(values.cbegin(), values.cend(), T());
   return sum / static_cast<int>(values.size());
}


PhysicsResult get_physical_results(const double temperature, const int warmup_runs) {
	magneto::IsingSystem system(1, temperature, 200);
	std::vector<magneto::PhysicalProperties> properties;
   for (int i = 1; i < warmup_runs; ++i) {
      system.wang_sweeps();
   }

	const int metropolis_sweets = 200;
	for (int i = 1; i < metropolis_sweets; ++i) {
		properties.emplace_back(get_properties(system));
		system.wang_sweeps();
	}
   const magneto::PhysicalProperties mean_properties = get_mean(properties);
   const double e_mean = mean_properties.energy;
   const double e_sq_mean = mean_properties.energy_squared;
   const double m_mean = mean_properties.magnetization;
   const double m_sq_mean = mean_properties.magnetization_sq;
	const auto L = system.get_L();
   const double T = system.get_temp().value();
   const double cv = (e_sq_mean - e_mean * e_mean) * L * L / (T * T);
   const double chi = (m_sq_mean - m_mean * m_mean) * L * L / T;
   return {e_mean, cv, m_mean, chi };
}


int main() {
   std::vector<double> temps;
   const double Tc = 2.26;
   for (double t = 0.8 * Tc; t < 1.2 * Tc; t += 0.1)
      temps.emplace_back(t);

   std::vector<PhysicsResult> results(temps.size()); // needs to be correct size. .reserve() not enough
   std::transform(
      std::execution::par_unseq,
      std::cbegin(temps),
      std::cend(temps),
      std::begin(results),
      [](const double t) {return get_physical_results(t, 50); }
   );
   for (int i = 0; i < temps.size(); ++i)
      std::cout << fmt::format("T: {}, E: {}, cv: {}, m: {}, chi: {}", temps[i], results[i].energy, results[i].cv, results[i].magnetization, results[i].chi) << std::endl;

 //  logger = spdlog::basic_logger_mt< spdlog::async_factory>("basic_logger", "log.txt");
 //  logger->set_level(spdlog::level::info);

	//set_console_cursor_visibility(false);
	//ProgressIndicator progress;
	//magneto::IsingSystem system(1, 2.26, 600);
	////magneto::IsingSystem system(1, "input_random.png", "input_bunny.png");
 //  std::vector<magneto::PropertySnapshot> properties;

	//const size_t random_buffer_size = system.get_L() * system.get_L();
	//auto get_random_buffer_baked = [&] {return get_random_buffer(random_buffer_size); };
	//std::vector<double> random_buffer = get_random_buffer_baked();
	//auto get_lattice_indices_baked = [&] {return get_lattice_indices(random_buffer_size, static_cast<int>(system.get_L())); };
	//IndexPairVector lattice_index_buffer = get_lattice_indices_baked();
	//
	//magneto::Output writer(system.get_L(), 1);
	//const int max_rng_threads = 2; // two threads computing random values saturate one thread of metropolis sweeps
 //  std::vector<std::future<IndexPairVectorReturn>> lattice_index_futures;
	//std::vector<std::future<UniformRandomBufferReturn>> future_random_buffers;

	//auto t0 = std::chrono::high_resolution_clock::now();

	//for (int i = 0; i < max_rng_threads; ++i) {
	//	future_random_buffers.emplace_back(std::async(std::launch::async, get_random_buffer_baked));
	//	lattice_index_futures.emplace_back(std::async(std::launch::async, get_lattice_indices_baked));
	//}

	//const int metropolis_sweets = 100;
	//for (int i = 1; i < metropolis_sweets; ++i) {
	//	writer.photograph(system.get_lattice());
 //     properties.emplace_back(get_properties(system));
 //     logger->info("E: {}\t E^2: {}", properties.back().energy, properties.back().energy_squared);

	//	system.wang_sweeps();
	//	//system.metropolis_sweeps(lattice_index_buffer, random_buffer);
 // //    refill_target_buffer_from_futures<std::vector<double>>(future_random_buffers, random_buffer, get_random_buffer_baked);
	//	//refill_target_buffer_from_futures<IndexPairVector>(lattice_index_futures, lattice_index_buffer, get_lattice_indices_baked);
	//	progress.set_progress(i, metropolis_sweets);
	//	progress.write_progress();
	//}
	//const double e_mean = get_energy_mean(properties);
	//const double e_sq_mean = get_energy_squared_mean(properties);
	//const double cv = (e_sq_mean - e_mean * e_mean) / system.get_temp().value();
	//std::cout << fmt::format("e_mean: {}, e_sq_mean: {}, cv: {}", e_mean, e_sq_mean, cv) << std::endl;

	//for (auto& future : future_random_buffers)
	//	future.get();
	//writer.make_movie();


	return 0;
}
