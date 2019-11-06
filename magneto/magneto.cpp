#include "system.h"
#include "Output.h"
#include "ProgressIndicator.h"
#include "windows.h"
#include "LatticeAlgorithms.h"

#include <execution>
#include <fstream>
#include <numeric>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/async.h"


std::shared_ptr<spdlog::logger> logger;


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
   double temp;
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


PhysicsResult get_physical_results(const double T, const int n, const int warmup_runs, const int L) {
	magneto::IsingSystem system(1, T, L);
	std::vector<magneto::PhysicalProperties> properties;
   for (int i = 1; i < warmup_runs; ++i) {
      system.wang_sweeps();
   }

	for (int i = 1; i < n; ++i) {
		properties.emplace_back(get_properties(system));
		system.wang_sweeps();
	}
   const magneto::PhysicalProperties mean_properties = get_mean(properties);
   const double e_mean = mean_properties.energy;
   const double e_sq_mean = mean_properties.energy_squared;
   const double m_mean = mean_properties.magnetization;
   const double m_sq_mean = mean_properties.magnetization_sq;
   const double cv = (e_sq_mean - e_mean * e_mean) * L * L / (T * T);
   const double chi = (m_sq_mean - m_mean * m_mean) * L * L / T;
   return {T, e_mean, cv, m_mean, chi };
}


/// <summary>Returns vector of n equidistant temperatures</summary>
std::vector<double> get_temps(const double tmin, const double tmax, const int n) {
   std::vector<double> temps;
   double temperature = tmin;
   const double temperature_step = (tmax - tmin) / n;
   for (int i = 0; i < n + 1; ++i) {
      temps.emplace_back(temperature);
      temperature += temperature_step;
   }
   return temps;
}


/// <summary>Returns theoretical value for critical temperature in 2D Ising model</summary>
double get_Tc() {
   return 2.0 / (log(1.0 + sqrt(2.0)));
}


void write_results(const std::vector<PhysicsResult>& results, const std::filesystem::path& path) {
   std::ofstream myfile(path);
   if (!myfile.is_open()) {
      logger->error("Couldn't open file {} for writing result.", path.string());
      return;
   }

   for(const PhysicsResult& result : results){
      myfile << fmt::format("{} \t {} \t {} \t {} \t {} \n", result.temp, result.energy, result.cv, result.magnetization, result.chi);
   }
   myfile.close();
}


void do_physics(
   const double tmin, const double tmax, const int steps,
   const int L, 
   const int n_warmup, const int n,
   const std::filesystem::path& path
) {
   if (tmin < 0) {
      logger->error("Negative temperature: Tmin={}", tmin);
      return;
   }
   if (tmax > tmin) {
      logger->error("TMax < TMin");
      return;
   }
   if (steps < 1) {
      logger->error("Too few temperature steps ({})", steps);
      return;
   }
   if (L < 1) {
      logger->error("System size invalid ({})", L);
      return;
   }
   const auto temps = get_temps(tmin, tmax, steps);

   std::vector<PhysicsResult> results(temps.size());
   std::transform(
      std::execution::par_unseq,
      std::cbegin(temps),
      std::cend(temps),
      std::begin(results),
      [&](const double t) {return get_physical_results(t, n, n_warmup, L); }
   );
   write_results(results, path);
}


void do_movie(const double t, const int L, const int J, const int n) {
   ProgressIndicator progress;
   magneto::IsingSystem system(J, t, L);
   //magneto::IsingSystem system(1, "input_random.png", "input_bunny.png");
   std::vector<magneto::PhysicalProperties> properties;
   
   magneto::Output writer(L, 1);
   magneto::Metropolis metropolis_alg(J, t, L, 2);
   for (int i = 1; i < n; ++i) {
   	writer.photograph(system.get_lattice());
      properties.emplace_back(get_properties(system));

      metropolis_alg.run(system.get_lattice_nc());
   	progress.set_progress(i, n);
   	progress.write_progress();
   }

   writer.make_movie();
}


int main() {
   logger = spdlog::basic_logger_mt< spdlog::async_factory>("basic_logger", "log.txt", true);
   logger->set_level(spdlog::level::info);
   set_console_cursor_visibility(false);

   //do_physics(0.01, 2.0*get_Tc(), 15, 200, 50, 300, "results.txt");
   do_movie(2.26, 500, 1, 30*10);

	return 0;
}
