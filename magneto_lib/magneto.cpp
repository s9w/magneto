#include "magneto.h"

#include "Job.h"
#include "IsingSystem.h"
#include "VisualOutput.h"
#include "ProgressIndicator.h"
#include "windows.h"
#include "LatticeAlgorithms.h"
#include "file_tools.h"
#include "physics_tools.h"

#include <execution>
#include <sstream>

#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/async.h>


/// <summary>Self-explanatory, but doesn't seem to work on powershell</summary>
void set_console_cursor_visibility(bool visibility) {
	HANDLE out = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_CURSOR_INFO     cursorInfo;
	GetConsoleCursorInfo(out, &cursorInfo);
	cursorInfo.bVisible = visibility; // set the cursor visibility
	cursorInfo.dwSize = 100;
	SetConsoleCursorInfo(out, &cursorInfo);
}


template<class TAlg, class TImager>
magneto::PhysicsResult get_physical_results(const double T, const magneto::Job& job) {
   spdlog::get("magneto_logger")->info("Starting computations for T={:<4.3f}, L={}", T, job.m_L);
   const int J = 1;
	magneto::IsingSystem system(J, T, job.m_L);

   // Initial warmup runs to bring the system into a realistic state
   {
      magneto::SW wang(J, T, job.m_L);
      for (unsigned int i = 1; i < job.m_start_runs; ++i) {
         wang.run(system.get_lattice_nc());
      }
   }

   // Main iterations
   TImager image_writer(job.m_L, job.m_image_mode, T);
   TAlg algorithm(J, T, job.m_L);
   std::vector<magneto::PhysicalProperties> properties;
	for (unsigned int i = 1; i < job.m_n; ++i) {
      image_writer.snapshot(system.get_lattice());
		properties.emplace_back(get_properties(system));
      algorithm.run(system.get_lattice_nc());
	}
   image_writer.end_actions();

   // compute results
   magneto::get_logger()->info("Finished computations for T={:<4.3f}, L={}", T, job.m_L);
   return get_physical_results(properties, job.m_L, T);
}


/// <summary>Returns vector of n equidistant temperatures</summary>
std::vector<double> get_temps(const double tmin, const double tmax, const int n) {
   std::vector<double> temps;
   double temperature = tmin;
   const double temperature_step = (tmax - tmin) / (n-1);
   for (int i = 0; i < n; ++i) {
      temps.emplace_back(temperature);
      temperature += temperature_step;
   }
   return temps;
}


void write_results(const std::vector<magneto::PhysicsResult>& results, const magneto::PhysicsConfig& physics_config) {
   std::string file_content;

   for (const magneto::PhysicsResult& result : results) {
      try {
         file_content += fmt::format(physics_config.m_format + "\n"
            , fmt::arg("T", result.temp)
            , fmt::arg("E", result.energy)
            , fmt::arg("cv", result.cv)
            , fmt::arg("M", result.magnetization)
            , fmt::arg("chi", result.chi)
         );
      }
      catch (const fmt::format_error& /*e*/) {
         magneto::get_logger()->error("Formatting string could not be parsed. Not writing results.");
      }
   }
   magneto::write_string_to_file(physics_config.m_outputfile, file_content);
}


std::shared_ptr<spdlog::logger> magneto::get_logger(std::vector<spdlog::sink_ptr> sinks) {
   const std::string logger_name = "magneto_logger";
   auto logger = spdlog::get(logger_name);
   if (!logger){
      if (!sinks.empty()){
         logger = std::make_shared<spdlog::logger>(logger_name,
            std::begin(sinks),
            std::end(sinks));
         spdlog::register_logger(logger);
      }
      else{
         std::vector<spdlog::sink_ptr> sinks;
         sinks.push_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
         sinks.push_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>("log.txt", true));
         logger = std::make_shared<spdlog::logger>(logger_name,
            std::begin(sinks),
            std::end(sinks));
         spdlog::register_logger(logger);
      }
   }
   logger->set_level(spdlog::level::info);
   spdlog::set_pattern("[%T] %^%v%$");
   return logger;
}


template<class TAlg, class TImager>
void run_job(const std::vector<double>& temps, const magneto::Job& job) {
   std::vector<magneto::PhysicsResult> results(temps.size());
   std::transform(
      std::execution::par_unseq,
      std::cbegin(temps),
      std::cend(temps),
      std::begin(results),
      [&](const double t) {return get_physical_results<TAlg, TImager>(t, job); }
   );
   write_results(results, job.m_physics_config);
}


template<class TAlg, class TImager>
void run_job(const magneto::Job& job) {
   if(job.m_temp_mode == magneto::TempStartMode::Image)
      magneto::get_logger()->error("Image-based temperatures currently not implemented");
   else if (job.m_temp_mode == magneto::TempStartMode::Many) {
      const auto stepped_temps = get_temps(job.m_t_min, job.m_t_max, job.m_temp_steps);
      run_job<TAlg, TImager>(stepped_temps, job);
   }
   else if (job.m_temp_mode == magneto::TempStartMode::Single)
      run_job<TAlg, TImager>({ job.m_t_single }, job);
}


template<class TAlg>
void run_job(const magneto::Job& job) {
   if (job.m_image_mode.m_mode == magneto::ImageOrMovie::Movie)
      run_job<magneto::Metropolis, magneto::MovieWriter>(job);
   else if (job.m_image_mode.m_mode == magneto::ImageOrMovie::Intervals)
      run_job<magneto::Metropolis, magneto::IntervalWriter>(job);
}


void run_job(const magneto::Job& job) {
   if (job.m_algorithm == magneto::Algorithm::Metropolis)
      run_job<magneto::Metropolis>(job);
   else
      run_job<magneto::SW>(job);
}


void magneto::start() {
   get_logger();
   set_console_cursor_visibility(false);

   const std::filesystem::path default_config_path = "magneto_config.json";
   const std::optional<Job> job = get_parsed_job(default_config_path);
   if (!job.has_value()) {
      get_logger()->error("No configuration file found at {}", default_config_path.string());
      return;
   }

   run_job(job.value());
}
