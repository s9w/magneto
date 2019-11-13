#include "magneto.h"

#include "Job.h"
#include "IsingSystem.h"
#include "VisualOutput.h"
#include "ProgressIndicator.h"
#include "windows.h"
#include "LatticeAlgorithms.h"
#include "file_tools.h"
#include "physics_tools.h"
#include "logging.h"

#include <execution>
#include <sstream>


/// <summary>Self-explanatory, but doesn't seem to work on powershell</summary>
void set_console_cursor_visibility(bool visibility) {
	HANDLE out = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_CURSOR_INFO     cursorInfo;
	GetConsoleCursorInfo(out, &cursorInfo);
	cursorInfo.bVisible = visibility; // set the cursor visibility
	cursorInfo.dwSize = 100;
	SetConsoleCursorInfo(out, &cursorInfo);
}


std::string get_temperature_string(const double T) {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(3) << T;
   std::string temp_string = stream.str();
   return stream.str();
}


std::string get_temperature_string(const magneto::LatticeDType& T) {
   return "image_temps";
}


std::unique_ptr<magneto::VisualOutput> get_visual_output(
   const magneto::ImageOrMovie& mode
   , const unsigned int L
   , const magneto::ImageMode& image_mode,
   const std::string& temp_string
) {
   if (mode == magneto::ImageOrMovie::Movie)
      return std::make_unique<magneto::MovieWriter>(L, image_mode, temp_string);
   else if (mode == magneto::ImageOrMovie::Intervals)
      return std::make_unique<magneto::IntervalWriter>(L, image_mode, temp_string);
   else if (mode == magneto::ImageOrMovie::Endimage)
      return std::make_unique<magneto::EndImageWriter>(L, image_mode, temp_string);
   else
      return std::make_unique<magneto::NullImageWriter>(L, image_mode, temp_string);
}


std::unique_ptr<magneto::LatticeAlgorithm> get_lattice_algorithm(
   const magneto::Algorithm& alg, 
   const magneto::LatticeDType& lattice_temps,
   const int L
) {
   const int J = 1;
   if (alg == magneto::Algorithm::Metropolis) {
         return std::make_unique<magneto::VariableMetropolis>(J, lattice_temps, L);
   }
   else {
      return std::make_unique<magneto::VariableSW>(J, lattice_temps, L);
   }
}


std::unique_ptr<magneto::LatticeAlgorithm> get_lattice_algorithm(
   const magneto::Algorithm& alg,
   const double T,
   const int L
) {
   const int J = 1;
   if (alg == magneto::Algorithm::Metropolis) {
      return std::make_unique<magneto::Metropolis>(J, T, L);
   }
   else {
      return std::make_unique<magneto::SW>(J, T, L);
   }
}


template<class TTemp>
void warmup_system(magneto::IsingSystem& system, const TTemp& T, const unsigned int runs) {
   const int L = static_cast<int>(system.get_lattice().size());
   auto alg = get_lattice_algorithm(magneto::Algorithm::SW, T, L);
   for (unsigned int i = 1; i < runs; ++i) {
      alg->run(system.get_lattice_nc());
   }
}


double get_t_representation_for_measurements(const double t) {
   return t;
}

double get_t_representation_for_measurements(const magneto::LatticeDType& /*t*/) {
   return 1.0;
}


template<class TTemp>
magneto::PhysicalProperties get_physical_properties(
   const TTemp T, 
   const magneto::Job& job
) {
   const std::string temp_string = get_temperature_string(T);
   std::unique_ptr<magneto::VisualOutput> visual_output(get_visual_output(job.m_image_mode.m_mode, job.m_L, job.m_image_mode, temp_string));
   std::unique_ptr<magneto::LatticeAlgorithm> algorithm(get_lattice_algorithm(job.m_algorithm, T, job.m_L));

   magneto::get_logger()->info("Starting computations for T={}, L={}", temp_string, job.m_L);
   const int J = 1;
	magneto::IsingSystem system(J, job.initial_spins);

   warmup_system(system, T, job.m_start_runs);

   // Main iterations
   std::vector<magneto::PhysicalMeasurement> measurements;
	for (unsigned int i = 1; i < job.m_n; ++i) {
      visual_output->snapshot(system.get_lattice());
		measurements.emplace_back(get_properties(system));
      algorithm->run(system.get_lattice_nc());
	}
   visual_output->snapshot(system.get_lattice(), true);
   visual_output->end_actions();

   // compute results
   magneto::get_logger()->info("Finished computations for T={}, L={}", temp_string, job.m_L);
   magneto::PhysicalProperties props;
   props.measurements = measurements;
   return { measurements, get_t_representation_for_measurements(T), job.m_L };
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


std::vector<magneto::PhysicalProperties> run_job_fixed_t(const magneto::Job& job, const std::vector<double>& temps) {
   std::vector<magneto::PhysicalProperties> properties(temps.size());
   std::transform(
      std::execution::par_unseq,
      std::cbegin(temps),
      std::cend(temps),
      std::begin(properties),
      [&](const double t) {return get_physical_properties(t, job); }
   );
   return properties;
}


void run_job(const magneto::Job& job, const std::variant<magneto::LatticeDType, std::vector<double>>& temp_variant) {
   struct V {
      V(const magneto::Job& job) : m_job(job) { }
      void operator()(const magneto::LatticeDType& T) {
         [[maybe_unused]] const magneto::PhysicalProperties properties = get_physical_properties(T, m_job);
      }
      void operator()(const std::vector<double>& T) {
         const std::vector<magneto::PhysicalProperties> properties = run_job_fixed_t(m_job, T);
         std::vector<magneto::PhysicsResult> results;
         for (const magneto::PhysicalProperties& prop : properties) {
            results.emplace_back(magneto::get_physical_results(prop));
         }
         write_results(results, m_job.m_physics_config);
      }
      magneto::Job m_job;
   };
   std::visit(V(job), temp_variant);
}


void magneto::start() {
   get_logger();
   set_console_cursor_visibility(false);

   const std::filesystem::path default_config_path = "magneto_config.json";
   const std::optional<JsonJob> parsed_job = get_parsed_job(default_config_path);
   if (!parsed_job.has_value()) {
      get_logger()->error("No configuration file found at {}", default_config_path.string());
      return;
   }

   const auto [job, T] = get_job(parsed_job.value());
   run_job(job, T);
}
