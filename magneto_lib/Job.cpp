#include "Job.h"
#include "file_tools.h"
#include "logging.h"
#include "IsingSystem.h"
#include "nlohmann/json.hpp"
#include <random>


namespace {
   bool is_equal(double a, double b) {
      constexpr double tol = 0.001;
      return fabs(a - b) < tol;
   }


   bool has_ending(std::string const& fullString, std::string const& ending) {
      if (fullString.length() >= ending.length()) {
         return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
      }
      else {
         return false;
      }
   }


   template<class T>
   void set_enum_from_key(const nlohmann::json& j, T& target_enum, const char* key, const std::vector<std::string>& associations) {
      if (!j.contains(key))
         return;
      std::string value = j.at(key).get<std::string>();
      for (int i = 0; i < associations.size(); ++i) {
         if(!(associations[i] == value))
            continue;
         target_enum = static_cast<T>(i);
         return;
      }
      magneto::get_logger()->warn("Value {} invalid for config file setting {}.", value, key);
   }


   template<class T>
   void write_value_from_json(const nlohmann::json& j, const char* key, T& target) {
      if (j.contains(key))
         j.at(key).get_to(target);
   }

   template<>
   void write_value_from_json(const nlohmann::json& j, const char* key, std::filesystem::path& target) {
      if (j.contains(key))
         target = j.at(key).get<std::string>();
   }


   /// <summary>Returns vector of n equidistant temperatures</summary>
   std::vector<double> get_temps(const magneto::TempStartMode& mode, const double tmin, const double tmax, const int n) {
      if (mode == magneto::TempStartMode::Single)
         return { tmin };
      std::vector<double> temps;
      double temperature = tmin;
      const double temperature_step = (tmax - tmin) / (n - 1);
      for (int i = 0; i < n; ++i) {
         temps.emplace_back(temperature);
         temperature += temperature_step;
      }
      return temps;
   }

} // namespace {}

void magneto::from_json(const nlohmann::json& j, magneto::JsonJob& job) {
   set_enum_from_key(j, job.m_spin_start_mode, "spin_start", {"random", "image"});
   set_enum_from_key(j, job.m_temp_mode, "temp", { "single", "range", "image" });
   set_enum_from_key(j, job.m_algorithm, "algorithm", { "metropolis", "SW" });
   set_enum_from_key(j, job.m_image_mode.m_mode, "image_output_mode", { "none", "endimage", "intervals", "movie" });
   write_value_from_json(j, "t_min", job.m_t_min);
   write_value_from_json(j, "t_max", job.m_t_max);
   write_value_from_json(j, "t", job.m_t_single);
   write_value_from_json(j, "t_steps", job.m_temp_steps);
   write_value_from_json(j, "t_image", job.m_temperature_image);
   write_value_from_json(j, "start_runs", job.m_start_runs);
   write_value_from_json(j, "L", job.m_L);
   write_value_from_json(j, "J", job.m_J);
   write_value_from_json(j, "iterations", job.m_n);
   write_value_from_json(j, "spin_start_image_path", job.m_spin_start_image_path);
   write_value_from_json(j, "image_intervals", job.m_image_mode.m_intervals);
   write_value_from_json(j, "image_path", job.m_image_mode.m_path);
   write_value_from_json(j, "fps", job.m_image_mode.m_fps);
   write_value_from_json(j, "physics_path", job.m_physics_config.m_outputfile);
   write_value_from_json(j, "physics_format", job.m_physics_config.m_format);
}

magneto::JsonJob magneto::get_parsed_job(const std::string& file_contents){
   nlohmann::json json = nlohmann::json::parse(file_contents);
   magneto::JsonJob job(json.get<JsonJob>());
   job.m_temperatures = get_temps(job.m_temp_mode, job.m_t_min, job.m_t_max, job.m_temp_steps);
   return job;
}


//magneto::Job
std::tuple<magneto::Job, std::variant<magneto::LatticeDType, std::vector<double>>>
magneto::get_job(const JsonJob& json_job){
   Job job;
   job.m_L = json_job.m_L;
   
   if (json_job.m_spin_start_mode == SpinStartMode::Random)
      job.initial_spins = get_randomized_system(json_job.m_L);
   else //SpinStartMode::Image
      job.initial_spins = get_lattice_from_png_file(json_job.m_spin_start_image_path);

   job.m_algorithm = json_job.m_algorithm;
   job.m_n = json_job.m_n;
   job.m_J = json_job.m_J;
   job.m_image_mode = json_job.m_image_mode;
   job.m_physics_config = json_job.m_physics_config;


   std::variant<magneto::LatticeDType, std::vector<double>> T;
   if (json_job.m_temp_mode == TempStartMode::Image)
      T = get_lattice_temps_from_png_file(json_job.m_temperature_image, json_job.m_t_min, json_job.m_t_max);
   else if (json_job.m_temp_mode == TempStartMode::Single)
      T = std::vector<double>{ json_job.m_t_min };
   else
      T = get_temps(json_job.m_temp_mode, json_job.m_t_min, json_job.m_t_max, json_job.m_temp_steps);

   return { job, T };
}


std::optional<magneto::JsonJob> magneto::get_parsed_job(const std::filesystem::path& path){
   const auto file_contents = get_file_contents(path);
   if (!file_contents.has_value())
      return std::nullopt;
   return get_parsed_job(file_contents.value());
}

bool magneto::operator==(const ImageMode& a, const ImageMode& b) {
   return std::tie(a.m_fps, a.m_intervals, a.m_mode, a.m_path) ==
      std::tie(b.m_fps, b.m_intervals, b.m_mode, b.m_path);
}
bool magneto::operator==(const PhysicsConfig& a, const PhysicsConfig& b) {
   return std::tie(a.m_outputfile, a.m_format) == std::tie(b.m_outputfile, b.m_format);
}


bool magneto::operator==(const JsonJob& a, const JsonJob& b) {
   // use the std::tie trick for most
   if (std::tie(a.m_spin_start_mode, a.m_spin_start_image_path, a.m_temperature_image, a.m_temp_mode
         , a.m_temp_steps, a.m_start_runs
         , a.m_L, a.m_n, a.m_algorithm, a.m_image_mode, a.m_physics_config)
      !=
      std::tie(b.m_spin_start_mode, b.m_spin_start_image_path, a.m_temperature_image, b.m_temp_mode
         , b.m_temp_steps, b.m_start_runs
         , b.m_L, b.m_n, b.m_algorithm, b.m_image_mode, b.m_physics_config))
   {
      return false;
   }

   // Handle doubles separately because of the floating point comparisons
   if (!(is_equal(a.m_t_single, b.m_t_single)))
      return false;
   if (!(is_equal(a.m_t_min, b.m_t_min)))
      return false;
   if (!(is_equal(a.m_t_max, b.m_t_max)))
      return false;
   return true;
}
