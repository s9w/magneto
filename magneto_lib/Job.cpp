#include "Job.h"
#include "nlohmann/json.hpp"
#include "file_tools.h"
#include "logging.h"

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

} // namespace {}

void magneto::from_json(const nlohmann::json& j, magneto::Job& job) {
   set_enum_from_key(j, job.m_spin_start_mode, "spin_start", {"random", "ones", "image"});
   set_enum_from_key(j, job.m_temp_mode, "temp", { "single", "range", "image" });
   set_enum_from_key(j, job.m_algorithm, "algorithm", { "metropolis", "SW" });
   set_enum_from_key(j, job.m_image_mode.m_mode, "image_output_mode", { "none", "endimage", "intervals", "movie" });
   write_value_from_json(j, "t_min", job.m_t_min);
   write_value_from_json(j, "t_max", job.m_t_max);
   write_value_from_json(j, "t", job.m_t_single);
   write_value_from_json(j, "t_steps", job.m_temp_steps);
   write_value_from_json(j, "start_runs", job.m_start_runs);
   write_value_from_json(j, "L", job.m_L);
   write_value_from_json(j, "iterations", job.m_n);
   write_value_from_json(j, "spin_start_image_path", job.m_spin_start_image_path);
   write_value_from_json(j, "image_intervals", job.m_image_mode.m_intervals);
   write_value_from_json(j, "image_path", job.m_image_mode.m_path);
   write_value_from_json(j, "fps", job.m_image_mode.m_fps);
   write_value_from_json(j, "physics_path", job.m_physics_config.m_outputfile);
   write_value_from_json(j, "physics_format", job.m_physics_config.m_format);
}

magneto::Job magneto::get_parsed_job(const std::string& file_contents){
   nlohmann::json json = nlohmann::json::parse(file_contents);
   return json.get<Job>();
}


std::optional<magneto::Job> magneto::get_parsed_job(const std::filesystem::path& path){
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


bool magneto::operator==(const Job& a, const Job& b) {
   // use the std::tie trick for most
   if (std::tie(a.m_spin_start_mode, a.m_spin_start_image_path, a.m_temp_mode
         , a.m_temp_steps, a.m_start_runs
         , a.m_L, a.m_n, a.m_algorithm, a.m_image_mode, a.m_physics_config)
      !=
      std::tie(b.m_spin_start_mode, b.m_spin_start_image_path, b.m_temp_mode
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
