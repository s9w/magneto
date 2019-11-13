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
   void set_enum_from_key(
      const nlohmann::json& j, T& target_enum, const char* key, const std::vector<std::string>& associations
   ) {
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
   std::vector<double> get_temps(
      const magneto::TempStartMode& mode, const double tmin, const double tmax, const int n
   ) {
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


   std::variant<magneto::LatticeDType, std::vector<double>>
   get_temp_variant(const magneto::JsonJob& job) {
      std::variant<magneto::LatticeDType, std::vector<double>> T;
      if (job.temp_mode == magneto::TempStartMode::Image)
         return magneto::get_lattice_temps_from_png_file(job.temperature_image, job.t_min, job.t_max);
      else if (job.temp_mode == magneto::TempStartMode::Single)
         return std::vector<double>{ job.t_min };
      else
         return get_temps(job.temp_mode, job.t_min, job.t_max, job.temp_steps);

   }


   /// <summary>Returns target system dimensions. Uses all the data supplied from the job and
   /// adheres to the prority defined in the documentation.</summary>
   std::pair<unsigned int, unsigned int> get_system_size(
      const magneto::JsonJob& json_job,
      const std::optional<magneto::LatticeType>& image_spin_state,
      const std::variant<magneto::LatticeDType, std::vector<double>>& t_variant
   ) {
      // User-supplied Lx and Ly
      if (json_job.Lx > 0 && json_job.Ly > 0)
         return { json_job.Lx, json_job.Ly };

      // User-supplied L
      if (json_job.L > 0)
         return { json_job.L, json_job.L };

      // Size taken from spin state image
      if (image_spin_state.has_value()) {
         return magneto::get_dimensions_of_lattice(image_spin_state.value());
      }

      // Size taken from temperature image
      if (std::holds_alternative<magneto::LatticeDType>(t_variant)) {
         return magneto::get_dimensions_of_lattice(std::get<magneto::LatticeDType>(t_variant));
      }

      // Fallback-default
      const unsigned int default_Lx = 500;
      const unsigned int default_Ly = 500;
      return { default_Lx , default_Ly };
   }


   std::optional<magneto::LatticeType> get_spin_state_from_job(const magneto::JsonJob& json_job) {
      if (json_job.spin_start_mode != magneto::SpinStartMode::Image)
         return std::nullopt;
      return magneto::get_spin_state_from_png(json_job.spin_start_image_path);
   }


} // namespace {}


void magneto::from_json(const nlohmann::json& j, magneto::JsonJob& job) {
   set_enum_from_key(j, job.spin_start_mode, "spin_start", {"random", "image"});
   set_enum_from_key(j, job.temp_mode, "temp", { "single", "range", "image" });
   set_enum_from_key(j, job.algorithm, "algorithm", { "metropolis", "SW" });
   set_enum_from_key(j, job.image_mode.m_mode, "image_output_mode", { "none", "endimage", "intervals", "movie" });
   write_value_from_json(j, "t_min", job.t_min);
   write_value_from_json(j, "t_max", job.t_max);
   write_value_from_json(j, "t", job.t_single);
   write_value_from_json(j, "t_steps", job.temp_steps);
   write_value_from_json(j, "t_image", job.temperature_image);
   write_value_from_json(j, "start_runs", job.start_runs);
   write_value_from_json(j, "L", job.L);
   write_value_from_json(j, "Lx", job.Lx);
   write_value_from_json(j, "Ly", job.Ly);
   write_value_from_json(j, "J", job.J);
   write_value_from_json(j, "iterations", job.n);
   write_value_from_json(j, "spin_start_image_path", job.spin_start_image_path);
   write_value_from_json(j, "image_intervals", job.image_mode.m_intervals);
   write_value_from_json(j, "image_path", job.image_mode.m_path);
   write_value_from_json(j, "fps", job.image_mode.m_fps);
   write_value_from_json(j, "physics_path", job.physics_config.m_outputfile);
   write_value_from_json(j, "physics_format", job.physics_config.m_format);
}


magneto::JsonJob magneto::get_parsed_job(const std::string& file_contents){
   nlohmann::json json = nlohmann::json::parse(file_contents);
   magneto::JsonJob job(json.get<JsonJob>());
   return job;
}


std::optional<magneto::LatticeType> get_corrected_image_spin_state(
   std::optional<magneto::LatticeType> image_spin_state,
   unsigned int Lx, 
   unsigned int Ly,
   const magneto::JsonJob& json_job
)
{
   if (!image_spin_state.has_value())
      return image_spin_state.value();
   const auto [x, y] = magneto::get_dimensions_of_lattice(image_spin_state.value());
   if(x == Ly && y == Ly)
      return image_spin_state.value();

   // There is an image-based spin state and its size is wrong
   return magneto::get_resized_data<magneto::LatticeType>(
      json_job.spin_start_image_path,
      Lx,
      Ly,
      [](const std::filesystem::path path) {return magneto::get_spin_state_from_png(path).value(); }
   );
}


std::variant<magneto::LatticeDType, std::vector<double>>
get_corrected_temp_state(
   const std::variant<magneto::LatticeDType, std::vector<double>>& t_variant,
   unsigned int Lx,
   unsigned int Ly,
   const magneto::JsonJob& job
) {
   if (!std::holds_alternative<magneto::LatticeDType>(t_variant))
      return t_variant;
   const auto [x, y] = magneto::get_dimensions_of_lattice(std::get<magneto::LatticeDType>(t_variant));
   if (x == Ly && y == Ly)
      return t_variant;

   // There is an image-based temp state and its size is wrong
   //return magneto::get_lattice_temps_from_png_file(job.temperature_image, job.t_min, job.t_max, Lx, Ly);
   return magneto::get_resized_data<magneto::LatticeDType>(
      job.temperature_image,
      Lx,
      Ly,
      [&](const std::filesystem::path path) {return magneto::get_lattice_temps_from_png_file(path, job.t_min, job.t_max); }
   );
}


//magneto::Job
std::tuple<magneto::Job, std::variant<magneto::LatticeDType, std::vector<double>>>
magneto::get_job(const JsonJob& json_job){
   Job job;

   // A word about System size: The user can supply it explicitly, he can also supply images that
   // can be used to implicitly infer the system size. Trouble arises when the explicit size doesn't
   // match the image sizes. Then, these are being resized to the system size and read again.
   std::optional<magneto::LatticeType> image_spin_state = get_spin_state_from_job(json_job);
   auto t = get_temp_variant(json_job);

   std::tie(job.m_Lx, job.m_Ly) = get_system_size(json_job, image_spin_state, t);

   // Resize image-based data if necessary
   image_spin_state = get_corrected_image_spin_state(image_spin_state, job.m_Lx, job.m_Ly, json_job);
   t = get_corrected_temp_state(t, job.m_Lx, job.m_Ly, json_job);
   
   if (json_job.spin_start_mode == SpinStartMode::Random)
      job.initial_spins = get_randomized_system(job.m_Lx, job.m_Ly);
   else //SpinStartMode::Image
      job.initial_spins = image_spin_state.value();

   job.m_algorithm = json_job.algorithm;
   job.m_n = json_job.n;
   job.m_start_runs = json_job.start_runs;
   job.m_J = json_job.J;
   job.m_image_mode = json_job.image_mode;
   job.m_physics_config = json_job.physics_config;

   return { job, t };
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
   if (std::tie(a.spin_start_mode, a.spin_start_image_path, a.temperature_image, a.temp_mode
         , a.temp_steps, a.start_runs
         , a.L, a.n, a.algorithm, a.image_mode, a.physics_config)
      !=
      std::tie(b.spin_start_mode, b.spin_start_image_path, a.temperature_image, b.temp_mode
         , b.temp_steps, b.start_runs
         , b.L, b.n, b.algorithm, b.image_mode, b.physics_config))
   {
      return false;
   }

   // Handle doubles separately because of the floating point comparisons
   if (!(is_equal(a.t_single, b.t_single)))
      return false;
   if (!(is_equal(a.t_min, b.t_min)))
      return false;
   if (!(is_equal(a.t_max, b.t_max)))
      return false;
   return true;
}
