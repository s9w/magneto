#include "Job.h"

#include "nlohmann/json.hpp"

namespace {
   bool is_equal(double a, double b) {
      constexpr double tol = 0.001;
      return fabs(a - b) < tol;
   }

} // namespace {}


void magneto::from_json(const nlohmann::json& j, magneto::Job& job) {
   // oh the joys of parsing json

   if (j.contains("spin_start")) {
      std::string spin_start = j.at("spin_start").get<std::string>();
      if (spin_start == "random")
         job.m_spinstart = RandomSpinStart();
      else if (spin_start == "ones")
         job.m_spinstart = OneSpinStart();
      else
         job.m_spinstart = ImageSpinStart{ spin_start };
   }

   if (j.contains("temp")) {
      auto x = j.at("temp");
      if (j.at("temp").is_number_float())
         job.m_temp = j.at("temp").get<double>();
      else {
         job.m_temp = ImageTempStart{ j.at("temp").get<std::string>() };
      }
   }

   if (j.contains("image_temp_min")) 
      j.at("image_temp_min").get_to(std::get<ImageTempStart>(job.m_temp).t_min);
   if (j.contains("image_temp_max"))
      j.at("image_temp_max").get_to(std::get<ImageTempStart>(job.m_temp).t_max);

   if (j.contains("start_runs"))
      j.at("start_runs").get_to(job.m_start_runs);
   if (j.contains("L"))
      j.at("L").get_to(job.m_L);
   if (j.contains("iterations"))
      j.at("iterations").get_to(job.m_n);

   if (j.contains("algorithm")) {
      if (j.at("algorithm").get<std::string>() == "metropolis")
         job.m_algorithm = Algorithm::Metropolis;
      else if (j.at("algorithm").get<std::string>() == "SW")
         job.m_algorithm = Algorithm::SW;
   }

   if (j.contains("image_mode")) {
      if (j.at("image_mode").get<std::string>() == "none")
         job.m_image_mode = None();
      else if (j.at("image_mode").get<std::string>() == "movie")
         job.m_image_mode = Movie();
      else if (j.at("image_mode").get<std::string>() == "intervals")
         job.m_image_mode = Intervals();
      else if (j.at("image_mode").get<std::string>() == "end")
         job.m_image_mode = EndImage();
   }

   if (j.contains("img_intervals")) {
      if (std::holds_alternative<Intervals>(job.m_image_mode))
         j.at("img_intervals").get_to(std::get<Intervals>(job.m_image_mode).m_between);
   }

   if (j.contains("image_path")) {
      if (std::holds_alternative<Movie>(job.m_image_mode))
         std::get<Movie>(job.m_image_mode).m_path = j.at("image_path").get<std::string>();
      else if (std::holds_alternative<Intervals>(job.m_image_mode))
         std::get<Intervals>(job.m_image_mode).m_path = j.at("image_path").get<std::string>();
      else if (std::holds_alternative<EndImage>(job.m_image_mode))
         std::get<Intervals>(job.m_image_mode).m_path = j.at("image_path").get<std::string>();
   }

   if (j.contains("physics_path"))
      job.m_physics_config.m_outputfile = j.at("physics_path").get<std::string>();
}

magneto::Job magneto::get_parsed_job(const std::string& file_contents){
   nlohmann::json json = nlohmann::json::parse(file_contents);
   Job job = json.get<Job>();
   return job;
}


bool magneto::operator==(const ImageSpinStart& a, const ImageSpinStart& b) {
   return a.m_path == b.m_path;
}
bool magneto::operator==(const RandomSpinStart& a, const RandomSpinStart& b) {
   return true;
}
bool magneto::operator==(const OneSpinStart& a, const OneSpinStart& b) {
   return true;
}
bool magneto::operator==(const None& a, const None& b) {
   return true;
}
bool magneto::operator==(const Movie& a, const Movie& b) {
   return std::tie(a.m_path, a.m_fps) == std::tie(b.m_path, b.m_fps);
}
bool magneto::operator==(const Intervals& a, const Intervals& b) {
   return std::tie(a.m_path, a.m_between) == std::tie(b.m_path, b.m_between);
}
bool magneto::operator==(const EndImage& a, const EndImage& b) {
   return a.m_path == b.m_path;
}
bool magneto::operator==(const PhysicsConfig& a, const PhysicsConfig& b) {
   return std::tie(a.m_outputfile, a.m_format, a.m_mode) == std::tie(b.m_outputfile, b.m_format, b.m_mode);
}

bool magneto::operator==(const ImageTempStart& a, const ImageTempStart& b){
   if (a.m_path != b.m_path)
      return false;
   if (!is_equal(a.t_min, b.t_min))
      return false;
   if (!is_equal(a.t_max, b.t_max))
      return false;
   return true;
}


bool magneto::operator==(const Job& a, const Job& b) {
   // use the std::tie trick for most
   if (std::tie(a.m_start_runs, a.m_L, a.m_n, a.m_spinstart, a.m_algorithm, a.m_image_mode, a.m_physics_config) !=
      std::tie(b.m_start_runs, b.m_L, b.m_n, b.m_spinstart, b.m_algorithm, b.m_image_mode, b.m_physics_config))
   {
      return false;
   }

   // temps extra because of the doubles
   if (a.m_temp.index() != b.m_temp.index())
      return false;
   if (std::holds_alternative<double>(a.m_temp)) {
      if (!is_equal(std::get<double>(a.m_temp), std::get<double>(a.m_temp)))
         return false;
   }
   else {
      if (!(std::get<ImageTempStart>(a.m_temp) == std::get<ImageTempStart>(b.m_temp)))
         return false;
   }
   return true;
}
