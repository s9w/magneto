#pragma once

#include "export_macro.h"

#include <nlohmann/json.hpp>

#include <filesystem>
#include <variant>
#include <optional>
#include "types.h"


namespace magneto {
   enum class Algorithm { Metropolis, SW };
   enum class SpinStartMode { Random, Image };
   enum class TempStartMode { Single, Many, Image };

   enum class ImageOrMovie { None, Endimage, Intervals, Movie };
   struct ImageMode {
      ImageOrMovie m_mode = ImageOrMovie::Endimage;
      unsigned int m_intervals = 10;
      unsigned int m_fps = 30;
      std::filesystem::path m_path = "magneto_images";
   };

   struct PhysicsConfig {
      std::filesystem::path m_outputfile = "magneto_results.txt";
      std::string m_format = "T: {T:<5.3f},\tEnergy: {E:<5.3f},\tcv: {cv:<5.3f}, mag: {M:<5.3f}, chi: {chi:<5.3f}";
   };
   

   struct JsonJob {
      SpinStartMode m_spin_start_mode = SpinStartMode::Random;
      std::filesystem::path m_spin_start_image_path;
      std::filesystem::path m_temperature_image;
      TempStartMode m_temp_mode = TempStartMode::Single;

      double m_t_single = 2.26;

      // this for image temp start and many temps
      double m_t_min = 2.0;
      double m_t_max = 2.5;

      std::vector<double> m_temperatures;

      // this only for many temps
      unsigned int m_temp_steps = 3;

      // How many Swendsen-Wang runs before anything is being recorded/computed
      unsigned int m_start_runs = 0;

      unsigned int m_L = 400;
      unsigned int m_n = 100;

      // Algorithm used for propagation (after the initial start runs)
      Algorithm m_algorithm = Algorithm::Metropolis;

      ImageMode m_image_mode;

      PhysicsConfig m_physics_config;
   };


   struct Job {
      // system start
      unsigned int m_L = 400;
      //std::variant<LatticeDType, std::vector<double>> T;
      LatticeType initial_spins;
      unsigned int m_start_runs = 0;

      // system evolution
      Algorithm m_algorithm = Algorithm::Metropolis;
      unsigned int m_n = 100;

      // output
      ImageMode m_image_mode;
      PhysicsConfig m_physics_config;
   };

   CLASS_DECLSPEC bool operator==(const ImageMode& a, const ImageMode& b);
   CLASS_DECLSPEC bool operator==(const PhysicsConfig& a, const PhysicsConfig& b);
   CLASS_DECLSPEC bool operator==(const JsonJob& a, const JsonJob& b);

   void from_json(const nlohmann::json& j, magneto::JsonJob& job);
   CLASS_DECLSPEC JsonJob get_parsed_job(const std::string& file_contents);
   CLASS_DECLSPEC std::tuple<Job, std::variant<LatticeDType, std::vector<double>>> get_job(const JsonJob& json_job);
   CLASS_DECLSPEC std::optional<JsonJob> get_parsed_job(const std::filesystem::path& path);

}
