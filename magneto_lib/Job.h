#pragma once

#include "export_macro.h"

#include "nlohmann/json.hpp"

#include <filesystem>
#include <variant>


namespace magneto {
   enum class Algorithm { Metropolis, SW };

   struct RandomSpinStart {};
   struct OneSpinStart {};
   struct ImageSpinStart { std::filesystem::path m_path; };   

   struct ImageTempStart {
      std::filesystem::path m_path;
      double t_min = 1.0;
      double t_max = 5.0;
   };

   struct None {};

   // movie from every iteration of the system
   struct Movie {
      std::filesystem::path m_path = "magneto_movie.mp4";
      unsigned int m_fps = 30;
   };

   // png files every n iterations of the system
   struct Intervals {
      std::filesystem::path m_path = "magneto_snapshots.png";
      unsigned int m_between = 10;
   };

   // Only one png file from the last system state
   struct EndImage { std::filesystem::path m_path; };

   enum class PhysicsMode { EveryStep, AtTheEnd };

   struct PhysicsConfig {
      std::filesystem::path m_outputfile = "magneto_results.txt";
      std::string m_format = "E: {E}";
      PhysicsMode m_mode = PhysicsMode::AtTheEnd;
   };
   

   struct Job {
      std::variant<RandomSpinStart, OneSpinStart, ImageSpinStart> m_spinstart = RandomSpinStart();
      std::variant<double, ImageTempStart> m_temp = 2.26;

      // How many Swendsen-Wang runs before anything is being recorded/computed
      unsigned int m_start_runs = 0;

      unsigned int m_L = 400;
      unsigned int m_n = 100;

      // Algorithm used for propagation (after the initial start runs)
      Algorithm m_algorithm = Algorithm::Metropolis;

      std::variant<None, Movie, Intervals, EndImage> m_image_mode = EndImage{"magneto_image.png"};

      PhysicsConfig m_physics_config;
   };

   CLASS_DECLSPEC bool operator==(const ImageSpinStart& a, const ImageSpinStart& b);
   CLASS_DECLSPEC bool operator==(const RandomSpinStart& a, const RandomSpinStart& b);
   CLASS_DECLSPEC bool operator==(const OneSpinStart& a, const OneSpinStart& b);
   CLASS_DECLSPEC bool operator==(const None& a, const None& b);
   CLASS_DECLSPEC bool operator==(const Movie& a, const Movie& b);
   CLASS_DECLSPEC bool operator==(const Intervals& a, const Intervals& b);
   CLASS_DECLSPEC bool operator==(const EndImage& a, const EndImage& b);
   CLASS_DECLSPEC bool operator==(const PhysicsConfig& a, const PhysicsConfig& b);
   CLASS_DECLSPEC bool operator==(const ImageTempStart& a, const ImageTempStart& b);
   CLASS_DECLSPEC bool operator==(const Job& a, const Job& b);

   void from_json(const nlohmann::json& j, magneto::Job& job);
   CLASS_DECLSPEC Job get_parsed_job(const std::string& file_contents);

}
