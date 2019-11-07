#include "pch.h"
#include <fstream>
#include "../magneto_lib/Job.h"

namespace {
   std::string get_file_contents(const std::filesystem::path& path) {
      const std::filesystem::path current = std::filesystem::current_path(); // magneto/x64/Release/
      const std::filesystem::path magneto_main = current.parent_path().parent_path(); // magneto/
      std::filesystem::path complete_path = magneto_main / "Tests" / path;
      std::ifstream filestream(complete_path);
      if (!filestream.is_open())
         return "";
      std::stringstream buffer;
      buffer << filestream.rdbuf();
      filestream.close();
      return buffer.str();
   }
}


class Jobs : public ::testing::Test {
public:
   Jobs() {
      job1 = magneto::get_parsed_job(get_file_contents("test_configs/config1.json"));
      job2 = magneto::get_parsed_job(get_file_contents("test_configs/config2.json"));
      job3 = magneto::get_parsed_job(get_file_contents("test_configs/config3.json"));
      empty_job = magneto::get_parsed_job(get_file_contents("test_configs/config_empty.json"));
   }

   void SetUp() {}
   void TearDown() {}

   magneto::Job job1;
   magneto::Job job2;
   magneto::Job job3;
   magneto::Job empty_job;
};


TEST_F(Jobs, EqualityOperators) {
   EXPECT_TRUE(magneto::RandomSpinStart() == magneto::RandomSpinStart());
   EXPECT_TRUE(magneto::OneSpinStart() == magneto::OneSpinStart());

   EXPECT_TRUE(magneto::ImageSpinStart{ "a" } == magneto::ImageSpinStart{"a"});
   EXPECT_FALSE(magneto::ImageSpinStart{ "a" } == magneto::ImageSpinStart{"b"});

   EXPECT_TRUE(magneto::ImageTempStart({ "a", 1.0, 1.0 }) == magneto::ImageTempStart({ "a", 1.0, 1.0 }));
   EXPECT_FALSE(magneto::ImageTempStart({ "a", 1.0, 1.0 }) == magneto::ImageTempStart({ "b", 1.0, 1.0 }));
   EXPECT_FALSE(magneto::ImageTempStart({ "a", 1.0, 1.0 }) == magneto::ImageTempStart({ "a", 1.1, 1.0 }));
   EXPECT_FALSE(magneto::ImageTempStart({ "a", 1.0, 1.0 }) == magneto::ImageTempStart({ "a", 1.0, 1.1 }));
   EXPECT_TRUE (magneto::ImageTempStart({ "a", 1.0, 1.0 }) == magneto::ImageTempStart({ "a", 1.0, 1.0001 }));

   EXPECT_TRUE(magneto::None() == magneto::None());

   EXPECT_TRUE (magneto::Movie({ "a", 1 }) == magneto::Movie({ "a", 1 }));
   EXPECT_FALSE(magneto::Movie({ "a", 1 }) == magneto::Movie({ "b", 1 }));
   EXPECT_FALSE(magneto::Movie({ "a", 1 }) == magneto::Movie({ "a", 2 }));

   EXPECT_TRUE(magneto::Intervals({ "a", 1 }) == magneto::Intervals({ "a", 1 }));
   EXPECT_FALSE(magneto::Intervals({ "a", 1 }) == magneto::Intervals({ "b", 1 }));
   EXPECT_FALSE(magneto::Intervals({ "a", 1 }) == magneto::Intervals({ "a", 2 }));

   EXPECT_TRUE(magneto::EndImage({ "a" }) == magneto::EndImage({ "a" }));
   EXPECT_FALSE(magneto::EndImage({ "a" }) == magneto::EndImage({ "b" }));

   EXPECT_TRUE(magneto::PhysicsConfig({ "a", "a", magneto::PhysicsMode::AtTheEnd }) == magneto::PhysicsConfig({ "a", "a", magneto::PhysicsMode::AtTheEnd }));
   EXPECT_FALSE(magneto::PhysicsConfig({ "a", "a", magneto::PhysicsMode::AtTheEnd }) == magneto::PhysicsConfig({ "b", "a", magneto::PhysicsMode::AtTheEnd }));
   EXPECT_FALSE(magneto::PhysicsConfig({ "a", "a", magneto::PhysicsMode::AtTheEnd }) == magneto::PhysicsConfig({ "a", "b", magneto::PhysicsMode::AtTheEnd }));
}


TEST_F(Jobs, compare_job_1) {
   magneto::Job target;
   target.m_spinstart = magneto::RandomSpinStart();
   target.m_temp = 2.0;
   target.m_start_runs = 50;
   target.m_L = 500;
   target.m_n = 300;
   target.m_algorithm = magneto::Algorithm::Metropolis;
   target.m_image_mode = magneto::Movie{ "movie.mp4", 30 };
   target.m_physics_config.m_outputfile = "data.txt";
   target.m_physics_config.m_format = "E: {E}";
   target.m_physics_config.m_mode = magneto::PhysicsMode::AtTheEnd;
   EXPECT_TRUE(job1 == target);

   target.m_temp = 2.1;
   EXPECT_FALSE(job1 == target);
}

TEST_F(Jobs, compare_job_2) {
   magneto::Job target;
   target.m_spinstart = magneto::OneSpinStart();
   target.m_temp = magneto::ImageTempStart{"temp_input.png", 1.0, 6.0};
   target.m_start_runs = 50;
   target.m_L = 500;
   target.m_n = 300;
   target.m_algorithm = magneto::Algorithm::Metropolis;
   target.m_image_mode = magneto::Intervals{ "img.png", 10 };
   target.m_physics_config.m_outputfile = "data.txt";
   target.m_physics_config.m_format = "E: {E}";
   target.m_physics_config.m_mode = magneto::PhysicsMode::AtTheEnd;
   EXPECT_TRUE(job2 == target);
}

