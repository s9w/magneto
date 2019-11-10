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
      //job1 = magneto::get_parsed_job(get_file_contents("test_configs/config1.json"));
      //job2 = magneto::get_parsed_job(get_file_contents("test_configs/config2.json"));
      //job3 = magneto::get_parsed_job(get_file_contents("test_configs/config3.json"));
      empty_job = magneto::get_parsed_job(get_file_contents("test_configs/config_empty.json"));
   }

   void SetUp() {}
   void TearDown() {}

   //magneto::Job job1;
   //magneto::Job job2;
   //magneto::Job job3;
   magneto::Job empty_job;
};


TEST_F(Jobs, EqualityOperators) {
   EXPECT_TRUE(magneto::PhysicsConfig({ "file", "{E}" }) == magneto::PhysicsConfig({ "file", "{E}" }));
   magneto::Job job1;
   magneto::Job job2;
   EXPECT_TRUE(job1==job2);
}

