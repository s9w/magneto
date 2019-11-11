#include "file_tools.h"
#include "logging.h"
#include <fstream>
#include <sstream>


std::optional<std::string> magneto::get_file_contents(const std::filesystem::path& path){
   const std::filesystem::path current = std::filesystem::current_path();
   std::filesystem::path complete_path = current / path;
   std::ifstream filestream(complete_path);
   if (!filestream.is_open())
      return std::nullopt;
   std::stringstream buffer;
   buffer << filestream.rdbuf();
   filestream.close();
   return buffer.str();
}


void magneto::write_string_to_file(const std::filesystem::path& path, const std::string& content){
   std::ofstream file(path);
   if (!file.is_open()) {
      magneto::get_logger()->error("Couldn't open file {} for writing result.", path.string());
      return;
   }
   file << content;
   file.close();
}
