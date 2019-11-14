#pragma once

#include <filesystem>
#include <optional>
#include "types.h"
#include <functional>

namespace magneto {
	std::optional<std::string> get_file_contents(const std::filesystem::path& path);
	void write_string_to_file(const std::filesystem::path& path, const std::string& content);
   std::optional<LatticeDType> get_lattice_temps_from_png_file(const std::filesystem::path& path, const double temp_min, const double temp_max);
   
   std::optional<LatticeType> get_spin_state_from_png(const std::filesystem::path& path);

   class FileResizer {
   public:
      FileResizer(const std::filesystem::path& path, const int new_x, const int new_y);
      ~FileResizer();
      std::filesystem::path get_temp_file() const;
   private:
      std::filesystem::path m_temp_path;

      static std::filesystem::path get_resized_image_path(const std::filesystem::path& original_path);
   };


   // Resize the image at the path to the specified size with some function
   template<class T>
   T get_resized_data(
      const std::filesystem::path& path,
      const int Lx,
      const int Ly,
      const std::function<T(const std::filesystem::path& path)>& f
   ) {
      FileResizer resizer(path, Lx, Ly);
      return f(resizer.get_temp_file());
   }

}
