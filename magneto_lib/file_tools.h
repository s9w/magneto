#pragma once

#include <filesystem>
#include <optional>
#include "types.h"
#include <functional>

namespace magneto {
   /// <summary>abc.png -> abc_resized.png</summary>
   std::filesystem::path get_resized_image_path(const std::filesystem::path& original_path);

	std::optional<std::string> get_file_contents(const std::filesystem::path& path);
	void write_string_to_file(const std::filesystem::path& path, const std::string& content);
   LatticeDType get_lattice_temps_from_png_file(const std::filesystem::path& path, const double temp_min, const double temp_max);
   
   std::optional<LatticeType> get_spin_state_from_png(const std::filesystem::path& path);


   // Resize the image at the path to the specified size with some function
   template<class T>
   T get_resized_data(
      const std::filesystem::path& path,
      const int Lx,
      const int Ly,
      const std::function<T(const std::filesystem::path& path)>& f
   ) {
      std::filesystem::path resized_path = get_resized_image_path(path);

      const std::string cmd = fmt::format(
         "{ffmpeg} -y -hide_banner -loglevel panic -i {input} -vf scale={x}:{y} {output}",
         fmt::arg("ffmpeg", "ffmpeg.exe"),
         fmt::arg("input", path.string()),
         fmt::arg("output", resized_path.string()),
         fmt::arg("x", Lx),
         fmt::arg("y", Ly)
      );
      system(cmd.c_str());

      const T resized_state = f(resized_path);
      if (!std::filesystem::remove(resized_path)) {
         magneto::get_logger()->error("Couldn't remove temporary file {}.", resized_path.string());
      }
      return resized_state;
   }

}
