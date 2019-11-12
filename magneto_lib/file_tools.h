#pragma once

#include <filesystem>
#include <optional>
#include "types.h"

namespace magneto {
	std::optional<std::string> get_file_contents(const std::filesystem::path& path);
	void write_string_to_file(const std::filesystem::path& path, const std::string& content);
   LatticeDType get_lattice_temps_from_png_file(const std::filesystem::path& path, const double temp_min, const double temp_max);
   LatticeType get_lattice_from_png_file(const std::filesystem::path& path);
}
