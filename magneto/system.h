#pragma once

#include <filesystem>
#include <vector>

using GridType = std::vector<std::vector<int>>;

GridType get_randomized_system(const int L);
void write_png(const GridType& grid, const std::filesystem::path& path);