#pragma once

#include <filesystem>
#include <vector>
#include <random>

using GridType = std::vector<std::vector<int>>;

GridType get_randomized_system(const int L);
void write_png(const GridType& grid, const std::filesystem::path& path);
int get_dE(const GridType& grid, int i, int j);
void metropolis_sweeps(GridType& grid, std::mt19937& rng, const std::vector<double>& exp_values, const int n = 1);