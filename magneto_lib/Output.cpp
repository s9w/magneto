#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/fmt.h>
#include <string>
#include <sstream>
#include <filesystem>

#include "Output.h"
#include "Job.h"

namespace {
	/// <summary>[-1,1] -> [0,255]</summary>
	int get_255_value_from_pm_one(const int value) {
		return static_cast<unsigned char>((value + 1) / 2 * 255);
	}


	void write_png(const magneto::LatticeType& grid, const std::filesystem::path& path) {
		std::vector<unsigned char> grid_png;
		const int L = static_cast<int>(grid.size());
		grid_png.reserve(L * L);
		for (const auto& row : grid)
			for (const int elem : row)
				grid_png.emplace_back(static_cast<unsigned char>(elem));

		stbi_write_png(path.string().c_str(), L, L, 1, grid_png.data(), L * 1);
	}


	std::string get_rounded_string(const double number) {
		std::stringstream stream;
		stream << std::fixed << std::setprecision(3) << number;
		std::string temp_string = stream.str();
		return stream.str();
	}


	std::string get_png_directory_name(const double T) {
		return fmt::format("temp_png_{}", get_rounded_string(T));
	}


	/// <summary>Transforms movie.mp4 into movie_2.266.mp4 to differentiate between movies of different temperatures</summary>
	std::filesystem::path get_movie_filename(const std::filesystem::path& base_name, const double T) {
		auto new_name = base_name.stem();
		new_name += "_";
		new_name += get_rounded_string(T);
		new_name += base_name.extension();
		return new_name;
	}
}


magneto::MovieWriter::MovieWriter(const size_t L, const magneto::ImageMode& image_mode, const double T, const int blend_frames /*= 1*/)
	: m_gridbuffer(std::vector<std::vector<int>>(L, std::vector<int>(L, 0)))
	, m_framecount(0)
	, m_blendframes(blend_frames)
	, m_png_counter(0)
   , m_mode(image_mode)
	, m_temp_directory_name(get_png_directory_name(T))
	, m_output_filename(get_movie_filename(m_mode.m_path, T))
{
   if (m_mode.m_mode == ImageOrMovie::None)
      return;
	
   if (m_mode.m_mode == ImageOrMovie::Movie) {
      clear_png_directory();
      std::filesystem::create_directory(m_temp_directory_name);
   }
}


void magneto::MovieWriter::snapshot(const LatticeType& grid, const bool last_frame){
   if (m_mode.m_mode == ImageOrMovie::None)
      return;

	const int L = static_cast<int>(grid.size());
	for (int i = 0; i < L; ++i) {
		for (int j = 0; j < L; ++j) {
			m_gridbuffer[i][j] += get_255_value_from_pm_one(grid[i][j]);
		}
	}
	m_framecount++;

	if (m_framecount == m_blendframes) {
		for (int i = 0; i < L; ++i) {
			for (int j = 0; j < L; ++j) {
				m_gridbuffer[i][j] /= m_blendframes;
			}
		}
		const std::string filename = fmt::format("{}\\image_{}.png", m_temp_directory_name, m_png_counter);
		write_png(m_gridbuffer, filename);
		m_framecount = 0;
		m_png_counter++;
		clear_buffer();
	}
}

void magneto::MovieWriter::end_actions(){
	make_movie();
}


void magneto::MovieWriter::make_movie() const{
   if (m_mode.m_mode == ImageOrMovie::None)
      return;
	const std::string cmd = fmt::format(
      "{} -y -hide_banner -loglevel panic -framerate {} -i {}\\image_%d.png -c:v libx264 {}", "ffmpeg.exe"
      , m_mode.m_fps, m_temp_directory_name, m_output_filename.string()
   );
	system(cmd.c_str());
	clear_png_directory();
}


void magneto::MovieWriter::clear_buffer(){
	const int L = static_cast<int>(m_gridbuffer.size());
	for (int i = 0; i < L; ++i) {
		for (int j = 0; j < L; ++j) {
			m_gridbuffer[i][j] = 0;
		}
	}
}


void magneto::MovieWriter::clear_png_directory() const{
	std::filesystem::remove_all(m_temp_directory_name);
}


magneto::IntervalWriter::IntervalWriter(const size_t L, const ImageMode& image_mode, const double T){
}

void magneto::IntervalWriter::snapshot(const LatticeType& grid, const bool last_frame){
}

void magneto::IntervalWriter::end_actions(){
}
