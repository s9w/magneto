#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/fmt.h>
#include <string>
#include <sstream>
#include <filesystem>

#include "VisualOutput.h"
#include "Job.h"
#include "logging.h"

namespace {
	/// <summary>[-1,1] -> [0,255]</summary>
	int get_255_value_from_pm_one(const int value) {
		return static_cast<unsigned char>((value + 1) / 2 * 255);
	}


   /// <summary>grid_buffer is already in [0,255] value range</summary>
	void write_png(const magneto::LatticeIType& grid_buffer, const std::filesystem::path& path) {
		std::vector<unsigned char> grid_png;
      const auto [Lx, Ly] = magneto::get_dimensions_of_lattice(grid_buffer);
		grid_png.reserve(Lx * Ly);
		for (const std::vector<int>& row : grid_buffer)
			for (const int elem : row)
				grid_png.emplace_back(static_cast<unsigned char>(elem));

      const int bpp = 1;
      const int pixel_row_stride = Lx * bpp;
		stbi_write_png(path.string().c_str(), Lx, Ly, 1, grid_png.data(), pixel_row_stride);
	}


	std::string get_rounded_string(const double number) {
		std::stringstream stream;
		stream << std::fixed << std::setprecision(3) << number;
		std::string temp_string = stream.str();
		return stream.str();
	}


	std::string get_png_directory_name(const std::string& temp_string) {
		return fmt::format("temp_png_{}", temp_string);
	}


	/// <summary>Transforms movie.mp4 into movie_2.266.mp4 to differentiate between movies of different temperatures</summary>
	std::filesystem::path get_movie_filename(const std::filesystem::path& base_name, const std::string temp_string) {
		auto new_name = base_name.stem();
		new_name += "_";
		new_name += temp_string;
		new_name += base_name.extension();
		return new_name;
	}


   /// <summary>Transforms image.png into image_2.266_{}.png</summary>
   std::string get_image_filename_pattern(const std::filesystem::path& base_name, const std::string& temp_string) {
      std::string new_name = base_name.stem().string();
      new_name += "_";
      new_name += temp_string;
      new_name += "_{}";
      new_name += base_name.extension().string();
      return new_name;
   }


	/// <summary>Fills the buffer with the grid data. Output is in [0,255] range.</summary>
	void add_grid_to_buffer(magneto::LatticeIType& buffer, const magneto::LatticeType& grid) {
      const auto [Lx, Ly] = magneto::get_dimensions_of_lattice(grid);
		for (unsigned int i = 0; i < Ly; ++i) {
			for (unsigned int j = 0; j < Lx; ++j) {
				buffer[i][j] += get_255_value_from_pm_one(grid[i][j]);
			}
		}
	}


	magneto::LatticeIType get_png_buffer_from_lattice(const magneto::LatticeType& grid) {
      const auto [Lx, Ly] = magneto::get_dimensions_of_lattice(grid);
		magneto::LatticeIType png_buffer(Ly, std::vector<int>(Lx));
		add_grid_to_buffer(png_buffer, grid);
		return png_buffer;
	}

} // namespace {}


magneto::MovieWriter::MovieWriter(const size_t Lx, const size_t Ly, const magneto::ImageMode& image_mode, const std::string& temp_string, const int blend_frames /*= 1*/)
	: m_framecount(0)
	, m_blendframes(blend_frames)
	, m_png_counter(0)
   , m_fps(image_mode.m_fps)
	, m_temp_directory_name(get_png_directory_name(temp_string))
	, m_output_filename(get_movie_filename(image_mode.m_path, temp_string))
	, m_buffer(Lx, Ly)
{
   clear_png_directory();
   std::filesystem::create_directory(m_temp_directory_name);
}


void magneto::MovieWriter::snapshot(const LatticeType& grid, const bool /*last_frame*/){
	m_buffer.add(grid);
	m_framecount++;

	if (m_framecount == m_blendframes) {
		const std::string filename = fmt::format("{}\\image_{}.png", m_temp_directory_name, m_png_counter);
		write_png(m_buffer.get_average(), filename);
		m_framecount = 0;
		m_png_counter++;
		m_buffer.clear();
	}
}

void magneto::MovieWriter::end_actions(){
	make_movie();
}


void magneto::MovieWriter::make_movie() const{
   get_logger()->info("Starting ffmpeg to write movie {}.", m_output_filename.string());
	const std::string cmd = fmt::format(
      "{} -y -hide_banner -loglevel panic -framerate {} -i {}\\image_%d.png -c:v libx264 {}", "ffmpeg.exe"
      , m_fps, m_temp_directory_name, m_output_filename.string()
   );
	system(cmd.c_str());
   get_logger()->info("Done writing movie {}.", m_output_filename.string());
	clear_png_directory();
}


void magneto::MovieWriter::clear_png_directory() const{
	std::filesystem::remove_all(m_temp_directory_name);
}


magneto::IntervalWriter::IntervalWriter(const size_t /*Lx*/, const size_t /*Ly*/, const ImageMode& image_mode, const std::string& temp_string)
   : m_framecount(0)
   , m_frame_intervals(image_mode.m_intervals)
   , m_fn_pattern(get_image_filename_pattern(image_mode.m_path, temp_string))
{}


void magneto::IntervalWriter::snapshot(const LatticeType& grid, const bool /*last_frame*/){
   ++m_framecount;

   if (m_framecount % m_frame_intervals == 0) {
		const std::string filename = fmt::format(m_fn_pattern, m_framecount);
      write_png(get_png_buffer_from_lattice(grid), filename);
   }
}


void magneto::IntervalWriter::end_actions(){}


magneto::TemporalAverageLattice::TemporalAverageLattice(const size_t Lx, const size_t Ly)
	: m_buffer(std::vector<std::vector<int>>(Ly, std::vector<int>(Lx, 0)))
	, m_recorded_frames(0)
{}


void magneto::TemporalAverageLattice::add(const LatticeType & grid){
	add_grid_to_buffer(m_buffer, grid);
	++m_recorded_frames;
}


magneto::LatticeIType magneto::TemporalAverageLattice::get_average(){
   const auto [Lx, Ly] = get_dimensions_of_lattice(m_buffer);
	for (unsigned int i = 0; i < Ly; ++i) {
		for (unsigned int j = 0; j < Lx; ++j) {
			m_buffer[i][j] /= m_recorded_frames;
		}
	}
	return m_buffer;
}


void magneto::TemporalAverageLattice::clear(){
   const auto [Lx, Ly] = get_dimensions_of_lattice(m_buffer);
	for (unsigned int i = 0; i < Ly; ++i) {
		for (unsigned int j = 0; j < Lx; ++j) {
			m_buffer[i][j] = 0;
		}
	}
   m_recorded_frames = 0;
}


magneto::EndImageWriter::EndImageWriter(const size_t /*Lx*/, const size_t /*Ly*/, const ImageMode& image_mode, const std::string& temp_string)
   : m_output_filename(get_movie_filename(image_mode.m_path, temp_string))
{}

void magneto::EndImageWriter::snapshot(const LatticeType& grid, const bool last_frame){
   if (!last_frame)
      return;
   write_png(get_png_buffer_from_lattice(grid), m_output_filename);
}

void magneto::EndImageWriter::end_actions()
{}



magneto::NullImageWriter::NullImageWriter(const size_t /*Lx*/, const size_t /*Ly*/, const ImageMode& /*image_mode*/, const std::string& /*temp_string*/)
{}

void magneto::NullImageWriter::snapshot(const LatticeType& /*grid*/, const bool /*last_frame*/)
{}

void magneto::NullImageWriter::end_actions()
{}
