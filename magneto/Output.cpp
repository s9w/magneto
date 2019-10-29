#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <string>
#include <filesystem>

#include "Output.h"

namespace {
	void write_png(const magneto::LatticeType& grid, const std::filesystem::path& path) {
		std::vector<unsigned char> grid_png;
		const int L = static_cast<int>(grid.size());
		grid_png.reserve(L * L);
		for (const auto& row : grid)
			for (const int elem : row)
				grid_png.emplace_back(elem);

		stbi_write_png(path.string().c_str(), L, L, 1, grid_png.data(), L * 1);
	}
}


magneto::Output::Output(const int L, const int blend_frames /*= 1*/)
	: m_gridbuffer(std::vector<std::vector<int>>(L, std::vector<int>(L, 0)))
	, m_framecount(0)
	, m_blendframes(blend_frames)
	, m_png_counter(0)
{
	clear_png_directory();
	std::filesystem::create_directory("png");
}


void magneto::Output::photograph(const LatticeType& grid){
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
		const std::string filename = std::string("png\\test_") + std::to_string(m_png_counter) + ".png";
		write_png(m_gridbuffer, filename);
		m_framecount = 0;
		m_png_counter++;
		clear_buffer();
	}
}


void magneto::Output::make_movie() const{
	const std::string ffmpeg_path = "C:\\Users\\swerhausen\\Dropbox\\magneto\\magneto\\ffmpeg.exe";
	const std::string cmd = ffmpeg_path + " -y -i png\\test_%d.png -c:v huffyuv test.mkv";
	system(cmd.c_str());
	//clear_png_directory();
}


void magneto::Output::clear_buffer(){
	const int L = static_cast<int>(m_gridbuffer.size());
	for (int i = 0; i < L; ++i) {
		for (int j = 0; j < L; ++j) {
			m_gridbuffer[i][j] = 0;
		}
	}
}


void magneto::Output::clear_png_directory() const{
	std::filesystem::remove_all("png");
}

unsigned char magneto::get_255_char_from_pm_one(const int value){
	return (value + 1) / 2 * 255;
}


unsigned short magneto::get_255_value_from_pm_one(const int value){
	return (value + 1) / 2 * 255;
}
