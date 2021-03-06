#include "file_tools.h"

#include "logging.h"

#include <fstream>
#include <sstream>

#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>


namespace {
   /// <summary>Converts [0,255] to [-1,1]. </summary>
   template<class TIn, class TOut>
   TOut get_pm1_from_255_value(const TIn value) {
      return static_cast<TOut>(value / 255 * 2 - 1);
   }


   // Making sure this will not be used in an unsigned char
   template<>
   unsigned char get_pm1_from_255_value(const unsigned char value) = delete;


   magneto::LatticeType get_lattice_from_monochrome_bitmap_data(
      unsigned char* png_data,
      const int Lx,
      const int Ly
   ) {
      magneto::LatticeType lattice(Ly, std::vector<char>(Lx));
      for (int i = 0; i < Ly; ++i) {
         for (int j = 0; j < Lx; ++j) {
            lattice[i][j] = get_pm1_from_255_value<unsigned char, char>(png_data[i * Lx + j]);
         }
      }
      return lattice;
   }


   magneto::LatticeType get_lattice_from_rgba_bitmap_data(
      unsigned char* png_data,
      const int Lx,
      const int Ly
   ) {
      magneto::LatticeType lattice(Ly, std::vector<char>(Lx));
      for (int i = 0; i < Ly; ++i) {
         for (int j = 0; j < Lx; ++j) {
            unsigned int value = 0;
            for (int k = 0; k < 4; ++k) {
               value += png_data[(i * Ly + j) * 4 + k];
            }
            value /= 4;
            lattice[i][j] = get_pm1_from_255_value<unsigned int, char>(value);
         }
      }
      return lattice;
   }


   magneto::LatticeType get_lattice_from_png_data(
      unsigned char* png_data, 
      const int bpp,
      const int Lx,
      const int Ly
   ) {
      if (bpp == 1)
         return get_lattice_from_monochrome_bitmap_data(png_data, Lx, Ly);
      else if (bpp == 4)
         return get_lattice_from_rgba_bitmap_data(png_data, Lx, Ly);
   }

} // namespace {}


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


std::optional<magneto::LatticeDType> magneto::get_lattice_temps_from_png_file(
   const std::filesystem::path& path, const double temp_min, const double temp_max
) {
   int x, y, bpp;
   unsigned char* image_data = stbi_load(path.string().c_str(), &x, &y, &bpp, 0);
   if (image_data == nullptr) {
      magneto::get_logger()->error("Couldn't open image {}", path.string());
      return std::nullopt;
   }
   magneto::LatticeDType temps(y, std::vector<double>(x, 2.26));
   const double temp_factor = temp_max - temp_min;
   for (int i = 0; i < y; ++i) {
      for (int j = 0; j < x; ++j) {
         int value = 0;
         for (int k = 0; k < bpp; ++k) {
            value += image_data[(i * x + j) * bpp + k];
         }
         value /= bpp;
         temps[i][j] = temp_min + value * 1.0 / 256 * temp_factor;
      }
   }
   stbi_image_free(image_data);
   return temps;
}


std::optional<magneto::LatticeType> magneto::get_spin_state_from_png(const std::filesystem::path& path) {
   int x, y, bpp;
   unsigned char* data = stbi_load(path.string().c_str(), &x, &y, &bpp, 0);
   if (data == nullptr) {
      magneto::get_logger()->error("Couldn't open image {}", path.string());
      return std::nullopt;
   }
   magneto::LatticeType lattice = get_lattice_from_png_data(data, bpp, x, y);
   stbi_image_free(data);
   return lattice;
}


magneto::FileResizer::FileResizer(const std::filesystem::path& path, const int new_x, const int new_y)
   : m_temp_path(get_resized_image_path(path))
{
   const std::string cmd = fmt::format(
      "{ffmpeg} -y -hide_banner -loglevel panic -i {input} -vf scale={x}:{y} {output}",
      fmt::arg("ffmpeg", "ffmpeg.exe"),
      fmt::arg("input", path.string()),
      fmt::arg("output", m_temp_path.string()),
      fmt::arg("x", new_x),
      fmt::arg("y", new_y)
   );
   system(cmd.c_str());
}


magneto::FileResizer::~FileResizer(){
   if (!std::filesystem::remove(m_temp_path)) {
      magneto::get_logger()->error("Couldn't remove temporary file {}.", m_temp_path.string());
   }
}


std::filesystem::path magneto::FileResizer::get_temp_file() const{
   return m_temp_path;
}

std::filesystem::path magneto::FileResizer::get_resized_image_path(const std::filesystem::path& original_path){
   std::filesystem::path new_path = original_path;
   std::filesystem::path new_filename = original_path.stem();
   new_filename += "_resized";
   new_filename += original_path.extension();

   new_path.replace_filename(new_filename);
   return new_path;
}
