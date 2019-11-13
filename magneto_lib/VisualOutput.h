#pragma once

#include "IsingSystem.h"
#include "Job.h"

namespace magneto {

   class VisualOutput {
   public:
      virtual void snapshot(const LatticeType& grid, const bool last_frame = false) = 0;
      virtual void end_actions() = 0;
   };


   class TemporalAverageLattice {
   public:
      TemporalAverageLattice(const size_t Lx, const size_t Ly);

      /// <summary>Expects a +-1 grid input</summary>
      void add(const LatticeType& grid);

      /// <summary>Returns the temporal average over the recorded data in [0,255] range</summary>
      LatticeIType get_average();
      void clear();

   private:
      LatticeIType m_buffer;
      int m_recorded_frames;
   };


   class MovieWriter : public VisualOutput {
   public:
      MovieWriter(const size_t Lx, const size_t Ly, const ImageMode& image_mode, const std::string& temp_string, const int blend_frames = 1);
      void snapshot(const LatticeType& grid, const bool last_frame = false);
      void end_actions();
      void make_movie() const;

   private:
      void clear_png_directory() const;

      int m_framecount;
      int m_blendframes;
      int m_png_counter;
      unsigned int m_fps;
      std::string m_temp_directory_name;
      std::filesystem::path m_output_filename;
      TemporalAverageLattice m_buffer;
   };


   class IntervalWriter : public VisualOutput {
   public:
      IntervalWriter(const size_t Lx, const size_t Ly, const ImageMode& image_mode, const std::string& temp_string);
      void snapshot(const LatticeType& grid, const bool last_frame = false);
      void end_actions();

   private:
      int m_framecount;
      int m_frame_intervals;
      std::string m_fn_pattern;
   };


   class EndImageWriter : public VisualOutput {
   public:
      EndImageWriter(const size_t Lx, const size_t Ly, const ImageMode& image_mode, const std::string& temp_string);
      void snapshot(const LatticeType& grid, const bool last_frame = false);
      void end_actions();

   private:
      std::filesystem::path m_output_filename;
   };


   class NullImageWriter : public VisualOutput {
   public:
      NullImageWriter(const size_t Lx, const size_t Ly, const ImageMode& image_mode, const std::string& temp_string);
      void snapshot(const LatticeType& grid, const bool last_frame = false);
      void end_actions();

   private:
   };

}

