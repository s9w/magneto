#pragma once

#include "IsingSystem.h"
#include "Job.h"

namespace magneto {

	class VisualOutputInterface {
		virtual void snapshot(const LatticeType& grid, const bool last_frame = false) = 0;
		virtual void end_actions() = 0;
	};


	class MovieWriter : public VisualOutputInterface{
	public:
		MovieWriter(const size_t L, const ImageMode& image_mode, const double T, const int blend_frames = 1);
		void snapshot(const LatticeType& grid, const bool last_frame = false);
		void end_actions();
		void make_movie() const;

	private:
		void clear_buffer();
		void clear_png_directory() const;

		LatticeType m_gridbuffer;
		int m_framecount;
		int m_blendframes;
		int m_png_counter;
      ImageMode m_mode;
		std::string m_temp_directory_name;
		std::filesystem::path m_output_filename;
	};


	class IntervalWriter : public VisualOutputInterface {
	public:
		IntervalWriter(const size_t L, const ImageMode& image_mode, const double T);
		void snapshot(const LatticeType& grid, const bool last_frame = false);
		void end_actions();
	};

}

