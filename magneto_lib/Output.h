#pragma once

#include "system.h"

namespace magneto {

	class Output
	{
	public:
		Output(const size_t L, const int blend_frames = 1);
		void photograph(const LatticeType& grid);
		void make_movie() const;
	private:
		void clear_buffer();
		void clear_png_directory() const;

		LatticeType m_gridbuffer;
		int m_framecount;
		int m_blendframes;
		int m_png_counter;
	};

}

