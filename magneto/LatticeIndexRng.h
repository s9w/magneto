#pragma once

#include <random>

namespace magneto {

	class LatticeIndexRng{
	public:
		LatticeIndexRng(const int L);
		int get();

	private:
		int m_max;
		std::mt19937_64 m_internal_rng;
		std::uniform_int_distribution<int> m_dist;
		int m_used;
		int m_current_number;
		int m_needed_bits;
	};

}

