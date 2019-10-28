#include <chrono>

#include "LatticeIndexRng.h"

namespace {
	int how_many_bits_needed(int number) {
		return static_cast<int>(std::ceil(std::log2(number + 1)));
	}
}

magneto::LatticeIndexRng::LatticeIndexRng(const int L)
	: m_internal_rng(std::chrono::system_clock::now().time_since_epoch().count())
	, m_max(L*L)
	, m_dist(0, L*L)
	, m_used(0)
	, m_current_number(0)
{
	m_current_number = m_internal_rng();
	m_needed_bits = how_many_bits_needed(L*L);
}

int magneto::LatticeIndexRng::get(){
	if (m_used >= 64 - m_needed_bits) {
		m_current_number = m_internal_rng();
		m_used = 0;
	}
	const int mask = (1 << m_needed_bits) - 1;
	int output_number = m_current_number & mask;
	m_current_number = m_current_number >> 2;
	m_used += 2;
	if (output_number > m_max - 1)
		return get();
	return output_number;
}


