#include <iostream>
#include "ProgressIndicator.h"


void ProgressIndicator::set_progress(const int new_progress_percent){
	m_progress_percent = new_progress_percent;
}


void ProgressIndicator::set_progress(const int new_progress, const int total){
	set_progress(static_cast<int>(std::round(100.0 * new_progress / total)));
}


void ProgressIndicator::write_progress() const{
	std::cout << "                    " << std::flush;
	std::cout << '\r' << std::flush;
	std::cout << m_progress_percent << " %" << std::flush;
}
