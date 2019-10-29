#pragma once

class ProgressIndicator{
public:
	void set_progress(const int new_progress_percent);
	void set_progress(const int new_progress, const int total);
	void write_progress() const;

private:
	int m_progress_percent = 0;
};

