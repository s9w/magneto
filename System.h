#include <vector>
#include <random>
#include "Config.h"
#include "helpers.h"

enum resultTypes { energy, mag, cv, chi, states, corrfun };

struct correlationPoint {
	unsigned int i, j;
	std::vector<double> corr_ab, corr_a;
};

class System{
public:
	System(Config, LabConfig&);
	void compute();
	Config cfg;
	std::vector<std::string> results;
private:
	void measure();
	std::vector<std::vector<int> > genRandomSystem(int seedOffset);
	std::vector<std::vector<int> > getRelaxedSys(int seedOffsets);
	std::vector<std::vector<int> > getFileState(std::string filename);
	void metropolis_sweeps(unsigned int);
	void wangRuns(const unsigned int n);
	std::mt19937 gen_metro;
	std::vector<double> exp_values;

	bool calc_e=false, calc_m=false, calc_cv=false, calc_chi=false, calc_states=false, calc_corrfun=false;
	std::vector<std::vector<int> > grid;
	double e_avg=0.0, e2_avg=0.0, m_avg=0.0, m2_avg=0.0;
	void recordResults();

	std::vector<correlationPoint> correlations;
	unsigned int corr_range, corr_count;
};
