#include <vector>
#include <random>
#include "Config.h"
#include "helpers.h"

enum resultTypes { energy, mag, cv, chi, states, corrfun };

struct CorrelationPoint {
	unsigned int i, j;
	std::vector<double> corr_ab, corr_a;
};

typedef std::vector<std::vector<int> > grid_type;

class System{
public:
	System(Config&, LabConfig&);
	void compute();
	const Config cfg;
	std::vector<std::string> results;
private:
	void measure();
	grid_type genRandomSystem(int seedOffset);
	grid_type getRelaxedSys(int seedOffsets);
	grid_type getFileState(std::string filename);
	grid_type& metropolis_sweeps(grid_type&, unsigned int);
	grid_type& wangRuns(grid_type&, const unsigned int n);
	void recordResults();

	std::mt19937 rng;
	std::vector<double> exp_values;
	bool calc_e=false, calc_m=false, calc_cv=false, calc_chi=false, calc_states=false, calc_corrfun=false;
	grid_type grid;
	double e_avg=0.0, e2_avg=0.0, m_avg=0.0, m2_avg=0.0;
	std::vector<CorrelationPoint> correlations;
	unsigned int corr_range, corr_count;
	std::vector<int> access_i;
	std::vector<int> access_j;
};
