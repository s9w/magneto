#include <map>
#include <vector>
#include "Config.h"

class System{
public:
	System(Config, LabConfig&);
	void compute();
	void measure();
	Config cfg;
	std::vector<std::string> results = {"", "", "", "", ""};
	std::vector<std::string> resultsStates;
private:
	bool calc_e=false, calc_m=false, calc_cv=false, calc_chi=false, calc_states=false, calc_corr=false;
	std::vector<std::vector<int> > grid;
	double e_avg=0.0, e2_avg=0.0, m_avg=0.0, m2_avg=0.0;
	std::vector<double> corr_ab, corr_a, corr_b;
	void recordResults();
};