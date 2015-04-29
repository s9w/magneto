#include <vector>
#include "Config.h"

class System{
public:
	System(Config);
	void compute();
	void measure();
	Config cfg;
	std::vector<std::string> results;
private:
	bool calc_e, calc_m, calc_cv, calc_chi;
	std::vector<std::vector<int> > grid;
	double e_avg=0.0, e2_avg=0.0, m_avg=0.0, m2_avg=0.0;
};
