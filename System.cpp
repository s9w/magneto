#include "System.h"
#include "algs.h"
#include "physics.h"

std::vector<std::vector<int> > getRelaxedSys(const unsigned int L, const double T, double J, unsigned int n1, int alg, int seedOffset=0);

template<typename T>
std::string to_string(T const & value) {
    std::stringstream ss;
    ss << value;
    return ss.str();
}

System::System(Config p_cfg) {
    cfg = p_cfg;
    if(!cfg.fileEnergy.empty())
        calc_e = true;
    if(!cfg.fileMag.empty())
        calc_m = true;
    if(!cfg.fileCv.empty())
        calc_cv = true;
    if(!cfg.fileChi.empty())
        calc_chi = true;
    grid = getRelaxedSys(cfg.L, cfg.T, cfg.J, cfg.n1, cfg.alg1);
}

void System::compute() {
    for(int i=0; i<cfg.n2; ++i) {
        measure();

        //  evolve
        if(cfg.alg2==0)
            wangRepeats(grid, cfg.T, cfg.n3, cfg.J);
        else
            metropolis_sweeps(grid, cfg.T, cfg.n3, cfg.J);
    }
    std::vector<std::string> results_temp;
    double result;
    if(calc_e)
        results_temp.push_back(to_string(e_avg/cfg.n2));
    if(calc_m)
        results_temp.push_back(to_string(m_avg/cfg.n2));
    if(calc_cv){
        result = 1.0*(e2_avg/cfg.n2 - e_avg/cfg.n2*e_avg/cfg.n2)*cfg.L*cfg.L/(cfg.T*cfg.T);
        results_temp.push_back(to_string(result));
    }
    if(calc_chi){
        result = 1.0*(m2_avg/cfg.n2 - m_avg/cfg.n2*m_avg/cfg.n2)*cfg.L*cfg.L/(cfg.T*cfg.T);
        results_temp.push_back(to_string(result));
    }
    results = results_temp;
}

void System::measure() {
    if (calc_e || calc_cv){
        double en = calc_E(grid);
        e_avg += en;
        if(calc_cv)
            e2_avg += en*en;
    }
    if (calc_m || calc_chi){
        double mag = calc_m_abs(grid);
        m_avg += mag;
        if(calc_chi)
            m2_avg += mag*mag;
    }
}
