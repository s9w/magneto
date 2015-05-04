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

System::System(Config p_cfg, LabConfig& labCfg) {
    cfg = p_cfg;
    if(!labCfg.fileEnergy.empty())
        calc_e = true;
    if(!labCfg.fileMag.empty())
        calc_m = true;
    if(!labCfg.fileCv.empty())
        calc_cv = true;
    if(!labCfg.fileChi.empty())
        calc_chi = true;
    if(!labCfg.fileStates.empty())
        calc_states = true;
    grid = getRelaxedSys(cfg.L, cfg.T, cfg.J, cfg.n1, cfg.alg1);
}

void System::compute() {
    if (calc_states)
        resultsStates.push_back(to_string(cfg.T) + "\n");

    for(int evolveStep=0; evolveStep<cfg.n2; ++evolveStep) {
        measure();

        // record results
        if(cfg.recordMain)
            recordResults();

        //  evolve
        if (cfg.alg2 == 0)
            wangRepeats(grid, cfg.T, cfg.n3, cfg.J);
        else
            metropolis_sweeps(grid, cfg.T, cfg.n3, cfg.J);
    }
    recordResults();
}

void System::recordResults() {
    double tempResult;
    std::string str="";
    if (calc_e)
        results[0] += to_string(e_avg / cfg.n2);

    if (calc_m)
        results[1] += to_string(m_avg / cfg.n2);

    if (calc_cv) {
        tempResult = 1.0*(e2_avg / cfg.n2 - e_avg / cfg.n2* e_avg / cfg.n2) * cfg.L* cfg.L / (cfg.T* cfg.T);
        results[2] += to_string(tempResult);
    }

    if (calc_chi) {
        tempResult = 1.0*(m2_avg / cfg.n2 - m_avg / cfg.n2* m_avg / cfg.n2) * cfg.L* cfg.L / (cfg.T* cfg.T);
        results[3] += to_string(tempResult);
    }

    if (calc_states) {
        unsigned int L = grid.size();
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                str += ((grid[i][j] == 1) ? "1" : "0");
                if (j < L - 1)
                    str += ",";
            }
            str += "\n";
        }
        resultsStates.push_back(str);
    }
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
