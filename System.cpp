#include "System.h"
#include "algs.h"
#include "physics.h"


template<typename T>
std::string to_string(T const & value) {
    std::stringstream ss;
    ss << value;
    return ss.str();
}

System::System(Config p_cfg, LabConfig& labCfg) {
    cfg = p_cfg;
    long long int seed1 = std::chrono::_V2::system_clock::now().time_since_epoch().count();
    gen_metro = std::mt19937(seed1);

    double beta = 1.0f/cfg.T;
    int buffer_offset = 8*abs(cfg.J);
    for(int dE=0; dE<(buffer_offset*2+1); ++dE)
        exp_values.push_back(exp(-(dE-buffer_offset)*beta));

    if(!labCfg.fileEnergy.empty())
        calc_e = true;
    if(!labCfg.fileMag.empty())
        calc_m = true;
    if(!labCfg.fileCv.empty())
        calc_cv = true;
    if(!labCfg.fileChi.empty())
        calc_chi = true;
    if(!labCfg.fileCorr.empty()) {
        calc_corr = true;
        corr_ab.assign(cfg.L/2-1, 0.0);
        corr_a.assign(cfg.L/2-1, 0.0);
        corr_b.assign(cfg.L/2-1, 0.0);
    }
    if(!labCfg.fileStates.empty())
        calc_states = true;
    grid = getRelaxedSys(0);
}

std::vector<std::vector<int> > System::genRandomSystem(int seedOffset){
    unsigned int L = cfg.L;
    long long int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed1 + seedOffset);
    std::uniform_int_distribution <int> dist(0,1);
    std::vector<std::vector<int> > grid(L, std::vector<int>(L));
    for (int i = 0; i < L; ++i){
        for (int j = 0; j < L; ++j){
            grid[i][j] = dist(generator)*2 - 1;
        }
    }
    return grid;
}

std::vector<std::vector<int> > System::getRelaxedSys(int seedOffset) {
    auto grid = genRandomSystem(seedOffset);
    if(cfg.alg1=="metro")
        metropolis_sweeps();
    else if(cfg.alg1=="sw")
        wangRepeats(grid, cfg.T, cfg.n1, cfg.J);
    else
        std::cerr << "unknown alg!" << std::endl;
    return grid;
}

void System::metropolis_sweeps() {
    std::uniform_int_distribution <int> dist_grid(0, cfg.L-1);
    std::uniform_real_distribution <double > dist_one(0.0, 1.0);

    int buffer_offset = 8*abs(cfg.J);
    int flipIdx1, flipIdx2;
    int dE;
    for (int i=0; i < cfg.L*cfg.L*cfg.n3; ++i){
        flipIdx1 = dist_grid(gen_metro);
        flipIdx2 = dist_grid(gen_metro);
        dE = cfg.J * calc_dE(grid, flipIdx1, flipIdx2, cfg.L);
        if (dE <= 0 || (dist_one(gen_metro) < exp_values[dE+buffer_offset]) )
            grid[flipIdx1][flipIdx2] *= -1;
    }
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
        if (cfg.alg2 == "metro")
            metropolis_sweeps();
        else if(cfg.alg2=="sw")
            wangRepeats(grid, cfg.T, cfg.n3, cfg.J);
        else
            std::cerr << "unknown alg!" << std::endl;
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
        tempResult = 1.0*(m2_avg / cfg.n2 - m_avg / cfg.n2* m_avg / cfg.n2) * cfg.L* cfg.L / cfg.T;
        results[3] += to_string(tempResult);
    }

    if (calc_corr) {
        for (int d = 0; d < cfg.L/2-1; d++) {
            corr_ab[d] = corr_ab[d]/cfg.n2;
            corr_a[d] = corr_a[d]/cfg.n2;
            corr_b[d] = corr_b[d]/cfg.n2;
        }
        tempResult = 0.0;
        for(int i=0; i<cfg.L/2-1; ++i)
            tempResult += corr_ab[i] - corr_a[i] * corr_b[i];
        tempResult = tempResult/cfg.T;
        results[4] += to_string(tempResult);
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
    if (calc_corr){
        for (int d = 0; d < cfg.L/2-1; d++) {
            corr_ab[d] += grid[0][0] * grid[0][d];
            corr_a[d] += grid[0][0];
            corr_b[d] += grid[0][d];
        }
    }
}
