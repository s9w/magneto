#include <sstream>
#include "System.h"
#include "physics.h"
#include <boost/algorithm/string.hpp>

System::System(Config& p_cfg, LabConfig& labCfg) :
        cfg(p_cfg){

    // precalc exp values
    int buffer_offset = 8*abs(cfg.J);
    for(int dE=0; dE < buffer_offset*2+1; ++dE)
        exp_values.push_back(exp(-(dE-buffer_offset)/cfg.T));

    results = std::vector<std::string>(6, "");
    corr_count = 500;
    corr_range = cfg.L/2;
    long long int seed1 = std::chrono::_V2::system_clock::now().time_since_epoch().count();
    rng = std::mt19937(seed1);
    std::uniform_int_distribution<unsigned int> dist_grid(0, cfg.L-1);

    for(int i=0; i< corr_count; ++i){
        correlations.push_back(CorrelationPoint());
        correlations.back().i = dist_grid(rng);
        correlations.back().j = dist_grid(rng);
        correlations.back().corr_a.assign(corr_range, 0.0);
        correlations.back().corr_ab.assign(corr_range, 0.0);
    }

    if(!labCfg.output_filenames[energy].empty())
        calc_e = true;
    if(!labCfg.output_filenames[mag].empty())
        calc_m = true;
    if(!labCfg.output_filenames[cv].empty())
        calc_cv = true;
    if(!labCfg.output_filenames[chi].empty())
        calc_chi = true;
    if(!labCfg.output_filenames[corrfun].empty())
        calc_corrfun = true;
    if(!labCfg.output_filenames[states].empty())
        calc_states = true;

    grid = getRelaxedSys(0);
}

grid_type System::genRandomSystem(int seedOffset){
    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed1 + seedOffset);
    std::uniform_int_distribution <int> dist(0,1);
    grid_type grid(cfg.L, std::vector<int>(cfg.L));
    for (int i = 0; i < cfg.L; ++i){
        for (int j = 0; j < cfg.L; ++j){
            grid[i][j] = dist(generator)*2 - 1;
        }
    }
    return grid;
}

grid_type System::getFileState(std::string filename){
    std::ifstream fileIn(filename);
    std::string line;
    std::vector<std::string> strs;
    grid_type grid(cfg.L, std::vector<int>(cfg.L));
    for (int i = 0; i < cfg.L; ++i){
        std::getline(fileIn, line);
        boost::split(strs, line, boost::is_any_of(","));
        for (int j = 0; j < cfg.L; ++j){
           grid[i][j] = strs[j]=="1"?1:-1;
        }
    }
    return grid;
}

grid_type System::getRelaxedSys(int seedOffset) {
    grid_type new_grid;
    if(cfg.initial == "random")
        new_grid = genRandomSystem(seedOffset);
    else
        new_grid = getFileState(cfg.initial);
    grid = new_grid;

    if(cfg.alg == "metro")
        grid = metropolis_sweeps(grid, cfg.n1);
    else if(cfg.alg == "sw")
        grid = wangRuns(grid, cfg.n1);
    else
        std::cerr << "unknown alg!" << std::endl;
    return grid;
}

grid_type& System::metropolis_sweeps(grid_type& lattice, unsigned int n) {
    std::uniform_int_distribution <int> dist_grid(0, cfg.L-1);
    std::uniform_real_distribution <double > dist_one(0.0, 1.0);

    int buffer_offset = 8*abs(cfg.J);
    int flipIdx1, flipIdx2;
    int dE;
    for (int i=0; i < cfg.L*cfg.L*n; ++i){
        flipIdx1 = dist_grid(rng);
        flipIdx2 = dist_grid(rng);
        dE = cfg.J * calc_dE(lattice, flipIdx1, flipIdx2, cfg.L);
        if (dE <= 0 || (dist_one(rng) < exp_values[dE+buffer_offset]) )
            lattice[flipIdx1][flipIdx2] *= -1;
    }
    return lattice;
}

grid_type& System::wangRuns(grid_type& lattice, const unsigned int n) {
    std::uniform_real_distribution<double> dist_one(0,1);
    const double freezeProbability = 1.0 - exp(-2.0f*cfg.J/cfg.T);

    grid_type discovered;
    grid_type doesBondNorth;
    grid_type doesBondEast;

    bool flipCluster;
    int x, y, nx, ny;
    for(int run=0; run<n; ++run) {
        discovered.assign(cfg.L, std::vector<int>(cfg.L, 0));
        doesBondNorth.assign(cfg.L, std::vector<int>(cfg.L, 0));
        doesBondEast.assign(cfg.L, std::vector<int>(cfg.L, 0));

        for (int i = 0; i < cfg.L; ++i) {
            for (int j = 0; j < cfg.L; ++j) {
                doesBondNorth[i][j] = dist_one(rng) < freezeProbability;
                doesBondEast[i][j] = dist_one(rng) < freezeProbability;
            }
        }


        for (int i = 0; i < cfg.L; ++i) {
            for (int j = 0; j < cfg.L; ++j) {
                if (!discovered[i][j]) {
                    flipCluster = dist_one(rng) < 0.5;
                    std::deque<std::tuple<int, int>> deq(1, std::make_tuple(i, j));
                    discovered[i][j] = 1;

                    while (!deq.empty()) {
                        x = std::get<0>(deq.front());
                        y = std::get<1>(deq.front());

                        nx = x;
                        ny = (y + 1) % cfg.L;
                        if (lattice[x][y] == lattice[nx][ny] && discovered[nx][ny] == 0 && doesBondNorth[x][y]) {
                            deq.push_back(std::make_tuple(nx, ny));
                            discovered[nx][ny] = 1;
                        }

                        nx = (x + 1) % cfg.L;
                        ny = y;
                        if (lattice[x][y] == lattice[nx][ny] && discovered[nx][ny] == 0 && doesBondEast[x][y]) {
                            deq.push_back(std::make_tuple(nx, ny));
                            discovered[nx][ny] = 1;
                        }

                        nx = x;
                        ny = (y - 1 + cfg.L) % cfg.L;
                        if (lattice[x][y] == lattice[nx][ny] && discovered[nx][ny] == 0 && doesBondNorth[x][ny]) {
                            deq.push_back(std::make_tuple(nx, ny));
                            discovered[nx][ny] = 1;
                        }

                        nx = (x - 1 + cfg.L) % cfg.L;
                        ny = y;
                        if (lattice[x][y] == lattice[nx][ny] && discovered[nx][ny] == 0 && doesBondEast[nx][y]) {
                            deq.push_back(std::make_tuple(nx, ny));
                            discovered[nx][ny] = 1;
                        }

                        if (flipCluster)
                            lattice[x][y] *= -1;
                        deq.pop_front();
                    }
                }
            }
        }
    }
    return lattice;
}

void System::compute() {
    if (calc_states)
        results[states] += to_string(cfg.T) + "\n";

    for(int evolveStep=0; evolveStep<cfg.n2; ++evolveStep) {
        measure();

        if(cfg.recordMain)
            recordResults();

        //  evolve
        if (cfg.alg == "metro")
            grid = metropolis_sweeps(grid, cfg.n3);
        else if(cfg.alg == "sw")
            grid = wangRuns(grid, cfg.n3);
        else
            std::cerr << "unknown alg!" << std::endl;
    }
    recordResults();
}

void System::recordResults() {
    double tempResult;
    if (calc_e)
        results[energy] += ", " + to_string(e_avg / cfg.n2);

    if (calc_m)
        results[mag] += ", " + to_string(m_avg / cfg.n2);

    if (calc_cv) {
        tempResult = 1.0*(e2_avg / cfg.n2 - e_avg / cfg.n2* e_avg / cfg.n2) * cfg.L* cfg.L / (cfg.T* cfg.T);
        results[cv] += ", " + to_string(tempResult);
    }

    if (calc_chi) {
        tempResult = 1.0*(m2_avg / cfg.n2 - m_avg / cfg.n2* m_avg / cfg.n2) * cfg.L* cfg.L / cfg.T;
        results[chi] += ", " + to_string(tempResult);
    }

    if (calc_corrfun) {
        const unsigned int d_limit = cfg.L/2;
        std::vector<double> sigma_ij = std::vector<double>(d_limit, 0.0);
        std::vector<double> sigma_i = std::vector<double>(d_limit, 0.0);
        std::vector<double> sigma_j = std::vector<double>(d_limit, 0.0);
        std::vector<double> G = std::vector<double>(d_limit, 0.0);

        for(auto &corr : correlations){
            for (int d = 0; d < corr_range; d++) {
                sigma_ij[d] += corr.corr_ab[d];
                sigma_i[d]  += corr.corr_a[0];
                sigma_j[d]  += corr.corr_a[d];
            }
        }

        for (int d = 0; d < d_limit; d++) {
            G[d] = sigma_ij[d]/(cfg.n2* corr_count) - sigma_i[d]/(cfg.n2* corr_count)*sigma_j[d]/(cfg.n2*corr_count);
            results[corrfun] += ", " + to_string(G[d]);
        }
    }

    if (calc_states) {
        std::string str="";
        unsigned int L = grid.size();
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                str += ((grid[i][j] == 1) ? "1" : "0");
                if (j < L - 1)
                    str += ",";
            }
            str += "\n";
        }
        results[states] += str+"\n";
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
    if (calc_corrfun){
        for(auto &corr : correlations){
            for (int d = 0; d < corr_range; d++) {
                corr.corr_a[d]  += grid[corr.i][(corr.j+d)%cfg.L];
                corr.corr_ab[d] += grid[corr.i][corr.j] * grid[corr.i][(corr.j+d)%cfg.L];
            }
        }
    }
}
