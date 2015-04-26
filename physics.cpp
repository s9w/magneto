//
// Created by admin on 23.03.2015.
//

#include "physics.h"
#include "algs.h"

double calc_E(std::vector<std::vector<int> > &grid, const int L) {
    double E = 0.0;
    for(int i = 0; i < L; ++i){
        for (int j = 0; j < L; ++j)
         E += -grid[i][j] * ( grid[i][(j+1)%L] + grid[(i+1)%L][j] );
    }
    return E/(L*L);
}


double avg_En(std::vector<std::vector<int> > grid, double beta, int avg_n) {
    unsigned int L = grid.size();
    double E_avg = 0.0;
    for(int i=0; i<avg_n; ++i) {
        E_avg += calc_E(grid, L);
        wangRun(grid, 1.0f/beta);
    }
    return E_avg/avg_n;
}


int calc_dE(std::vector<std::vector<int> >& grid, int idx1, int idx2, const int L) {
    return 2*grid[idx1][idx2] * ( grid[idx1][(idx2+1)%L] + grid[(idx1+1)%L][idx2] + grid[idx1][(idx2-1+L)%L] + grid[(idx1-1+L)%L][idx2] );
}

double calc_dE_B(std::vector<std::vector<int> >& grid, int idx1, int idx2, const int L, float B) {
    return 2.0*grid[idx1][idx2] * ( grid[idx1][(idx2+1)%L] + grid[(idx1+1)%L][idx2] + grid[idx1][(idx2-1+L)%L] + grid[(idx1-1+L)%L][idx2] ) + 2.0*B*grid[idx1][idx2];
}


double calc_m_abs(std::vector<std::vector<int> > &grid) {
    unsigned int L = grid.size();
    int m = 0;
    for (int i = 0; i < L; ++i){
        for (int j = 0; j < L; ++j)
            m += grid[i][j];
    }
    return abs(m)*1.0f/(L*L);
}

double calc_m(std::vector<std::vector<int> > &grid) {
    unsigned int L = grid.size();
    int m = 0;
    for (int i = 0; i < L; ++i){
        for (int j = 0; j < L; ++j)
            m += grid[i][j];
    }
    return m*1.0f/(L*L);
}


double avg_m(std::vector<std::vector<int> > grid, double beta, int avg_n) {
    double m_avg = 0.0;
    for(int i=0; i<avg_n; ++i) {
        m_avg += calc_m_abs(grid);
//		metropolis_sweeps(grid, L, beta, 10);
        wangRun(grid, 1.0f/beta);
    }
    return m_avg/avg_n;
}


double avg_m2(std::vector<std::vector<int> > grid, double beta, double procTime, double B) {
    auto t0 = std::chrono::high_resolution_clock::now();
    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsed_seconds;

    unsigned int L = grid.size();
    double m_avg = 0.0;
    int avg_n=0;
    while(true) {
        t1 = std::chrono::high_resolution_clock::now();
        elapsed_seconds = (std::chrono::duration_cast <std::chrono::milliseconds > (t1-t0).count())*0.001f;
        if(elapsed_seconds > procTime)
            break;

        m_avg += calc_m(grid);
        metropolis_sweeps_B(grid, L, beta, 15, B);
        avg_n++;
    }
    return m_avg/avg_n;
}


double suszep_bvar(std::vector<std::vector<int> > grid, double beta, double procTime){
    double dB = 0.006;
    return avg_m2(grid, beta, procTime/2.0, dB)-avg_m2(grid, beta, procTime/2.0, 0.0);
}


double corr_len(std::vector<std::vector<int> > grid, double beta, double procTime){
    auto t0 = std::chrono::high_resolution_clock::now();
    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsed_seconds;

    unsigned int L = grid.size();
    unsigned int corr_range = L/2-1;
    std::vector<double> corr_ab(corr_range, 0);
    std::vector<double> corr_a(corr_range, 0);
    std::vector<double> corr_b(corr_range, 0);
    int avg_n=0;
    while(true) {
        t1 = std::chrono::high_resolution_clock::now();
        elapsed_seconds = (std::chrono::duration_cast <std::chrono::milliseconds > (t1-t0).count())*0.001f;
        if(elapsed_seconds > procTime)
            break;

        for (int d = 0; d < corr_range; d++) {
            corr_ab[d] += grid[0][0] * grid[0][d];
            corr_a[d] += grid[0][0];
            corr_b[d] += grid[0][d];
        }
        metropolis_sweeps(grid, L, beta, 15);
        avg_n++;
    }

    for (int d = 0; d < corr_range; d++) {
        corr_ab[d] = corr_ab[d]/avg_n;
        corr_a[d] = corr_a[d]/avg_n;
        corr_b[d] = corr_b[d]/avg_n;
    }

    double corr_len = 0.0;
    for(int i=0; i<corr_range; ++i) {
        corr_len += corr_ab[i] - corr_a[i] * corr_b[i];
    }
    return corr_len;
}


double avg_Z(std::vector<std::vector<int> > grid, double beta, int avg_n) {
    unsigned int L = grid.size();
    double Z_avg = 0.0;
    double E=0.0;
    for(int i=0; i<avg_n; ++i) {
        E = 1.0*calc_E(grid, L)/(L*L);
        Z_avg += exp( -beta*E );
//		wangRun(grid, 1.0f/beta);
        metropolis_sweeps(grid, L, beta, 10);
    }
    return Z_avg/(avg_n);
}


double avg_cv(std::vector<std::vector<int> > grid, double beta, double procTime) {
    auto t0 = std::chrono::high_resolution_clock::now();
    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsed_seconds;

    unsigned int L = grid.size();
    double E_avg = 0.0;
    double E2_avg = 0.0;
    double E;
    int avg_n=0;
    while(true) {
        t1 = std::chrono::high_resolution_clock::now();
        elapsed_seconds = (std::chrono::duration_cast <std::chrono::milliseconds > (t1-t0).count())*0.001f;
        if(elapsed_seconds > procTime)
            break;

        E = calc_E(grid, L);
        E_avg += E;
        E2_avg += E*E;
		metropolis_sweeps(grid, L, beta, 15);
        avg_n++;
    }
    E_avg = E_avg/(avg_n);
    E2_avg = E2_avg/(avg_n);
    double T = 1.0f/beta;
//    std::cout << "avg_n: " << avg_n << std::endl;
    return 1.0*(E2_avg - E_avg*E_avg)*L*L/(T*T);
}


double avg_chi(std::vector<std::vector<int> > grid, double beta, double procTime) {
    auto t0 = std::chrono::high_resolution_clock::now();
    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsed_seconds;

    unsigned int L = grid.size();
    double M_avg = 0.0;
    double M2_avg = 0.0;
    double M;
    int avg_n=0;
    while(true) {
        t1 = std::chrono::high_resolution_clock::now();
        elapsed_seconds = (std::chrono::duration_cast <std::chrono::milliseconds > (t1-t0).count())*0.001f;
        if(elapsed_seconds > procTime)
            break;
        M = calc_m_abs(grid)*L*L;
        M_avg += M;
        M2_avg += M*M;
//        wangRun(grid, 1.0f/beta);
        metropolis_sweeps(grid, L, beta, 15);
        avg_n++;
    }
    M_avg = M_avg /(avg_n);
    M2_avg = M2_avg /(avg_n);
    double T = 1.0f/beta;
//	std::cout << M_avg * M_avg << ", " << M2_avg << ", erg: " << (M2_avg - M_avg * M_avg) << ", T: " << T << std::endl;
    return 1.0*(M2_avg - M_avg * M_avg)/(L*L*T);
}
