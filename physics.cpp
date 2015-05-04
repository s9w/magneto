#include "physics.h"

double calc_E(std::vector<std::vector<int> > &grid) {
    unsigned int L = grid.size();
    double E = 0.0;
    for(int i = 0; i < L; ++i){
        for (int j = 0; j < L; ++j)
            E += -grid[i][j] * ( grid[i][(j+1)%L] + grid[(i+1)%L][j] );
    }
    return E/(L*L);
}


int calc_dE(std::vector<std::vector<int> >& grid, int idx1, int idx2, const int L) {
    return 2 * grid[idx1][idx2] * (
        grid[idx1][(idx2 + 1) % L] +
        grid[(idx1 + 1) % L][idx2] +
        grid[idx1][(idx2 - 1 + L) % L] +
        grid[(idx1 - 1 + L) % L][idx2]);
}


double calc_m_abs(std::vector<std::vector<int> >& grid) {
    unsigned int L = grid.size();
    int m = 0;
    for (int i = 0; i < L; ++i){
        for (int j = 0; j < L; ++j)
            m += grid[i][j];
    }
    return abs(m)*1.0f/(L*L);
}


double corr_len(std::vector<std::vector<int> >& grid, double T, int avg_n, int evolveRuns,
                const std::function<void(std::vector<std::vector<int> >&, double, int)> &evolve){
    unsigned int L = grid.size();
    unsigned int corr_range = L/2-1;
    std::vector<double> corr_ab(corr_range, 0);
    std::vector<double> corr_a(corr_range, 0);
    std::vector<double> corr_b(corr_range, 0);
    for(int i=0; i<avg_n; ++i) {
        for (int d = 0; d < corr_range; d++) {
            corr_ab[d] += grid[0][0] * grid[0][d];
            corr_a[d] += grid[0][0];
            corr_b[d] += grid[0][d];
        }
        evolve(grid, T, evolveRuns);
    }

    for (int d = 0; d < corr_range; d++) {
        corr_ab[d] = corr_ab[d]/avg_n;
        corr_a[d] = corr_a[d]/avg_n;
        corr_b[d] = corr_b[d]/avg_n;
    }

    double corr_len = 0.0;
    for(int i=0; i<corr_range; ++i)
        corr_len += corr_ab[i] - corr_a[i] * corr_b[i];
    return corr_len/T;
}


double avg_Z(std::vector<std::vector<int> >& grid, double T, int avg_n, int evolveRuns,
             const std::function<void(std::vector<std::vector<int> >&, double, int)> &evolve) {
    unsigned int L = grid.size();
    double Z_avg = 0.0;
    double E = 0.0;
    double beta = 1.0f/T;
    for(int i=0; i<avg_n; ++i) {
        E = 1.0*calc_E(grid)/(L*L);
        Z_avg += exp( -beta*E );
        evolve(grid, T, evolveRuns);
    }
    return Z_avg/(avg_n);
}
