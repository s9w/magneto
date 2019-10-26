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

