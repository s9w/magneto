//
// Created by admin on 23.03.2015.
//

#include "algs.h"
#include "physics.h"



double metropolis_sweeps(std::vector<std::vector<int> >& grid, const int L, double beta, int n) {
    long long int seed1 = std::chrono::_V2::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed1);
    std::uniform_int_distribution <int> dist_grid(0,L-1);
    std::uniform_real_distribution <double > dist_one(0.0, 1.0);

    int dE_sweep = 0;
    int flipIdx1, flipIdx2;
    int dE;
    for (int i = 0; i < L*L*n; ++i){
        flipIdx1 = dist_grid(generator);
        flipIdx2 = dist_grid(generator);
        dE = calc_dE(grid, flipIdx1, flipIdx2, L);
        if (dE <= 0 || (dist_one(generator)<exp(-dE*beta)) ){
            grid[flipIdx1][flipIdx2] *= -1;
            dE_sweep += dE;
        }
    }
    return dE_sweep*1.0f/(L*L*n);
}


void wangRun(std::vector<std::vector<int> >& grid, double T) {
    long long int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed1);
    std::uniform_real_distribution<double> dist(0,1);
    unsigned int L = grid.size();
    std::vector<std::vector<int> > discovered(L, std::vector<int>(L, 0));
    std::vector<std::vector<int> > doesBondNorth(L, std::vector<int>(L, 0));
    std::vector<std::vector<int> > doesBondEast(L, std::vector<int>(L, 0));
    double freezeProbability = 1.0 - exp(-2.0f/T);

    for (int i = 0; i < L; ++i){
        for (int j = 0; j < L; ++j) {
            doesBondNorth[i][j] = dist(gen)<freezeProbability;
            doesBondEast[i][j] = dist(gen)<freezeProbability;
        }
    }

    bool flipCluster;
    int x, y, nx, ny;
    for (int i = 0; i < L; ++i){
        for (int j = 0; j < L; ++j) {
            if(! discovered[i][j]){
                flipCluster = dist(gen)<0.5;
                std::deque<std::tuple<int, int>> deq( 1, std::make_tuple(i,j) );
                discovered[i][j] = 1;

                while(!deq.empty()){
                    x = std::get<0>(deq.front());
                    y = std::get<1>(deq.front());

                    nx = x;
                    ny = (y+1)%L;
                    if(grid[x][y]==grid[nx][ny] && discovered[nx][ny]==0 && doesBondNorth[x][y]){
                        deq.push_back(std::make_tuple(nx, ny));
                        discovered[nx][ny] = 1;
                    }

                    nx = (x+1)%L;
                    ny = y;
                    if(grid[x][y]==grid[nx][ny] && discovered[nx][ny]==0 && doesBondEast[x][y]){
                        deq.push_back(std::make_tuple(nx, ny));
                        discovered[nx][ny] = 1;
                    }

                    nx = x;
                    ny = (y-1+L)%L;
                    if(grid[x][y]==grid[nx][ny] && discovered[nx][ny]==0 && doesBondNorth[x][ny]){
                        deq.push_back(std::make_tuple(nx, ny));
                        discovered[nx][ny] = 1;
                    }

                    nx = (x-1+L)%L;
                    ny = y;
                    if(grid[x][y]==grid[nx][ny] && discovered[nx][ny]==0 && doesBondEast[nx][y]){
                        deq.push_back(std::make_tuple(nx, ny));
                        discovered[nx][ny] = 1;
                    }

                    if(flipCluster)
                        grid[x][y] *= -1;
                    deq.pop_front();
                }
            }
        }
    }
}
