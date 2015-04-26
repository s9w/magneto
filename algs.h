//
// Created by admin on 23.03.2015.
//

#ifndef _IPYNB_ALGS_H_
#define _IPYNB_ALGS_H_

#include <iostream>
#include <vector>
#include <array>
#include <random>
#include <chrono>
#include <fstream>
#include <functional>
#include <sstream>
#include <deque>
#include <iomanip>

double metropolis_sweeps(std::vector<std::vector<int> >& grid, const int L, double beta, int n);
double metropolis_sweeps_B(std::vector<std::vector<int> >& grid, const int L, double beta, int n, double B);
void wangRun(std::vector<std::vector<int> >& grid, double T);

#endif //_IPYNB_ALGS_H_
