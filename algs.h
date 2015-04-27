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

void metropolis_sweeps(std::vector<std::vector<int> >& grid, double T, int n);
void wangRun(std::vector<std::vector<int> >& grid, double T);
void wangRepeats(std::vector<std::vector<int> >& grid, double T, int n);

#endif //_IPYNB_ALGS_H_
