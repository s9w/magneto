//
// Created by admin on 23.03.2015.
//

#ifndef _IPYNB_PHYSICS_H_
#define _IPYNB_PHYSICS_H_

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

double calc_E(    std::vector<std::vector<int> >& grid);
int    calc_dE(   std::vector<std::vector<int> >& grid, int idx1, int idx2, const int L);
double calc_m_abs(std::vector<std::vector<int> >& grid);

#endif //_IPYNB_PHYSICS_H_
