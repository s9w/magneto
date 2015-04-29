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
double avg_En(    std::vector<std::vector<int> >& grid, double T, int avg_n, int,
                  const std::function<void(std::vector<std::vector<int> >&, double, int)> &f, double J);
double corr_len(  std::vector<std::vector<int> >& grid, double T, int avg_n, int,
                  const std::function<void(std::vector<std::vector<int> >&, double, int)> &f, double J);
double avg_cv(    std::vector<std::vector<int> >& grid, double T, int avg_n, int,
                  const std::function<void(std::vector<std::vector<int> >&, double, int)> &f, double J);
double avg_chi(   std::vector<std::vector<int> >& grid, double T, int avg_n, int,
                  const std::function<void(std::vector<std::vector<int> >&, double, int)> &f, double J);
double avg_m(     std::vector<std::vector<int> >& grid, double T, int avg_n, int,
                  const std::function<void(std::vector<std::vector<int> >&, double, int)> &f, double J);
double avg_Z(     std::vector<std::vector<int> >& grid, double T, int avg_n, int,
                  const std::function<void(std::vector<std::vector<int> >&, double, int)> &f, double J);

#endif //_IPYNB_PHYSICS_H_
