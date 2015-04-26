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
#include "physics.h"
#include "algs.h"


template <typename T>
std::string to_string(T const & value){
   std::stringstream ss;
   ss << value;
   return ss.str();
}


void printGrid(std::vector<std::vector<int> >& grid, const int L){
	for (int i = 0; i < L; ++i){
		for (int j = 0; j < L; ++j)
			std::cout << ((grid[i][j] == 1)?"1":" ");
		std::cout << std::endl;
	}
}


std::vector<std::vector<int> > genRandomSystem(const unsigned int L, int seedOffset){
	long long int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed1+seedOffset);
	std::uniform_int_distribution <int> dist(0,1);
	std::vector<std::vector<int> > grid(L, std::vector<int>(L));
	for (int i = 0; i < L; ++i){
		for (int j = 0; j < L; ++j){
			grid[i][j] = dist(generator)*2 - 1;
		}
	}
	return grid;
}


std::vector<std::vector<int> > getRelaxedSys(const unsigned int L, const double T, int seedOffset=0){
	auto grid = genRandomSystem(L, seedOffset);
	for(int i=0; i<500; ++i)
		wangRun(grid, T);
	return grid;
}

std::vector<std::vector<int> > getRelaxedSys_B(const unsigned int L, const double T, double B, int seedOffset=0){
	auto grid = genRandomSystem(L, seedOffset);
	for(int i=0; i<500; ++i)
		metropolis_sweeps_B(grid, L, 1.0/T, 15, B);

	return grid;
}


void writeFunOfT(unsigned int L, double Tmin, double Tmax, int steps, const std::function<double(std::vector<std::vector<int> >, double, double)> &fp, std::string filename, double procTime){
	std::vector<double> result_vec(steps, 0);
	std::string progStr = std::string(steps, '.');
	std::cout << progStr.c_str() << std::endl;

	std::vector<double> T_vec;
	double dT = (Tmax - Tmin)/(steps-1);
	for(int i=0; i<steps; ++i)
		T_vec.push_back(Tmin + i*dT);

#pragma omp parallel for num_threads(3)
	for(int i=0; i<steps; ++i){
		auto grid = getRelaxedSys(L, T_vec[i]);
		result_vec[i] = fp(grid, 1.0/ T_vec[i], procTime);
		std::cout << ".";
	}
	std::cout << std::endl;

	// write results
	std::ofstream fileOut(filename);
	for (int i = 0; i < steps; ++i)
		fileOut << T_vec[i] << ", " << result_vec[i] << std::endl;
	fileOut.close();
}


std::vector<double> corr_fun(std::vector<std::vector<int> > grid, double beta, int avg_n){
	unsigned int L = grid.size();
	unsigned int corr_range = L/2-1;
	std::vector<double> corr_ab(corr_range, 0);
	std::vector<double> corr_a(corr_range, 0);
	std::vector<double> corr_b(corr_range, 0);
	double a2=0.0;
	for(int k=0; k<avg_n; ++k) {
		for (int d = 0; d < corr_range; d++) {
			corr_ab[d] += grid[0][0] * grid[0][d];
			corr_a[d] += grid[0][0];
			corr_b[d] += grid[0][d];
		}
		a2 += grid[0][0];
		metropolis_sweeps(grid, L, beta, 5);
//		wangRun(grid, 1.0f/beta);
	}

	for (int d = 0; d < corr_range; d++) {
		corr_ab[d] = corr_ab[d]/avg_n;
		corr_a[d] = corr_a[d]/avg_n;
		corr_b[d] = corr_b[d]/avg_n;

	}
	a2 = a2/avg_n;

	std::vector<double> corr(corr_range);
	double term1, term2;
	for(int i=0; i<corr_range; ++i) {
		term1 = corr_ab[i];
		term2 = corr_a[i] * corr_b[i];

//		if(fabs(1.0f/beta-2.2<0.01))
//			std::cout << "T: " << 1.0f/beta << ", d: " << i << ", x: " << std::setw(10) << term1 <<  ", "  << std::setw(10) << term2 << ", " << corr_a[i] << ", " << corr_b[i] << ", a2: " << a2 << std::endl;
		corr[i] = (term1-term2);
	}
	return corr;
}


//void writeCorr(std::vector<config>& configs, int avg_n){
//	std::string progStr = std::string(configs.size(), '.');
//	std::cout << progStr.c_str() << std::endl;
//	int omp_get_thread_num();
//	std::vector<std::vector<double> > results_vec(configs.size());
//	#pragma omp parallel for num_threads(3)
//	for(int i=0; i<configs.size(); ++i){
//		auto grid = getRelaxedSys(configs[i].L, configs[i].T, configs[i].warmupTime, omp_get_thread_num());
//		results_vec[i] = corr_fun(grid, 1.0f/configs[i].T, avg_n);
//		std::cout << ".";
//	}
//
//	// write result header
//	std::ofstream fileOut("ising_corr0.txt");
//	fileOut << "#dist";
//	for (auto& cfg : configs)
//		fileOut << ", " << cfg.T;
//	fileOut << std::endl;
//
//	// then results
//	for(int iDist=0; iDist<results_vec[0].size(); ++iDist) {
//		fileOut << iDist;
//		for (int i = 0; i < configs.size(); ++i)
//			fileOut << ", " << results_vec[i][iDist];
//		fileOut << std::endl;
//	}
//	fileOut.close();
//}

void writeEndState(int L, double T, std::string filename) {
	auto grid = getRelaxedSys(L, T);
//	metropolis_sweeps(grid, L, 1.0/T, 100);

	std::ofstream fileOut(filename);
	for (int i = 0; i < L; ++i){
		for (int j = 0; j < L; ++j) {
			fileOut << ((grid[i][j] == 1) ? "1" : "0");
			if(j<L-1)
				fileOut << ",";
		}
		fileOut << std::endl;
	}
	fileOut.close();
}

int main(){
	double Tc = 2.0f/log(1+sqrt(2));
	int steps, avg_n, L;
	double Tmin, Tmax;

	// Energy
//	std::vector<unsigned int> L_vec = {40, 80};
//	Tmin = 0.1;
//	Tmax = 2*Tc;
//	steps = 30;
//	avg_n = 3000;
//	for(int i=0; i<L_vec.size(); ++i)
//		writeFunOfT(L_vec[i], Tmin, Tmax, steps, avg_En, "ising_e"+ to_string(i)+".txt", avg_n);


	// Magnetization
//	std::vector<unsigned int> L_vec = {40, 80};
//	Tmin = 0.1;
//	Tmax = 2*Tc;
//	steps = 30;
//	avg_n = 1000;
//	for(int i=0; i<L_vec.size(); ++i)
//		writeFunOfT(L_vec[i], Tmin, Tmax, steps, avg_m, "ising_m"+to_string(i)+".txt", avg_n);

	// Magnetization with external B-field
	std::vector<unsigned int> L_vec = {32};
	Tmin = 0.3;
	Tmax = 1.3*Tc;
	steps = 15;
	double procTime = 200.0;
	for(int i=0; i<L_vec.size(); ++i)
		writeFunOfT(L_vec[i], Tmin, Tmax, steps, suszep_bvar, "ising_hvar"+ to_string(i)+".txt", procTime/steps);

	// Correlation length
//	std::vector<unsigned int> L_vec = {32};
//	Tmin = 0.0;
//	Tmax = 2.0*Tc;
//	steps = 24;
//	double procTime = 200.0;
//	for(int i=0; i<L_vec.size(); ++i)
//		writeFunOfT(L_vec[i], Tmin, Tmax, steps, corr_len, "ising_corrlen"+ to_string(i)+".txt", procTime/steps);


	// Heat capacity
//	std::vector<unsigned int> L_vec = {16, 80};
//	Tmin = 2.0;
//	Tmax = 2.6;
//	steps = 24;
//	double procTime = 420.0;
//	for(int i=0; i<L_vec.size(); ++i)
//		writeFunOfT(L_vec[i], Tmin, Tmax, steps, avg_cv, "ising_cv"+ to_string(i)+".txt", procTime/steps);


	// Susceptibility
//	std::vector<unsigned int> L_vec = {16, 32, 64};
//	Tmin = 2.0;
//	Tmax = 2.6;
//	steps = 15;
//	double procTime = 120.0;
////	avg_n = 5000;
//	for(int i=0; i<L_vec.size(); ++i)
//		writeFunOfT(L_vec[i], Tmin, Tmax, steps, avg_chi, "ising_chi"+ to_string(i)+".txt", procTime/steps);



//	 Korrelation
//	steps = 33;
//	totaltime = 60.0f;
//	L = 60;
//	cfg0 = {L, 2.0, totaltime/steps};
//	cfg1 = {L, 2.6, totaltime/steps};
//	configs = createConfigs(cfg0, cfg1, steps);
//	writeCorr(configs, 10000);

//	const unsigned int L = 50;
//	const double T = 2.0;
//	auto grid = getRelaxedSys(L, T, 7);
//	auto corr = corr_fun(grid, L, 1.0f/T, 500);
//	std::ofstream fileOut("ising_corr0.txt");
//	for (int i = 0; i < corr.size(); ++i)
//		fileOut << i << ", " << corr[i] << std::endl;
//	fileOut.close();


	// states on different Temperatues
//	L=200;
//	std::vector<double> T_vec = {0.8*Tc, 0.9*Tc, 0.95*Tc, 1.0*Tc, 1.1*Tc, 2.0*Tc};
//	for(int i=0; i<T_vec.size(); ++i)
//		writeEndState(L, T_vec[i], "ising_state"+ to_string(i)+".txt");

	return 0;
}

