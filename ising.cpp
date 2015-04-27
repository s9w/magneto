#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <sstream>
#include <omp.h>
#include "physics.h"
#include "algs.h"


template <typename T>
std::string to_string(T const & value){
   std::stringstream ss;
   ss << value;
   return ss.str();
}


std::vector<std::vector<int> > genRandomSystem(const unsigned int L, int seedOffset){
	long long int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed1 + seedOffset);
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


std::vector<double> getTemps(double TMin, double TMax, unsigned int TSteps){
	std::vector<double> T_vec;
	double dT = (TMax - TMin)/(TSteps-1);
	if(TSteps==1)
		dT=0;
	for(int i=0; i<TSteps; ++i)
		T_vec.push_back(TMin + i*dT);
	return T_vec;
}


void compIsing(unsigned int L, double TMin, double TMax, unsigned int TSteps, std::string measurement,
			   std::string filename, int avgN, int algReps, const std::function<void(std::vector<std::vector<int> >&, double, int)> &fp){
	std::vector<double> result_vec(TSteps, 0);
	std::string progressStr = std::string(TSteps, '.');
	std::cout << progressStr.c_str() << std::endl;
	auto T_vec = getTemps(TMin, TMax, TSteps);


	#pragma omp parallel for
	for(int i=0; i<TSteps; ++i){
		auto grid = getRelaxedSys(L, T_vec[i]);
		if(measurement == "energy")
			result_vec[i] = avg_En(grid, T_vec[i], avgN, algReps, fp);
		else if(measurement == "mag")
			result_vec[i] = avg_m(grid, T_vec[i], avgN, algReps, fp);
		else if(measurement == "cv")
			result_vec[i] = avg_cv(grid, T_vec[i], avgN, algReps, fp);
		else if(measurement == "chi")
			result_vec[i] = avg_chi(grid, T_vec[i], avgN, algReps, fp);
		else if(measurement == "corrLen")
			result_vec[i] = corr_len(grid, T_vec[i], avgN, algReps, fp);
		else{
			std::cerr << "Invalid measurement!" << std::endl;
		}
		std::cout << ".";
	}
	std::cout << std::endl;

	// write results
	std::ofstream fileOut(filename+".txt");
	for (int i = 0; i < TSteps; ++i)
		fileOut << T_vec[i] << ", " << result_vec[i] << std::endl;
	fileOut.close();
}


void writeEndState(unsigned int L, double T, std::string filename) {
	auto grid = getRelaxedSys(L, T);

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


int main(int argc, char* argv[]){
	unsigned int avgN=500, L=30, TSteps=10, threadN=3, algReps =5;
	double TMin=0.1, TMax=4.53;
	std::string measureStr = "energy", filename = "out", algStr="metr";
	auto evolveFP = metropolis_sweeps;

	if (argc < 2) {
		std::cerr << "Too few arguments. See documentation." << std::endl;
		return 1;
	}

	std::string arg, param;
	for (int i = 1; i < argc; ++i) {
		arg = argv[i];

		param = "-L=";
		if(arg.compare(0, param.length(), param) == 0)
			L = atoi(arg.substr(param.length()).c_str());

		param = "-TSteps=";
		if(arg.compare(0, param.length(), param) == 0)
			TSteps = atoi(arg.substr(param.length()).c_str());

		param = "-avgN=";
		if(arg.compare(0, param.length(), param) == 0)
			avgN = atoi(arg.substr(param.length()).c_str());

		param = "-algN=";
		if(arg.compare(0, param.length(), param) == 0)
			algReps = atoi(arg.substr(param.length()).c_str());

		param = "-threads=";
		if(arg.compare(0, param.length(), param) == 0)
			threadN = atoi(arg.substr(param.length()).c_str());
		omp_set_num_threads(threadN);

		param = "-TMin=";
		if(arg.compare(0, param.length(), param) == 0)
			TMin = atof(arg.substr(param.length()).c_str());

		param = "-TMax=";
		if(arg.compare(0, param.length(), param) == 0)
			TMax = atof(arg.substr(param.length()).c_str());

		param = "-measure=";
		if(arg.compare(0, param.length(), param) == 0)
			measureStr = arg.substr(param.length());

		param = "-alg=";
		if(arg.compare(0, param.length(), param) == 0)
			algStr = arg.substr(param.length());
		if(algStr == "sw")
			evolveFP = wangRepeats;

		param = "-o=";
		if(arg.compare(0, param.length(), param) == 0)
			filename = arg.substr(param.length());
	}

	if(measureStr == "states"){
		auto T_vec = getTemps(TMin, TMax, TSteps);
		for(int i=0; i<T_vec.size(); ++i)
			writeEndState(L, T_vec[i], filename+to_string(i)+".txt");
	}
	else
		compIsing(L, TMin, TMax, TSteps, measureStr, filename, avgN, algReps, evolveFP);

	return 0;
}
