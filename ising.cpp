#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <sstream>
#include <omp.h>

#include "physics.h"
#include "algs.h"
#include "System.h"


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




std::vector<std::vector<int> > getRelaxedSys(const unsigned int L, const double T, double J, unsigned int n1, int alg, int seedOffset=0) {
    auto grid = genRandomSystem(L, seedOffset);
    if(alg==0)
        wangRepeats(grid, T, n1, J);
    else
        metropolis_sweeps(grid, T, n1, J);
    return grid;
}



//void writeEndState(unsigned int L, double T, std::string filename) {
//	auto grid = getRelaxedSys(L, T, J);
//
//	std::ofstream fileOut(filename);
//	for (int i = 0; i < L; ++i){
//		for (int j = 0; j < L; ++j) {
//			fileOut << ((grid[i][j] == 1) ? "1" : "0");
//			if(j<L-1)
//				fileOut << ",";
//		}
//		fileOut << std::endl;
//	}
//	fileOut.close();
//}
//
//
//void writeMovie(unsigned int L, double T, std::string filename, int frames,
//				const std::function<void(std::vector<std::vector<int> >&, double, int)>& algFP, int algN) {
//	auto grid = getRelaxedSys(L, T);
//	std::ofstream fileOut(filename);
//	for(int frame=0; frame<frames; ++frame) {
//		for (int i = 0; i < L; ++i) {
//			for (int j = 0; j < L; ++j) {
//				fileOut << ((grid[i][j] == 1) ? "1" : "0");
//				if (j < L - 1)
//					fileOut << ",";
//			}
//			fileOut << std::endl;
//		}
//		fileOut << std::endl;
//		algFP(grid, T, algN);
//	}
//	fileOut.close();
//}


std::vector<double> getTemps(double TMin, double TMax, unsigned int TSteps){
    std::vector<double> T_vec;
    double dT = (TMax - TMin)/(TSteps-1);
    if(TSteps==1)
        dT=0;
    for(int i=0; i<TSteps; ++i)
        T_vec.push_back(TMin + i*dT);
    return T_vec;
}

template<typename T>
std::string to_string(T const & value);

int main(int argc, char* argv[]){
	unsigned int TSteps=10;
	double TMin=0.1, TMax=4.53;
	Config cfg;

	if (argc < 2) {
		std::cerr << "Too few arguments. See documentation." << std::endl;
		return 1;
	}

	std::string arg, param;
	for (int i = 1; i < argc; ++i) {
		arg = argv[i];

		param = "-L=";
		if(arg.compare(0, param.length(), param) == 0)
			cfg.L = atoi(arg.substr(param.length()).c_str());

		param = "-N2=";
		if(arg.compare(0, param.length(), param) == 0)
			cfg.n2 = atoi(arg.substr(param.length()).c_str());

		param = "-N3=";
		if(arg.compare(0, param.length(), param) == 0)
			cfg.n3 = atoi(arg.substr(param.length()).c_str());

		param = "-threads=";
		if(arg.compare(0, param.length(), param) == 0)
			cfg.threadCount = atoi(arg.substr(param.length()).c_str());
		omp_set_num_threads(cfg.threadCount);

		param = "-J=";
		if(arg.compare(0, param.length(), param) == 0)
			cfg.J = atof(arg.substr(param.length()).c_str());

        param = "-alg1=";
        if(arg.compare(0, param.length(), param) == 0)
            cfg.alg1 = (arg.substr(param.length())=="metro") ? 0 : 1;

		param = "-alg2=";
		if(arg.compare(0, param.length(), param) == 0)
			cfg.alg2 = (arg.substr(param.length())=="metro") ? 0 : 1;

		param = "-en=";
		if(arg.compare(0, param.length(), param) == 0)
			cfg.fileEnergy = arg.substr(param.length());

        param = "-mag=";
        if(arg.compare(0, param.length(), param) == 0)
            cfg.fileMag = arg.substr(param.length());

        param = "-cv=";
        if(arg.compare(0, param.length(), param) == 0)
            cfg.fileCv = arg.substr(param.length());

        param = "-chi=";
        if(arg.compare(0, param.length(), param) == 0)
            cfg.fileChi = arg.substr(param.length());

		// Temperatures
		param = "-TSteps=";
		if(arg.compare(0, param.length(), param) == 0)
			TSteps = atoi(arg.substr(param.length()).c_str());

		param = "-TMin=";
		if(arg.compare(0, param.length(), param) == 0)
			TMin = atof(arg.substr(param.length()).c_str());

		param = "-TMax=";
		if(arg.compare(0, param.length(), param) == 0)
			TMax = atof(arg.substr(param.length()).c_str());
	}
	auto temps = getTemps(TMin, TMax, TSteps);
    std::vector<System> systems;
    for(double T : temps){
		cfg.T = T;
        systems.push_back(System(cfg));
    }

    std::string progressStr = std::string(TSteps, '.');
    std::cout << progressStr.c_str() << std::endl;
    #pragma omp parallel for
    for(unsigned int i=0; i<systems.size(); ++i){
        systems[i].compute();
        std::cout << ".";
    }

    std::vector<std::string> filenames = {cfg.fileEnergy, cfg.fileMag, cfg.fileCv, cfg.fileChi};
    std::ofstream fileOut;
    for(unsigned int i=0; i<4; ++i){
        if(! filenames[i].empty()){
            fileOut.open(filenames[i]+".txt");
            for (auto& sys : systems)
                fileOut << to_string(sys.cfg.T) << ", " << sys.results[i] << std::endl;
            fileOut.close();
        }
    }



//	if(measureStr == "states"){
//		auto T_vec = getTemps(TMin, TMax, TSteps);
//		for(int i=0; i<T_vec.size(); ++i)
//			writeEndState(L, T_vec[i], filename+to_string(i)+".txt");
//	}
//	else if(measureStr == "movie"){
//		writeMovie(L, TMin, filename+".txt", frameCount, algFP, algReps);
//	}
//	else

	return 0;
}
