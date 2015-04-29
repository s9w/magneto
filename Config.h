//
// Created by admin on 29.04.2015.
//

#ifndef ICING_CONFIG_H
#define ICING_CONFIG_H

#include <string>

struct Config{
	unsigned int n1=50, n2=500, n3=5, L=30, threadCount=3;
//	std::vector<double> T_vec;
	double J = 1.0, T=0.0;
    int alg1=0, alg2=1;
	std::string fileEnergy = "", fileMag="", fileCv="", fileChi="";
//	auto algEvolve = metropolis_sweeps;
};

#endif //ICING_CONFIG_H
