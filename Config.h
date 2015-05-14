//
// Created by admin on 29.04.2015.
//

#ifndef ICING_CONFIG_H
#define ICING_CONFIG_H

#include <string>

struct Config{
	unsigned int n1=50, n2=500, n3=5, L=30, threadCount=3;
	int J = 1;
	double T=0.0;
    std::string alg1="metro", alg2="metro", initial="random";
	bool recordMain = false;
};

struct LabConfig{
	double TMin = 0.1, TMax = 4.53;
	bool normalDist = false;
	unsigned int TSteps = 9;
	std::vector<std::string> output_filenames = std::vector<std::string>(6);
};

#endif //ICING_CONFIG_H
