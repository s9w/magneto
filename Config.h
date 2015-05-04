//
// Created by admin on 29.04.2015.
//

#ifndef ICING_CONFIG_H
#define ICING_CONFIG_H

#include <string>

struct Config{
	unsigned int n1=50, n2=500, n3=5, L=30, threadCount=3;
	double J = 1.0, T=0.0;
    int alg1=0, alg2=1;
	bool recordMain = false;
};

struct LabConfig{
	double TMin, TMax;
	unsigned int TSteps;
	std::string fileEnergy = "", fileMag="", fileCv="", fileChi="", fileStates="";
};

#endif //ICING_CONFIG_H
