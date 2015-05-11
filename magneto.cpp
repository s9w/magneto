#include <iostream>
#include <chrono>
#include <vector>
#include <random>
#include <fstream>
#include <sstream>
#include <omp.h>
#include "physics.h"
#include "System.h"
#include "helpers.h"

void checkParam(int argc, char* argv[], Config& cfg, LabConfig& labCfg){
    std::string key, value;
    int eqPos;
    std::string argument;

    for (int i = 1; i < argc; ++i) {
        argument = argv[i];
        eqPos = argument.find("=");
        key = argument.substr(1, eqPos-1);
        value = argument.substr(eqPos+1);

        if (key == "L")
            cfg.L = atoi(value.c_str());
        else if (key == "N1")
            cfg.n1 = atoi(value.c_str());
        else if (key == "N2")
            cfg.n2 = atoi(value.c_str());
        else if (key == "N3")
            cfg.n3 = atoi(value.c_str());
        else if (key == "threads")
            cfg.threadCount = atoi(value.c_str());
        else if (key == "J")
            cfg.J = atoi(value.c_str());
        else if (key == "alg1")
            cfg.alg1 = value;
        else if (key == "alg2")
            cfg.alg2 = value;
        else if (key == "record")
            cfg.recordMain = value == "main";
        else if (key == "dist")
            labCfg.normalDist = value == "normal";
        else if (key == "en")
            labCfg.fileEnergy = value;
        else if (key == "mag")
            labCfg.fileMag = value;
        else if (key == "cv")
            labCfg.fileCv = value;
        else if (key == "chi")
            labCfg.fileChi = value;
        else if (key == "corr")
            labCfg.fileCorr = value;
        else if (key == "states")
            labCfg.fileStates = value;
        else if (key == "TSteps")
            labCfg.TSteps = atoi(value.c_str());
        else if (key == "TMin")
            labCfg.TMin = atof(value.c_str());
        else if (key == "TMax")
            labCfg.TMax = atof(value.c_str());
        else
            std::cerr << "Unknown parameter: " << argument << std::endl;
    }
}

int main(int argc, char* argv[]){
    if (argc < 2) {
        std::cerr << "Too few arguments. See documentation." << std::endl;
        return 1;
    }

    Config cfg;
    LabConfig labCfg;
    checkParam(argc, argv, cfg, labCfg);
    omp_set_num_threads(cfg.threadCount);

    std::vector<double> temps;
    if(labCfg.normalDist)
        temps = normalSpace(labCfg.TMin, labCfg.TMax, labCfg.TSteps, 2.269, 1.0);
    else
	    temps = linspace(labCfg.TMin, labCfg.TMax, labCfg.TSteps);
    std::vector<System> systems;
    for(double T : temps){
		cfg.T = T;
        systems.push_back(System(cfg, labCfg));
    }

    std::string progressStr = std::string(labCfg.TSteps, '.');
    std::cout << progressStr.c_str() << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for
    for(unsigned int i=0; i<systems.size(); ++i){
        systems[i].compute();
        std::cout << ".";
    }
    std::cout << std::endl;

    std::vector<std::string> filenames = {labCfg.fileEnergy, labCfg.fileMag, labCfg.fileCv, labCfg.fileChi, labCfg.fileCorr};
    std::ofstream fileOut;

    for(unsigned int i=0; i<filenames.size(); ++i){
        if(! filenames[i].empty()){
            fileOut.open(filenames[i]+".txt");
            for (auto& sys : systems)
                fileOut << to_string(sys.cfg.T) << ", " << sys.results[i] << std::endl;
            fileOut.close();
        }
    }

    if(!labCfg.fileStates.empty()){
        for (int i=0; i<systems.size(); ++i){
            fileOut.open(labCfg.fileStates+to_string(i)+".txt");
            for(int j=0; j<systems[i].resultsStates.size(); ++j)
                fileOut << systems[i].resultsStates[j] << std::endl;
            fileOut.close();
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    float secs = (std::chrono::duration_cast <std::chrono::milliseconds > (t1-t0).count())*0.001f;
    std::cout << "runtime: " << secs << "s" << std::endl;

	return 0;
}
