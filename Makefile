all: ising.exe

ising.exe: *.cpp
	g++ -fopenmp -march=native -std=c++11 -O3 ising.cpp System.cpp physics.cpp algs.cpp -o ising.exe

nbstart:
	ipython notebook --profile my