all: magneto.exe

magneto.exe: *.cpp *.h
	g++ -fopenmp -march=native -std=c++11 -O3 magneto.cpp System.cpp physics.cpp algs.cpp -o magneto.exe

nbstart:
	ipython notebook --profile my