all: magneto.exe

magneto.exe: *.cpp *.h
	g++ -fopenmp -march=native -std=c++11 -O3 helpers.cpp magneto.cpp System.cpp physics.cpp -o magneto.exe

debug: *.cpp *.h
	g++ -fopenmp -g -std=c++11 helpers.cpp magneto.cpp System.cpp physics.cpp -o magneto.exe

testing.exe: *.cpp *.h
	g++ -std=c++11 helpers.cpp testing.cpp -o testing.exe

nbstart:
	ipython notebook --profile my