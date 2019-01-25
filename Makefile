all:
	g++ -std=c++11 -Wall -O3 -fopenmp main.cpp -o trace

clean:
	rm -f trace out.ppm
