CXX = g++
GMPHOME = /usr/local
GSLHOME = /usr/local
CXXFLAGS = -Wall -ansi -pedantic -std=c++11 -O3 -funroll-loops -pipe

free_fermions: main.o
	$(CXX) main.o -o free_fermions -L$(GMPHOME)/lib -lmpfr -lgmp -L$(GSLHOME)/lib -lgsl -lgslcblas

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -I$(GMPHOME)/include -I$(GSLHOME)/include -c main.cpp

clean:
	rm -f *.o free_fermions
