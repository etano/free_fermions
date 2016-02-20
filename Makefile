CXX = g++
GMPHOME = /usr/local
CXXFLAGS = -Wall -ansi -pedantic -std=c++11 -O3 -funroll-loops -pipe

free_fermions: main.o
	$(CXX) main.o -o free_fermions -L$(GMPHOME)/lib -lmpfr -lgmp

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -I$(GMPHOME)/include -c main.cpp

clean:
	rm -f *.o free_fermions
