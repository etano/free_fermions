.PHONY: all clean

.DEFAULT: all

CXX = g++

GMPFLAGS = -I/usr/local/include -L/usr/local/lib -lmpfr -lgmp
CXXFLAGS = -Wall -ansi -pedantic -std=c++11 -O3 -funroll-loops -pipe $(GMPFLAGS)

TARGETS = free_fermions

all: $(TARGETS)

clean:
	rm -f $(TARGETS)

$(TARGETS) : %: main.cpp *.hpp utils/*.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
