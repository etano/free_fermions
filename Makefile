.PHONY: all clean

.DEFAULT: all

CXX = g++

CXXFLAGS = -Wall -ansi -pedantic -std=c++11 -O3 -funroll-loops -pipe -I/usr/local/include -L/usr/local/lib -lmpfr -lgmp

TARGETS = free_fermions

all: $(TARGETS)

clean:
	rm -f $(TARGETS)

$(TARGETS) : %: main.cpp *.hpp utils/*.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
