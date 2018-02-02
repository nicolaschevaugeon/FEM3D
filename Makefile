CXX = g++
CXXFLAGS = -Wall -g -O3 -fopenmp --std=c++11 
all : main
main : main.cc gmsh_import.cc hyperelastic.h fem_test.h fem.h tensor2_9cm.h gmsh_import.h Makefile
	$(CXX) $(CXXFLAGS) main.cc gmsh_import.cc -o main

test_omp : test_omp.cc
	$(CXX) $(CXXFLAGS) test_omp.cc -o test_omp

