CXX = g++
CXXFLAGS = -Wall -g -O3 -fopenmp --std=c++11 
all : main test_omp test_one_elem
main : main.cc gmsh_import.cc hyperelastic.h fem.h tensor2_9cm.h gmsh_import.h Makefile
	$(CXX) $(CXXFLAGS) main.cc gmsh_import.cc -o main

test_omp : test_omp.cc
	$(CXX) $(CXXFLAGS) test_omp.cc -o test_omp


test_one_elem : test_one_elem.cc hyperelastic.h fem.h tensor2_9cm.h Makefile
	$(CXX) $(CXXFLAGS) test_one_elem.cc -o test_one_elem


