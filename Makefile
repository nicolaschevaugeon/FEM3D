CXX = g++
CXXFLAGS = -Wall -g -O3 -fopenmp --std=c++11 
all : main test_omp test_one_elem test_law test_tensor
main : main.cc gmsh_import.cc hyperelastic.h fem.h tensor2_9cm.h gmsh_import.h Makefile
	$(CXX) $(CXXFLAGS) main.cc gmsh_import.cc -o main

test_omp : test_omp.cc
	$(CXX) $(CXXFLAGS) test_omp.cc -o test_omp

test_one_elem : test_one_elem.cc hyperelastic.h fem.h tensor2_9cm.h Makefile
	$(CXX) $(CXXFLAGS) test_one_elem.cc -o test_one_elem

test_law : test_law.cc hyperelastic.h tensor2_9cm.h Makefile
	$(CXX) $(CXXFLAGS) test_law.cc -o test_law

test_tensor : test_tensor.cc tensor2_9cm.h Makefile
	$(CXX) $(CXXFLAGS) test_tensor.cc -o test_tensor

clean :
	rm -rf test_tensor test_law test_omp main
