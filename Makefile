all : a.out
a.out : main.cc gmsh_import.cc Makefile
	g++ -Wall -O3 -fopenmp -std=c++11 main.cc gmsh_import.cc
