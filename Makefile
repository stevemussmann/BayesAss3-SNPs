all : BA3

BA3 : BA3.h main.cpp
	g++ -O3 -Wall main.cpp -o BA3 -lboost_program_options -lgsl -lgslcblas

clean :
	rm *.o BA3
