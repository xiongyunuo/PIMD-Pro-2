CPPX = g++
CPP_FLAG = -O2 -std=c++11

PIMD : main.cpp random_x.cpp simulation_x.cpp statistic_x.cpp
	$(CPPX) -o PIMD main.cpp random_x.cpp simulation_x.cpp statistic_x.cpp $(CPP_FLAG)

.PHONY : clean
clean :
	rm PIMD