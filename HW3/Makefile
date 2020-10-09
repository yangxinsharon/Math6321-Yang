# Makefile for HW3
#
# Sharon Yang
# SMU Mathematics
# Fall 2020

# compiler & flags
CXX = g++
#CXXFLAGS = -O -std=c++11
CXXFLAGS = -O0 -g -std=c++11
LIBS = -larmadillo

# executable targets
all : prob3.exe

prob3.exe : prob3.cpp GeneralizedMidpoint.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(INCS)


# utilities
clean :
	\rm -rf *.txt *.exe *~ *.dSYM
