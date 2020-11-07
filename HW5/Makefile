# Makefile for Prob4 HW5
#
# Sharon Yang
# Math 6321 @ SMU
# Fall 2020

# compiler & flags
CXX = g++
#CXXFLAGS = -O -std=c++11
CXXFLAGS = -O0 -g -std=c++11
INCS = -I ../shared
LIBS = -larmadillo

# executable targets
all : prob4.exe

prob4.exe : prob4.cpp irk.cpp adapt_rkf.cpp	erk4.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(INCS)


# utilities
clean :
	\rm -rf *.txt *.exe *~ *.dSYM
