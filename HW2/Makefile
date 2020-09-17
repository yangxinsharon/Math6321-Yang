# Makefile for forward Euler examples
#
# Daniel R. Reynolds
# SMU Mathematics
# Fall 2020

# compiler & flags
CXX = g++
#CXXFLAGS = -O -std=c++11
CXXFLAGS = -O0 -g -std=c++11
INCS = -I ../shared
LIBS = -larmadillo

# makefile targets
all : prob1.exe prob2.exe

prob1.exe : prob1.cpp fwd_euler.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(INCS)

prob2.exe : prob2.cpp fwd_euler.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(INCS)


# utilities
clean :
	\rm -rf *.txt *.exe *~ *.dSYM
