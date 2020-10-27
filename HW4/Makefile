# Makefile for Runge-Kutta-Fehlberg
#
# Sharon Yang
# SMU Mathematics
# Fall 2020

# compiler & flags
CXX = g++
#CXXFLAGS = -O -std=c++11
CXXFLAGS = -O0 -g -std=c++11
LIBS = -larmadillo

# makefile targets
all : prob2.exe prob3.exe

prob2.exe : prob2.cpp adapt_rkf.cpp adapt_euler.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(INCS)

prob3.exe : prob3.cpp adapt_rkf.cpp erk4.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(INCS)


# utilities
clean :
	\rm -rf *.txt *.exe *~ *.dSYM
