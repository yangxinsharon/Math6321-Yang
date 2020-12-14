# Makefile for Project Stability-Optimized Explicit Runge-Kutta Methods
# Fourth-Order Low-Storage Runge-Kutta with 12, 13, 14 stages
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
all : BasicTest.exe Brusselator1D.exe ArtificialTest.exe

BasicTest.exe : BasicTest.o LSRK12.o LSRK13.o LSRK14.o erk4.o fwd_euler.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

Brusselator1D.exe : Brusselator1D.o LSRK12.o LSRK13.o LSRK14.o erk4.o fwd_euler.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

ArtificialTest.exe : ArtificialTest.o LSRK12.o LSRK13.o LSRK14.o erk4.o fwd_euler.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^

# utilities
clean :
	\rm -rf *.o *.exe *~ *.dSYM
