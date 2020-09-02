
output: prob2.o 
	g++ prob2.o -o output

prob2.o: prob2.cpp
	g++ -c prob2.cpp

clean:
	rm *.o output


# # compiler & flags
# CXX = g++
# CXXFLAGS = -O -std=c++11
# #CXXFLAGS = -O0 -g -std=c++11
# LIBS = -larmadillo

# # makefile targets
# all : prob2.exe

# prob2.exe : prob2.cpp
# 	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# clean :
# 	\rm -f *.o x.txt T*.txt *.png *.pdf

# realclean : clean
# 	\rm -f *.exe *~


# ####### End of Makefile #######




# Makefile for forward Euler examples
#
# Daniel R. Reynolds
# SMU Mathematics
# Fall 2020

# # compiler & flags
# CXX = g++
# #CXXFLAGS = -O -std=c++11
# CXXFLAGS = -O0 -g -std=c++11
# LIBS = -larmadillo

# # makefile targets
# all : driver_fwd_euler.exe driver_fwd_euler_system.exe \
#       driver_adapt_euler.exe driver_adapt_euler_system.exe

# driver_fwd_euler.exe : driver_fwd_euler.cpp fwd_euler.cpp
# 	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# driver_fwd_euler_system.exe : driver_fwd_euler_system.cpp fwd_euler.cpp
# 	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# driver_adapt_euler.exe : driver_adapt_euler.cpp adapt_euler.cpp
# 	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# driver_adapt_euler_system.exe : driver_adapt_euler_system.cpp adapt_euler.cpp
# 	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)


# # utilities
# clean :
# 	\rm -rf *.txt *.exe *~ *.dSYM
