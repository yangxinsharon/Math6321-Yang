#!/usr/bin/env python
#
# Plotting script for 1D Brusselator
#
# Sharon Yang
# Math 6321 @ SMU
# Fall 2020

# imports
from numpy import *
from matplotlib.pyplot import *

# First, load our saved data files
t = loadtxt('brusselator_tspan.txt')
# h(0) = 0.1
Y12_0 = loadtxt('Brusselator_Y12_0.txt')
Y13_0 = loadtxt('Brusselator_Y13_0.txt')
Y14_0 = loadtxt('Brusselator_Y14_0.txt')
Yerk4_0 = loadtxt('Brusselator_Yerk4_0.txt')
Yfe_0 = loadtxt('Brusselator_Yfe_0.txt')

## h(1) = 0.05
#Y12_1 = loadtxt('Brusselator_Y12_1.txt')
#Y13_1 = loadtxt('Brusselator_Y13_1.txt')
#Y14_1 = loadtxt('Brusselator_Y14_1.txt')
#Yerk4_1 = loadtxt('Brusselator_Yerk4_1.txt')
#Yfe_1 = loadtxt('Brusselator_Yfe_1.txt')
#
## h(2) = 0.01
#Y12_2 = loadtxt('Brusselator_Y12_2.txt')
#Y13_2 = loadtxt('Brusselator_Y13_2.txt')
#Y14_2 = loadtxt('Brusselator_Y14_2.txt')
#Yerk4_2 = loadtxt('Brusselator_Yerk4_2.txt')
#Yfe_2 = loadtxt('Brusselator_Yfe_2.txt')

# LSRK 12
# create plot 
figure(1)
plot(t[:],Y12_0[0,:],label = "species X", color='red')
plot(t[:],Y12_0[1,:],label = "species Y", color='blue')
xlabel('$time$')
ylabel('$concentration$')
legend(loc="upper right")
title('LSRK12')
savefig('Bru_Y12.png')

## create plot 
#figure(2)
#plot(t[:],Y12_1[0,:])
#plot(t[:],Y12_1[1,:])
#xlabel('$time$')
#ylabel('$concentration$')
#title('LSRK12 at h = 0.05')
#savefig('Bru_Y12_1.png')
#
## create plot 
#figure(3)
#plot(t[:],Y12_2[0,:])
#plot(t[:],Y12_2[1,:])
#xlabel('$time$')
#ylabel('$concentration$')
#title('LSRK12 at h = 0.01')
#savefig('Bru_Y12_2.png')

# LSRK 13
# create plot 
figure(2)
plot(t[:],Y13_0[0,:],label = "species X", color='red')
plot(t[:],Y13_0[1,:],label = "species Y", color='blue')
xlabel('$time$')
ylabel('$concentration$')
legend(loc="upper right")
title('LSRK13')
savefig('Bru_Y13.png')

## create plot 
#figure(5)
#plot(t[:],Y13_1[0,:])
#plot(t[:],Y13_1[1,:])
#xlabel('$time$')
#ylabel('$concentration$')
#title('LSRK13 at h = 0.05')
#savefig('Bru_Y13_1.png')
#
## create plot 
#figure(6)
#plot(t[:],Y13_2[0,:])
#plot(t[:],Y13_2[1,:])
#xlabel('$time$')
#ylabel('$concentration$')
#title('LSRK13 at h = 0.01')
#savefig('Bru_Y13_2.png')


# LSRK 14
# create plot 
figure(3)
plot(t[:],Y14_0[0,:],label = "species X", color='red')
plot(t[:],Y14_0[1,:],label = "species Y", color='blue')
xlabel('$time$')
ylabel('$concentration$')
legend(loc="upper right")
title('LSRK14')
savefig('Bru_Y14.png')

# create plot 
#igure(8)
#lot(t[:],Y14_1[0,:])
#lot(t[:],Y14_1[1,:])
#label('$time$')
#label('$concentration$')
#itle('LSRK14 at h = 0.05')
#avefig('Bru_Y14_1.png')

# create plot 
#igure(9)
#lot(t[:],Y14_2[0,:])
#lot(t[:],Y14_2[1,:])
#label('$time$')
#label('$concentration$')
#itle('LSRK14 at h = 0.01')
#avefig('Bru_Y14_2.png')

figure(4)
plot(t[:],Yerk4_0[0,:],label = "species X", color='red')
plot(t[:],Yerk4_0[1,:],label = "species Y", color='blue')
xlabel('$time$')
ylabel('$concentration$')
legend(loc="upper right")
title('ERK4')
savefig('Bru_Yerk4.png')


figure(5)
plot(t[:],Yfe_0[0,:],label = "species X", color='red')
plot(t[:],Yfe_0[1,:],label = "species Y", color='blue')
xlabel('$time$')
ylabel('$concentration$')
legend(loc="upper right")
title('FE')
savefig('Bru_Yfe.png')


figure(6)
plot(Y12_0[0,:],Y12_0[1,:],label = "LSRK(12,4)", color='green')
plot(Y13_0[0,:],Y13_0[1,:],label = "LSRK(13,4)", color='blue')
plot(Y14_0[0,:],Y14_0[1,:],label = "LSRK(14,4)", color='red')
plot(Yerk4_0[0,:],Yerk4_0[1,:],label = "ERK4", color='c')
plot(Yfe_0[0,:],Yfe_0[1,:],label = "FE", color='k')
xlabel('$concentration X$')
ylabel('$concentration Y$')
legend(loc="upper right")
title('Phase space plot')
savefig('Bru_phase.png')



# output observations
#print("\n  The first fixed step size result that is at all qualitatively correct\n  is the 10000-step run; however the orbits are 'tilted' a bit, and\n  these don't properly line up until the 20000-step run.  In the  adaptive,\n  10000-step and 20000-step runs, the plot is a bit 'kinky' at the right\n  end, but this is due to the coarse set of output times that we store.\n  It's also shocking that the adaptive run obtains the 'correct' answer\n  using only 402 steps (334 successful, 68 failed), which is a tiny\n  fraction of the work required with the fixed-step runs.")

# display all plots; these can be interacted with using the mouse
show()
