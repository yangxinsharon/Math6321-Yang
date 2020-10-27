#!/usr/bin/env python
#
# Plotting script for homework 4, problem 3.
#
# Sharon Yang
# Math 6321 @ SMU
# Fall 2020

# imports
from numpy import *
from matplotlib.pyplot import *

# First, load our saved data files
Y0 = loadtxt('AdaptRKF.txt')
Y1 = loadtxt('ERK4_100.txt')
Y2 = loadtxt('ERK4_1000.txt')
Y3 = loadtxt('ERK4_10000.txt')
Y4 = loadtxt('ERK4_20000.txt')

# create plots for Runge-Kutta-Fehlberg
figure(1,figsize=(10,4))
subplot(1,2,1)
plot(Y0[0,:],Y0[2,:])
xlabel('$u1$')
ylabel('$u2$')
title('AdaptRKF')
savefig('figure1.png')

# create plots for steps = 100
figure(2,figsize=(10,4))
subplot(1,2,1)
plot(Y1[0,:],Y1[2,:])
xlabel('$u1$')
ylabel('$u2$')
title('ERK4 with 100 steps')
savefig('figure2.png')

# create plots for steps = 1000
figure(3,figsize=(10,4))
subplot(1,2,1)
plot(Y2[0,:],Y2[2,:])
xlabel('$u1$')
ylabel('$u2$')
title('ERK4 with 1000 steps')
savefig('figure3.png')

# create plots for steps = 10000
figure(4,figsize=(10,4))
subplot(1,2,1)
plot(Y3[0,:],Y3[2,:])
xlabel('$u1$')
ylabel('$u2$')
title('ERK4 with 10000 steps')
savefig('figure4.png')

# create plots for steps = 20000
figure(5,figsize=(10,4))
subplot(1,2,1)
plot(Y4[0,:],Y4[2,:])
xlabel('$u1$')
ylabel('$u2$')
title('ERK4 with 20000 steps')
savefig('figure5.png')

print("Using the classical Runge-Kutta method, it shows a spiral and a limit cycle for 100 and 1000 steps. Therefore, 10000 steps are needed. ")
# display all plots; these can be interacted with using the mouse
show()

