#!/usr/bin/env python
#
# Plotting script for homework 5, problem 4.
#
# Sharon Yang
# Math 6321 @ SMU
# Fall 2020

# imports
from numpy import *
from matplotlib.pyplot import *
from prettytable import PrettyTable

# First, load our saved data files
t = loadtxt('tspan.txt')
H_devi_ERK4 = loadtxt('H_devi_ERK4.txt')
H_devi_RIA3 = loadtxt('H_devi_RIA3.txt')
H_devi_GL2 = loadtxt('H_devi_GL2.txt')
H_devi_RKF = loadtxt('H_devi_RKF.txt')

# create plots
figure(1)
plot(t,H_devi_ERK4[:,0])
xlabel('$t$')
ylabel('$H(t)-H(0)$')
title('Energy error using ERK4 method with h=2.0')
savefig('ERK4_2.0.png')

figure(2)
plot(t,H_devi_ERK4[:,1])
xlabel('$t$')
ylabel('$H(t)-H(0)$')
title('Energy error using ERK4 method with h=0.2')
savefig('ERK4_0.2.png')

figure(3)
plot(t,H_devi_RIA3[:,0])

xlabel('$t$')
ylabel('$H(t)-H(0)$')
title('Energy error using RIA3 method with h=2.0')
savefig('RIA3_2.0.png')

figure(4)
plot(t,H_devi_RIA3[:,1])
xlabel('$t$')
ylabel('$H(t)-H(0)$')
title('Energy error using RIA3 method with h=0.2')
savefig('RIA3_0.2.png')

figure(5)
plot(t,H_devi_GL2[:,0])
xlabel('$t$')
ylabel('$H(t)-H(0)$')
title('Energy error using GL2 method with h=2.0')
savefig('GL2_2.0.png')

figure(6)
plot(t,H_devi_GL2[:,1])
xlabel('$t$')
ylabel('$H(t)-H(0)$')
title('Energy error using GL2 method with h=0.2')
savefig('GL2_0.2.png')

figure(7)
plot(t,H_devi_RKF[:,0])
xlabel('$t$')
ylabel('$H(t)-H(0)$')
title('Energy error using RKF method with h=2.0')
savefig('RKF_2.0.png')

figure(8)
plot(t,H_devi_RKF[:,1])
xlabel('$t$')
ylabel('$H(t)-H(0)$')
title('Energy error using RKF method with h=0.2')
savefig('RKF_0.2.png')

a=np.max(abs(H_devi_ERK4),axis=0)
b=np.max(abs(H_devi_RIA3),axis=0)
c=np.max(abs(H_devi_GL2),axis=0)
d=np.max(abs(H_devi_RKF),axis=0)
field_names = ('h','ERK4','RIA3','GL2','rtol','RKF')
table = PrettyTable(field_names=field_names)
table.add_row(['2.0',a[0],b[0],c[0],'e-4',d[0]])
table.add_row(['0.2',a[1],b[1],c[1],'e-6',d[1]])
print(table)

# output observations
print("\n GL2 is better than RIA3, and RIA3 is better than ERK4.\n The adaptive RKF meets the relative tolerances giving good approximation compared with other methods when h =2.")
# display all plots; these can be interacted with using the mouse
#show()
