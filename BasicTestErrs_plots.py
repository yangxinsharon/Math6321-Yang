#!/usr/bin/env python
#
# Plotting script for basic test error comparison
#
# Sharon Yang
# Math 6321 @ SMU
# Fall 2020

# imports
from numpy import *
from matplotlib.pyplot import *
from prettytable import PrettyTable

# First, load our saved data files
h = loadtxt('basic_h.txt')
err12 = loadtxt('basic_err_12.txt')
err13 = loadtxt('basic_err_13.txt')
err14 = loadtxt('basic_err_14.txt')
err_erk4 = loadtxt('basic_err_erk4.txt')
err_fe = loadtxt('basic_err_fe.txt')
conv = loadtxt('basic_conv.txt')


# create plot 
figure(1)
plot(h[:],err12[:],label = "LSRK(12,4)", color='green',marker='^')
plot(h[:],err13[:],label = "LSRK(13,4)", color='blue',marker='o')
plot(h[:],err14[:],label = "LSRK(14,4)", color='red',marker='s')
plot(h[:],err_erk4[:],label = "ERK4", color='c',marker='+')
plot(h[:],err_fe[:],label = "Forward Euler", color='k',marker='.')
xlabel('$h$')
ylabel('$Total	Error$')
xscale('log')
yscale('log')
legend(loc="upper left")
title('Total error for the basic test system at t_end = 1.4')
savefig('BasicTestErrs.png')


#field_names = ('h','LSRK(12,4)','LSRK(13,4)','LSRK(14,4)','ERK4')
#table = PrettyTable(field_names=field_names)
table = PrettyTable()
table.add_column('h',conv[:,0])
table.add_column('LSRK(12,4)',conv[:,1])
table.add_column('LSRK(13,4)',conv[:,2])
table.add_column('LSRK(14,4)',conv[:,3])
table.add_column('ERK4',conv[:,4])
table.add_column('Forward Euler',conv[:,5])
print(table)

data = table.get_string()
with open('basic_table.txt', 'w') as file:
    file.write(data)

# display all plots; these can be interacted with using the mouse
show()
