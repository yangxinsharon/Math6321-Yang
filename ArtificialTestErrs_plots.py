#!/usr/bin/env python
#
# Plotting script for artificial test error comparison
#
# Sharon Yang
# Math 6321 @ SMU
# Fall 2020

# imports
from numpy import *
from matplotlib.pyplot import *
from prettytable import PrettyTable


# First, load our saved data files
h = loadtxt('artificial_h.txt')
err12 = loadtxt('artificial_err_12.txt')
err13 = loadtxt('artificial_err_13.txt')
err14 = loadtxt('artificial_err_14.txt')
err_erk4 = loadtxt('artificial_err_erk4.txt')
err_fe = loadtxt('artificial_err_fe.txt')

conv12 = loadtxt('artificial_conv_12.txt')
conv13 = loadtxt('artificial_conv_13.txt')
conv14 = loadtxt('artificial_conv_14.txt')
conv_erk4 = loadtxt('artificial_conv_erk4.txt')
conv_fe= loadtxt('artificial_conv_fe.txt')

# create table
# lambda = 100
table = PrettyTable()
table.title = 'Artificial Test Errors for lambda = 100'
table.add_column('h',h[:])
table.add_column('LSRK(12,4)',err12[:,0])
table.add_column('LSRK(13,4)',err13[:,0])
table.add_column('LSRK(14,4)',err14[:,0])
table.add_column('ERK4',err_erk4[:,0])
table.add_column('Forward Euler',err_fe[:,0])

print(table)

data100 = table.get_string()
with open('artificial_table_100.txt', 'w') as file:
    file.write(data100)


# lambda = 200
table = PrettyTable()
table.title = 'Artificial Test Errors for lambda = 200'
table.add_column('h',h[:])
table.add_column('LSRK(12,4)',err12[:,1])
table.add_column('LSRK(13,4)',err13[:,1])
table.add_column('LSRK(14,4)',err14[:,1])
table.add_column('ERK4',err_erk4[:,1])
table.add_column('Forward Euler',err_fe[:,1])

print(table)

data200 = table.get_string()
with open('artificial_table_100.txt', 'w') as file:
    file.write(data200)



# lambda = 400
table = PrettyTable()
table.title = 'Artificial Test Errors for lambda = 400'
table.add_column('h',h[:])
table.add_column('LSRK(12,4)',err12[:,2])
table.add_column('LSRK(13,4)',err13[:,2])
table.add_column('LSRK(14,4)',err14[:,2])
table.add_column('ERK4',err_erk4[:,2])
table.add_column('Forward Euler',err_fe[:,2])

print(table)

data300 = table.get_string()
with open('artificial_table_100.txt', 'w') as file:
    file.write(data300)
# display all plots; these can be interacted with using the mouse
show()
