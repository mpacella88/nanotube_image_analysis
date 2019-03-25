from math import sqrt
import sys
import numpy as np 
import matplotlib.pyplot as plt 
import math
import os
import glob 
from matplotlib import cm 


#quick script to create sub-plots of histograms for each timepoint network growth time series data 
data_2 = open("1_network_seeds_cluster_size.dat")
data_2_network_sizes=data_2.read().split('\n')
data_2.close()
data_2 = []
for network in data_2_network_sizes:
	if network == '':
		continue 
	data_2.append(float(network))
print data_2

data_5 = open("6_network_seeds_cluster_size.dat")
data_5_network_sizes=data_5.read().split('\n')
data_5.close()
data_5 = []
for network in data_5_network_sizes:
	if network == '':
		continue 
	data_5.append(float(network))
print data_5

data_12 = open("25_network_seeds_cluster_size.dat")
data_12_network_sizes=data_12.read().split('\n')
data_12.close()
data_12 = []
for network in data_12_network_sizes:
	if network == '':
		continue 
	data_12.append(float(network))
print data_12

data_24 = open("49_network_seeds_cluster_size.dat")
data_24_network_sizes=data_24.read().split('\n')
data_24.close()
data_24 = []
for network in data_24_network_sizes:
	if network == '':
		continue 
	data_24.append(float(network))
print data_24

'''data_48 = open("23_network_seeds_cluster_size.dat")
data_48_network_sizes=data_48.read().split('\n')
data_48.close()
data_48 = []
for network in data_48_network_sizes:
	if network == '':
		continue 
	data_48.append(float(network))
print data_48'''

'''data_72 = open("72_network_seeds_cluster_size.dat")
data_72_network_sizes=data_72.read().split('\n')
data_72.close()
data_72 = []
for network in data_72_network_sizes:
	if network == '':
		continue 
	data_72.append(float(network))
print data_72'''

fig, axes = plt.subplots(nrows=3, ncols=2)
ax0, ax1, ax2, ax3, ax4, ax5 = axes.flatten()

#ax0.hist(data_2, bins = range(0,15000,1000))
ax0.hist(data_2)
ax0.set_title('1hrs')
ax0.x_scale('log')
ax0.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#ax1.hist(data_5, bins = range(0,15000,1000))
ax1.hist(data_5)
ax1.set_title('2hrs')
ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#ax2.hist(data_12, bins = range(0,15000,1000))
ax2.hist(data_12)
ax2.set_title('4hrs')
ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#ax3.hist(data_24, bins = range(0,15000,1000))
ax3.hist(data_24)
ax3.set_title('8hrs')
ax3.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#ax4.hist(data_48, bins = range(0,15000,1000))
'''ax4.hist(data_48)
ax4.set_title('23hrs')
ax4.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax4.set_xlabel("cluster size (sq pixels)")'''
'''ax5.hist(data_72, bins = range(0,15000,1000))
ax5.set_title('72hrs')
ax5.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax5.set_xlabel("cluster size (sq pixels)")'''

fig.tight_layout()
plt.savefig("hist_subplots_cluster_different_scales.pdf")











