# -*- coding: utf-8 -*-
from __future__ import unicode_literals 
#script to generate data files needed for plotting
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats 
from scipy.interpolate import UnivariateSpline
import seaborn as sns 
import scikits.statsmodels.tools.tools as sm
import math

flow_rates = ['0', '1', '05', '2','4']
flow_rates_float = [0, 0.5, 1, 2, 4]
flow_rate_angle_files = ["0.dat", "05.dat", "1.dat", "2.dat", "4.dat"]


angle_bins = [i for i in range(0,90,5)]
colors = ['black', 'red', 'blue', 'green', 'yellow', 'magenta', 'orange', 'cyan', 'brown', 'purple', 'pink', 'gray', "tan"]

'''#this produces a bar plot, does not look very good! 
for i in range(len(flow_rates)):
	flow_rate = flow_rates[i]
	flow_rate_angle_file = flow_rate_angle_files[i]
	#angles = []
	with open(flow_rate_angle_file) as f:
		content = f.read().replace('\n', ' ').replace('\r', ' ')
		split_content = content.split(' ')
		print split_content
	print "printing content!!!!"
	#for x in content:
		#print x.strip()
	#for x in split_content:
	#	print x.strip()
	#	float(x.strip())
	angles = [float(x.strip()) for x in split_content if x != '']

	#generate bar plot/histogram
	print angles
	print angle_bins
	plt.hist(angles, bins = angle_bins, alpha = .5, color = colors[i], label = flow_rates[i], cumulative = True)


plt.xlabel('angle (degrees)')
plt.ylabel('counts')
plt.title('angle distributions')
plt.legend()
plt.savefig('angle_distributions.pdf')
plt.close()



#this will produce a smoothed probability distribution curve 
for i in range(len(flow_rates)):
	flow_rate = flow_rates[i]
	flow_rate_angle_file = flow_rate_angle_files[i]
	#angles = []
	with open(flow_rate_angle_file) as f:
		content = f.read().replace('\n', ' ').replace('\r', ' ')
		split_content = content.split(' ')
		print split_content
	print "printing content!!!!"
	#for x in content:
		#print x.strip()
	#for x in split_content:
	#	print x.strip()
	#	float(x.strip())
	angles = [float(x.strip()) for x in split_content if x != '']

	#generate bar plot/histogram
	print angles
	print angle_bins
	#density = gaussian_kde(expt_A_data_timepoint)
	#f = UnivariateSpline(lengths, expt_A_data_timepoint, s=.05)
	#plt.hist(angles, bins = angle_bins, alpha = .5, color = colors[i], label = flow_rates[i])
	sns.set_style('whitegrid')
	plot = sns.kdeplot(np.array(angles), bw=0.5, color = colors[i], label = flow_rates[i])
fig = plot.get_figure()
fig.savefig("angles_density.pdf")

'''

#also need to produce a CDF curve 
taps_mean_list = []
taps_sd_list = []
for i in range(len(flow_rates)):
	flow_rate = flow_rates[i]
	flow_rate_angle_file = flow_rate_angle_files[i]
	#angles = []
	with open(flow_rate_angle_file) as f:
		content = f.read().replace('\n', ' ').replace('\r', ' ')
		split_content = content.split(' ')
		print split_content
	print "printing content!!!!"


	taps = [float(x.strip()) for x in split_content if x != '']

	print taps
	n, min_max, mean, var, skew, kurt = scipy.stats.describe( taps )
	taps_mean_list.append(mean)
	taps_sd_list.append(math.sqrt(var))

flow_rates_glass_float = [ 0.5, 1, 2, 4]
taps_glass_mean_list = [21.28, 27.4, 27.38, 18.94]
taps_glass_sd_list = [3.45, 2.6, 2.71, 3.07]
	
	
plt.figure()
plt.errorbar(flow_rates_float, taps_mean_list, yerr = taps_sd_list, label = 'cell')
plt.errorbar(flow_rates_glass_float, taps_glass_mean_list, yerr = taps_glass_sd_list, label = 'glass')
plt.xlim(-.5, 5)
plt.xlabel('flow rate µL/s')
plt.ylabel('mean normalized total area projection (pixels)')
plt.title('normalized total area projection')
#plt.legend( loc=4 )
plt.savefig('tap_analysis.pdf')
plt.close()


