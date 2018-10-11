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

#flow_rates = ['0', '1', '05', '2','4']
#flow_rates_float = [0, 0.5, 1, 2, 4]
#flow_rate_angle_files = ["0.dat", "05.dat", "1.dat", "2.dat", "4.dat"]


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



shear_stress_glass_float = [ .05283, .10566, .21132, .42264, .63396, .84528, 1.0566, 1.26792, 1.69056, 2.1132]
taps_glass_mean_list = [52.29, 36.012, 20.57, 14.96, 11.06, 8.238, 5.24, 3.832, 2.21, 1.82]
taps_glass_sd_list = [2.64, 1.68, 1.35, 1.02, .667, .68, .626, .474, .452, .647]

shear_stress_cell_float = [ .05283, .10566, .21132, .42264, .63396, .84528, 1.0566, 1.26792]
taps_cell_mean_list = [54.52, 49.75, 32.88, 26.45, 24.08, 17.783, 9.625, 9.54 ]
taps_cell_sd_list = [3.53, 1.89, 2.29, 1.52, 1.46, .997, 1.16, .744 ]

shear_stress_model_float = [.21132, .42264, .63396, .84528, 1.0566, 1.26792, 1.69056, 2.1132]
taps_model_mean_list = [79.47, 58.421, 41.38, 29.36, 20.38, 14.57, 10.432, 7.542]
	
plt.figure()
plt.errorbar(shear_stress_cell_float, taps_cell_mean_list, yerr = taps_cell_sd_list, label = 'cell')
plt.errorbar(shear_stress_glass_float, taps_glass_mean_list, yerr = taps_glass_sd_list, label = 'glass')
plt.plot(shear_stress_model_float, taps_model_mean_list, label = 'model')
plt.xlim(0, 2.5)
plt.xlabel('shear stress (dyn*s/cm^2)')
plt.ylabel('mean total angle')
plt.title('mean angle visual inspection')
plt.legend( loc=1 )
plt.savefig('sisi_z_projection_visual_analysis.pdf')
plt.close()


