# -*- coding: utf-8 -*-
from __future__ import unicode_literals 
#script to generate data files needed for plotting
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.interpolate import UnivariateSpline
import seaborn as sns 
import scikits.statsmodels.tools.tools as sm
import scipy.stats 

flow_rates_glass = ['0', '1', '2', '4','8','12','16','20','24','32']
flow_rate_angle_files_glass = ["0.dat", "1.dat", "2.dat", "4.dat", "8.dat", "12.dat", "16.dat", "20.dat", "24.dat", "32.dat"]

flow_rates_cell = ['0', '1', '2', '4', '6', '8', '10', '12', '16' ,'20' ]
flow_rate_angle_files_cell = ["0.dat", "1.dat", "2.dat", "4.dat", "6.dat", "8.dat", "10.dat", "12.dat", "16.dat", "20.dat"]


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

#gathering angle data for glass 
flow_rates_float_glass = []
variance_list = []
for i in range(len(flow_rates_glass)):
	flow_rate_glass = flow_rates_glass[i]
	flow_rate_angle_file_glass = flow_rate_angle_files_glass[i]
	flow_rates_float_glass.append(float(flow_rate_glass))

	#angles = []
	with open("glass_surface_data/glass_data_3_27/"+flow_rate_angle_file_glass) as f:
		content = f.read().replace('\n', ' ').replace('\r', ' ')
		split_content = content.split(' ')
		print split_content
	print "printing content!!!!"

	#for x in content:
		#print x.strip()
	#for x in split_content:
	#	print x.strip()
	#	float(x.strip())
	angles_glass = [float(x.strip()) for x in split_content if x != '']

	#generate bar plot/histogram
	print angles_glass

	n, min_max, mean, var, skew, kurt = scipy.stats.describe( angles_glass )

	print "variance is: ", var 

	variance_list.append(var)

plt.plot(flow_rates_float_glass, variance_list, color = "blue", label = "glass")

#gathering angle data for cell 
flow_rates_float_cell = []
variance_list = []
for i in range(len(flow_rates_cell)):
	flow_rate_cell = flow_rates_cell[i]
	flow_rate_angle_file_cell = flow_rate_angle_files_cell[i]
	flow_rates_float_cell.append(float(flow_rate_cell))

	#angles = []
	with open("cell_surface_data/sis_data_3_27_cell_data/"+flow_rate_angle_file_cell) as f:
		content = f.read().replace('\n', ' ').replace('\r', ' ')
		split_content = content.split(' ')
		print split_content
	print "printing content!!!!"

	#for x in content:
		#print x.strip()
	#for x in split_content:
	#	print x.strip()
	#	float(x.strip())
	angles_cell = [float(x.strip()) for x in split_content if x != '']

	#generate bar plot/histogram
	print angles_cell

	n, min_max, mean, var, skew, kurt = scipy.stats.describe( angles_cell )

	print "variance is: ", var 

	variance_list.append(var)

plt.plot(flow_rates_float_glass, variance_list, color = "red", label = "cell")
plt.xlabel('flow rate (µL/s)')
plt.ylabel('variance (degrees squared)')
plt.legend()




plt.savefig('angle_variance_vs_flow_rate.pdf')
	#density = gaussian_kde(expt_A_data_timepoint)
	#f = UnivariateSpline(lengths, expt_A_data_timepoint, s=.05)
	#plt.hist(angles, bins = angle_bins, alpha = .5, color = colors[i], label = flow_rates[i])
	#sns.set_style('whitegrid')
	#plot = sns.kdeplot(np.array(angles), bw=0.5, color = colors[i], label = flow_rates[i])
'''plt.step(ecdf_data.x, ecdf_data.y, color = colors[i], label = flow_rates[i]+" µL/s")
plt.xlabel('angle (degrees)')
plt.ylabel('ecdf (step function)')
plt.title('angle distribution ecdf')
plt.legend( loc=4 )
plt.savefig('msad_angle_vs_flow_rate.pdf')
plt.close()'''


