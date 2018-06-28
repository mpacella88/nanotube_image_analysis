# -*- coding: utf-8 -*-
from __future__ import unicode_literals 
#script to generate data files needed for plotting
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.interpolate import UnivariateSpline
import seaborn as sns 
import scikits.statsmodels.tools.tools as sm

flow_rates = ['0', '2', '4', '8', '12', '16', '20', '24', '32', ]
flow_rate_angle_files = ["0.dat", "2.dat", "4.dat", "8.dat", "12.dat", "16.dat", "20.dat", "24.dat", '32.dat']

angle_bins = [i for i in range(0,90,5)]
colors = ['black', 'red', 'blue', 'green', 'yellow', 'magenta', 'orange', 'cyan', 'brown', 'purple', 'pink']

#this produces a bar plot, does not look very good! 
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



'''#this will produce a smoothed probability distribution curve 
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
fig.savefig("angles_density.pdf")'''

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
	angles_random_sign = []
	print angles
	j = 1
	for angle in angles:
		if j%2 == 0:
			angles_random_sign.append(-1.0 * angle)
		else:
			angles_random_sign.append(1.0 * angle)

		j+=1


	#generate bar plot/histogram
	print angles_random_sign
	print angle_bins
	#density = gaussian_kde(expt_A_data_timepoint)
	#f = UnivariateSpline(lengths, expt_A_data_timepoint, s=.05)
	#plt.hist(angles, bins = angle_bins, alpha = .5, color = colors[i], label = flow_rates[i])
	sns.set_style('whitegrid')
	plot = sns.kdeplot(np.array(angles_random_sign), bw=0.5, color = colors[i], label = flow_rates[i] + ' µ  L/s')
	#sns.xlabel("angle (degrees)")
	#sns.xlim(0.0, 90.0)
	plot.set(xlim = (0.0, 90.0))
	plot.set(xlabel = "angles (degrees)")
	plot.set(ylabel = "density")
fig = plot.get_figure()
#fig.xlabel("angle (degrees)")
#fig.xlim(0.0, 90.0)
fig.savefig("angles_half_density.pdf")



#also need to produce a CDF curve 
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
	ecdf_data = sm.ECDF(angles)
	print 'ecdf is: ', ecdf_data
	#density = gaussian_kde(expt_A_data_timepoint)
	#f = UnivariateSpline(lengths, expt_A_data_timepoint, s=.05)
	#plt.hist(angles, bins = angle_bins, alpha = .5, color = colors[i], label = flow_rates[i])
	#sns.set_style('whitegrid')
	#plot = sns.kdeplot(np.array(angles), bw=0.5, color = colors[i], label = flow_rates[i])
	plt.step(ecdf_data.x, ecdf_data.y, color = colors[i])
plt.xlabel('angle (degrees)')
plt.ylabel('ecdf (step function)')
plt.title('angle distribution ecdf')
plt.legend()
plt.savefig('angle_distributions_ecdf.pdf')
plt.close()


