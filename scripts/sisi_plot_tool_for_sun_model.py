# -*- coding: utf-8 -*-
from __future__ import unicode_literals 
#script to generate data files needed for plotting
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.interpolate import UnivariateSpline
import seaborn as sns 
import scikits.statsmodels.tools.tools as sm
import csv

flow_rates = ['0', '2', '4', '6', '8', '10', '12', '16', '20' ]
flow_rate_angle_files = ["0.dat", "2.dat", "4.dat", "6.dat", "8.dat", "10.dat", "12.dat", "16.dat", "20.dat"]

#flow_rates = ['0', '2', '4', '8', '12', '16', '20' ]
#flow_rate_angle_files = ["0.dat", "2.dat", "4.dat", "8.dat", "12.dat", "16.dat", "20.dat"]

angle_bins = [i for i in range(0,90,5)]
colors = ['black', 'red', 'blue', 'green', 'yellow', 'magenta', 'orange', 'cyan', 'brown', 'purple', 'pink', 'gray', "tan"]
'''
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

'''

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
'''
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
'''

for flow_rate in flow_rates:
	if flow_rate == '14' or flow_rate == '18':
		continue 
	r = csv.DictReader(open('sun_model.csv', "rb"), dialect = 'excel')
	angles_flow = []
	for row in r:
		if float(row['phi']) <=90.0:
			print "adding data point to cumulative sum: ", row['phi']
			angles_flow.append(float(row[flow_rate]))
	angles = np.arange(1.5, 91.5, 3.0)
	cdf = np.cumsum(angles_flow, dtype = 'float')
	cdf_sum = sum(angles_flow)
	cdf = cdf/cdf_sum 

	plt.plot(angles, cdf, label = flow_rate+ ' µL/s', color = colors[flow_rates.index(flow_rate)])

'''#also need to produce a CDF curve 
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
	'''
plt.xlabel('angle (degrees)')
plt.ylim([0.0,1.0])
plt.ylabel('ecdf (step function)')
plt.title('model angle distribution ecdf')
plt.legend(loc = 4)
plt.savefig('model_angle_distributions_ecdf_matching_cell.pdf')
plt.close()


