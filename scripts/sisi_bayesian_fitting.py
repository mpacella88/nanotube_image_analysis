# -*- coding: utf-8 -*-
from __future__ import unicode_literals 
#script to generate data files needed for plotting
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.interpolate import UnivariateSpline
import seaborn as sns 
import scikits.statsmodels.tools.tools as sm
import pymc3 as pm
import os

flow_rates = ['test']
flow_rate_angle_files = ['120 1_cell1_5NT.dat']

angle_bins = [i for i in range(0,90,5)]
colors = ['black', 'red', 'blue', 'green', 'yellow', 'magenta', 'orange', 'cyan', 'brown', 'purple', 'pink']

#using pymc3 and bayesian fitting to determine distribution for shear stress
#prior will be uniform in some range determined by the true shear stress
#data will be the list of measured angles from sisi
#likelihood fxn will be sean sun's model with all parameters other than shear stress determined

#getting our data (flow rates and corresponding angle files defined at the beginning
#of this script)
flow_rate_angle_file_list = os.listdir('dat_files')
for flow_rate_angle_file in flow_rate_angle_file_list:
	#flow_rate = flow_rates[i]
	#flow_rate_angle_file = flow_rate_angle_files[i]
	#angles = []
	with open(flow_rate_angle_file) as f:
		content = f.read().replace('\n', ' ').replace('\r', ' ')
		split_content = content.split(' ')
		print split_content

	angles = [float(x.strip()) for x in split_content if x != '']
	model = pm.Model()

	with model:
	    group1_std = pm.Uniform('group1_std', lower = 0, upper = 90)

	with model:
	    group1 = pm.HalfNormal('group1', sd = group1_std, observed = angles)

	with model:
	    trace = pm.sample(2500, njobs=1)   

	#pm.traceplot(trace)
	pm.plot_posterior(trace, varnames=['group1_std'])

	#plt.show()
	plt.savefig(flow_rate_angle_file + '_posterior.png')
#likelihood fxn
# theta_i ~ seans_model(shear stress)
# for now theta_i ~ half_normal(variance)



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
	print 'ecdf x is: ', ecdf_data.x
	print 'ecdf y is: ', ecdf_data.y 
	#density = gaussian_kde(expt_A_data_timepoint)
	#f = UnivariateSpline(lengths, expt_A_data_timepoint, s=.05)
	#plt.hist(angles, bins = angle_bins, alpha = .5, color = colors[i], label = flow_rates[i])
	#sns.set_style('whitegrid')
	#plot = sns.kdeplot(np.array(angles), bw=0.5, color = colors[i], label = flow_rates[i])
	plt.step(ecdf_data.x, ecdf_data.y, color = colors[i], label = flow_rates[i]+" µL/s")
plt.xlabel('angle (degrees)')
plt.ylabel('ecdf (step function)')
plt.title('angle distribution ecdf')
plt.legend( loc=4 )
plt.savefig('angle_distributions_ecdf.pdf')
plt.close()

'''
