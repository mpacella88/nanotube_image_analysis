import numpy as np
import matplotlib.pyplot as plt 
import scipy.stats

timepoints = ['1', '2', '3', '4']
timepoints_float = [1, 2, 3, 4]
timepoint_data_files = ["1.dat", "2.dat", "3.dat", "4.dat"]


#angle_bins = [i for i in range(0,90,5)]
#colors = ['black', 'red', 'blue', 'green', 'yellow', 'magenta', 'orange', 'cyan', 'brown', 'purple', 'pink', 'gray', "tan"]

for i in range(len(timepoints)):
	timepoint = timepoints[i]
	timpoint_data_file = timepoint_data_files[i]
	with open(timpoint_data_file) as f:
		content = f.read().replace('\n', ' ')
	split_network_sizes = content.split(' ')
	network_sizes = [float(x.strip()) for x in split_network_sizes if x != '']

	print network_sizes
	network_sizes_weighted_by_seed = []
	for network_size in network_sizes:
		for seed in range(int(network_size)):
			network_sizes_weighted_by_seed.append(network_size)
	print network_sizes_weighted_by_seed

	#make a histogram of the weighted network sizes
	plt.hist(network_sizes_weighted_by_seed)
	plt.xlabel("network size")
	plt.ylabel("count (weighted by number of seeds in network)")
	plt.savefig("network_size_weighted_by_seed_"+timepoint+".pdf")

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
plt.close()'''