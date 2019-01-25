from math import sqrt
from skimage import data
from skimage.feature import blob_dog, blob_log, blob_doh
from skimage.color import rgb2gray
from skimage import io 
import sys

import matplotlib.pyplot as plt
from skimage.data import camera
from skimage.filters import roberts, sobel, scharr, prewitt

def count_seeds(filename, thresh = .01):
	"""count the number of seeds in a grayscale image using the laplacian
	of gaussians method""" 
	#test

	image_gray = io.imread(filename)
	image = io.imread(filename)
	#image_gray = rgb2gray(image)

	blobs_log = blob_log(image_gray, max_sigma=30, num_sigma=10, threshold=thresh)
	#blobs_log = blob_log(image_gray, max_sigma=1, num_sigma=10, threshold=.1)

	print "number of seeds: "+str(len(blobs_log))

	# Compute radii in the 3rd column.
	blobs_log[:, 2] = blobs_log[:, 2] * sqrt(2)


	blobs_list = [blobs_log]
	colors = ['yellow']
	titles = ['Laplacian of Gaussian']
	sequence = zip(blobs_list, colors, titles)

	fig, ax = plt.subplots(figsize=(9, 3))
	ax.set_aspect('equal')
	#ax = axes.ravel()

	for idx, (blobs, color, title) in enumerate(sequence):
	    #ax[idx].set_title(title)
	    ax.imshow(image, interpolation='nearest')
	    for blob in blobs:
	        y, x, r = blob
	        c = plt.Circle((x, y), r, color=color, linewidth=2, fill=False)
	        ax.add_patch(c)
	    #ax[idx].set_axis_off()

	plt.tight_layout()
	plt.show()

	print blobs_log 

def calc_distance(endpoint1, endpoint2):
	#simple distance calculation
	distance_squared = (endpoint1[0]-endpoint2[0]) * (endpoint1[0]-endpoint2[0]) + (endpoint1[1]-endpoint2[1]) * (endpoint1[1]-endpoint2[1])
	distance = sqrt(distance_squared)

	return distance

def detect_capping_percentage(cap_filename, seed_filename):
	#simple method to detect the percentage of seeds that have a cap in close proximity

	#seed and cap must be closer than this distance to be counted as attached
	cutoff_distance = 2

	#read in seed image
	image_seed = io.imread(seed_filename)

	#read in caps image 
	image_cap = io.imread(cap_filename)

	seed_blobs = blob_log(image_seed, max_sigma=30, num_sigma=10, threshold=.01)
	cap_blobs = blob_log(image_cap, max_sigma=30, num_sigma=10, threshold=.01)

	seed_coordinates = seed_blobs[:,0:2]
	cap_coordinates = cap_blobs[:,0:2]

	#loop over each seed_coordinates and compute distance to nearest cap coordinate
	#count as attached if distance is less than cutoff distance

	attached_seed_coordinates = []

	for seed in seed_coordinates:
		for cap in cap_coordinates:
			distance = calc_distance(seed, cap)
			if distance <= cutoff_distance:
				print "attached cap detected"
				attached_seed_coordinates.append(seed)
				break

	attached_seed_percentage = float(len(attached_seed_coordinates))/float(len(seed_coordinates))

	print "attached seed percentage is: ", attached_seed_percentage
	print "total seeds: ", len(seed_coordinates)
	print "total seeds with caps is: ", len(attached_seed_coordinates)
	print "total caps: ", len(cap_coordinates)



	

	



