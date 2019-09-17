from math import sqrt
from skimage import data
from skimage.feature import blob_dog, blob_log, blob_doh
from skimage.color import rgb2gray
from skimage import io 
import sys
import numpy as np 
import matplotlib.pyplot as plt
from skimage.data import camera
from skimage.filters import roberts, sobel, scharr, prewitt
from scipy.spatial.distance import euclidean
from scipy.spatial import Delaunay, ConvexHull
import networkx as nx
import math
import os
import glob 
from scipy.interpolate import splprep, splev
from scipy import ndimage as ndi
from skimage.feature import canny
from skimage.morphology import closing, disk
from skimage import img_as_uint, img_as_bool, img_as_float
from skimage.draw import circle 
from skimage.measure import label, regionprops, find_contours
from scipy import ndimage
from matplotlib import cm 

#the way this will work, all abs will be included for calculating the binding fraction with real seeds, 
#for finding the simulated "null" binding fraction based on the density of abs on the membrane I will
#use sisi's optimal membrane thickness to first cull the number of abs counted as being on the membrane 

#only count seeds that are within proximity of membrane 

def count_seeds(filename, blob_log_thresh = .01):
	"""count the number of seeds in a grayscale image using the laplacian
	of gaussians method"""

	image_gray = io.imread(filename)
	image = io.imread(filename)
	#image_gray = rgb2gray(image)

	blobs_log = blob_log(image_gray, max_sigma=30, num_sigma=10, threshold= blob_log_thresh)
	#blobs_log = blob_log(image_gray, max_sigma=1, num_sigma=10, threshold=.1)

	print "number of seeds: "+str(len(blobs_log))

	# Compute radii in the 3rd column.
	blobs_log[:, 2] = blobs_log[:, 2] * sqrt(2)


	blobs_list = [blobs_log]
	colors = ['yellow']
	titles = ['Laplacian of Gaussian']
	sequence = zip(blobs_list, colors, titles)

	fig, ax = plt.subplots(figsize=(3, 3))
	ax.set_aspect('equal')
	#ax = axes.ravel()

	for idx, (blobs, color, title) in enumerate(sequence):
	    #ax[idx].set_title(title)
	    ax.imshow(image, interpolation='nearest')
	    for blob in blobs:
	        y, x, r = blob
	        if r <= 2.0:
	        	c = plt.Circle((x, y), r, color=color, linewidth=2, fill=False)
	        	ax.add_patch(c)
	    #ax[idx].set_axis_off()

	#plt.tight_layout()
	plt.savefig("plots/"+filename+"_blobs.pdf")
	plt.clf()

	print blobs_log 

def calc_distance(endpoint1, endpoint2):
	#simple distance calculation
	distance_squared = (endpoint1[0]-endpoint2[0]) * (endpoint1[0]-endpoint2[0]) + (endpoint1[1]-endpoint2[1]) * (endpoint1[1]-endpoint2[1])
	distance = sqrt(distance_squared)

	return distance



def detect_ab_seed_binding_fraction(ab_filename, seed_filename):
	#simple method to detect the percentage of seeds bound to antibodies presented on a cell membrane 

	#seed and cap must be closer than this distance to be counted as attached
	cutoff_distance = 5

	#read in seed image
	image_seed = io.imread(seed_filename)

	#read in caps image 
	image_ab = io.imread(ab_filename)

	#count_seeds(seed_filename, .005)
	#count_seeds(ab_filename, .005)

	seed_blobs = blob_log(image_seed, max_sigma=30, num_sigma=10, threshold=.005)
	ab_blobs = blob_log(image_ab, max_sigma=30, num_sigma=10, threshold=.005)

	#seed_coordinates = seed_blobs[:,0:2]
	#ab_coordinates = cap_blobs[:,0:2]

	#loop over each seed_coordinates and compute distance to nearest ab coordinate
	#count as attached if distance is less than cutoff distance

	attached_seed_coordinates = []
	ab_coords_list = []
	seed_coords_list = []

	for seed in seed_blobs:
		if seed[2] <= 2.0:#only consider blobs with small radii, these are the seeds/abs
			seed_coordinates = seed[0:2]
			seed_coords_list.append(seed_coordinates)

	for ab in ab_blobs:
		if ab[2] <= 2.0:#only consider blobs with small radii, these are the seeds/abs
			ab_coordinates = ab[0:2]
			ab_coords_list.append(ab_coordinates)



	for seed_coords in seed_coords_list:
		for ab_coords in ab_coords_list:
			distance = calc_distance(seed_coords, ab_coords)
			if distance <= cutoff_distance:
				print "attached seed detected"
				attached_seed_coordinates.append(seed_coords)
				break
	if len(seed_coords_list) == 0:
		attached_seed_percentage = 0.0
	else:
		attached_seed_percentage = float(len(attached_seed_coordinates))/float(len(seed_coords_list))

	print seed_coords_list
	print "attached seed percentage is: ", attached_seed_percentage
	print "total seeds: ", len(seed_coords_list)
	print "total seeds with abs is: ", len(attached_seed_coordinates)
	print "total abs: ", len(ab_coords_list)

	return attached_seed_percentage
def report_binding_fraction(seed_coords, ab_coords, cutoff_distance = 5.0):
	#simple method for calculating the binding fraction given an array of seed coordinates and an array
	#of antibody coordinates
	attached_seed_coordinates = []
	for seed_coord in seed_coords:
		for ab_coord in ab_coords:
			distance = calc_distance(seed_coord, ab_coord)
			#print "distance is: ", distance 
			if distance <= cutoff_distance:
				#print "attached seed detected"
				#print "distance is: ", distance 
				attached_seed_coordinates.append(seed_coord)
				break
	if len(seed_coords) == 0:
		attached_seed_percentage = 0.0
	else:
		attached_seed_percentage = float(len(attached_seed_coordinates))/float(len(seed_coords))

	#print seed_coords
	#print "attached seed percentage is: ", attached_seed_percentage
	#print "total seeds: ", len(seed_coords)
	#print "total seeds with abs is: ", len(attached_seed_coordinates)
	#print "total abs: ", len(ab_coords)

	return attached_seed_percentage

def report_binding_fraction_real_seeds(seed_coords, ab_coords, cutoff_distance = 5.0):
	#simple method for calculating the binding fraction given an array of seed coordinates and an array
	#of antibody coordinates
	print "reporting binding fractions for real seeds"
	attached_seed_coordinates = []
	plt.clf()
	#fig, ax = plt.subplots(figsize=(3, 3))
	#ax.set_aspect('equal')
	for seed_coord in seed_coords:
		print "seed coords are: ", seed_coord[1], seed_coord[0]
	for seed_coord in seed_coords:
		print "seed coords are: ", seed_coord[1], seed_coord[0]
		for ab_coord in ab_coords:
			distance = calc_distance(seed_coord, ab_coord)
			print "distance is: ", distance
			print "ab coords are: ", ab_coord[1], ab_coord[0]
			#c = plt.Circle((ab_coord[1], ab_coord[0]), radius = 3, color='yellow', linewidth=2, fill=False)
	        	#ax.add_patch(c) 
			if distance <= cutoff_distance:
				#print "attached seed detected"
				print "distance is: ", distance 
				attached_seed_coordinates.append(seed_coord)
				break
	plt.savefig("plots/distance_test.pdf")
	if len(seed_coords) == 0:
		attached_seed_percentage = 0.0
	else:
		attached_seed_percentage = float(len(attached_seed_coordinates))/float(len(seed_coords))

	#print seed_coords
	#print "attached seed percentage is: ", attached_seed_percentage
	#print "total seeds: ", len(seed_coords)
	#print "total seeds with abs is: ", len(attached_seed_coordinates)
	#print "total abs: ", len(ab_coords)

	return attached_seed_percentage

def count_abs_in_membrane(ab_filename, seed_filename):
	#goal here is to calculate the number of abs present in the membrane
	image_ab = io.imread(ab_filename)
	ab_blobs = blob_log(image_ab, max_sigma=30, num_sigma=10, threshold=.005)
	ab_coords_list = []

	blobs_list = [ab_blobs]
	colors = ['yellow']
	titles = ['Laplacian of Gaussian']
	sequence = zip(blobs_list, colors, titles)
	fig, ax = plt.subplots(figsize=(3, 3))
	ax.set_aspect('equal')
	#ax = axes.ravel()





	#processing seed channel to outline cell and determine membrane area
	image_seed = io.imread(seed_filename)
	seed_blobs = blob_log(image_seed, max_sigma=30, num_sigma=10, threshold=.005)
	seed_coords_list = []



	#now we will perform edge detection on the seed image to find the perimeter of the cell
	#perfoming edge detection and morphological filling
	#hard coded threshold will be 4334
	mask = image_seed < 10 #originally 5000
	image_seed[mask] = 0
	edges_open = canny(image_seed, 5, 2, 25) #originally 2,1,25 last param can go up to 500 for improved performance, must lower for poorer images
	#edges_open = canny(image, 2) #originally 2,1,25
	selem = disk(30)#originally 5
	edges = closing(edges_open, selem) 
	edges = img_as_uint(edges)
  
  	#now we will fill in the closed cell and then re-detect edges to find the outline of the cell 
	fill_cell = ndi.binary_fill_holes(edges)

	#we only want to analyze the largest closed object in the image 
	label_image = label(fill_cell)
	largest_area = 0.0
	for region in regionprops(label_image):
		area = region.area
		if area > largest_area:
			cell_region = region 
			largest_area = area 
		else:
			continue
	#now cell region should contain just the outlined cell 
	img_cell = np.zeros((512,512), dtype=np.uint8)
	cell_coords_list = cell_region.coords.tolist()
	for coord in cell_coords_list:
		img_cell[coord[0], coord[1]] = 1

	print img_cell 
	#instead of edge detection we will apply a kernel to find the edge pixels in img_cell 

	#now we make a kernel to compute the endpoints of the skeletonized image
	kernel = np.uint8([[1,  1, 1], [1, 10, 1], [1,  1, 1]])

	#now we convolve the kernel with the skeletonized image 
	convolved_skeleton = ndimage.convolve(img_cell, kernel, mode='constant', cval = 1)

	#now produce an output mask with only pixels with value 11, these are the endpoints
	edge_mask = np.zeros_like(convolved_skeleton)
	#print edge_mask 
	#print convolved_skeleton
	
	edge_mask[np.where(convolved_skeleton != 18)] = 1
	edge_mask[np.where(convolved_skeleton == 0)] = 0
	#try performing edge detection on just largest region
	#edges_2 = canny(img_cell, 5, 2, 25)
	#print edges_2
	#edges_2 = canny(fill_cell, 5, 2, 25)

	#gathering the indices of the edge pixels from the binary image
	row_indices, column_indices = np.nonzero(edge_mask)
	print "length of row indices: ", len(row_indices)
	edge_indices = []
	for i in range(len(row_indices)):
		if row_indices[i] == 0 or row_indices[i] == 511:
			continue
		if column_indices[i] == 0 or column_indices[i] == 511:
			continue
		edge_indices.append([row_indices[i], column_indices[i]])

	#io.imsave("edge_test.pdf", img_as_uint(edges))
	#ax.imshow(img_as_uint(fill_tubes))
	#plt.savefig("plots/"+ab_filename+"processed.pdf")

	ax.imshow(img_cell, interpolation='nearest', cmap=cm.Greys, alpha = .6)
	ax.set_yticklabels([])
	ax.set_xticklabels([])
	#now for each edge pixel draw a circle of a cutoff radius
	#this should match the cutoff distance for considering an ab to be part of the membrane
	membrane_radius = 1 #should be 1
	img_membrane = np.zeros((512,512), dtype=np.uint8)
	print "number of edge points: ", len(edge_indices)
	for edge_point in edge_indices:
		#print edge_point
		rr, cc = circle(edge_point[0], edge_point[1], membrane_radius)
		img_membrane[rr, cc] = 1

		c = plt.Circle((edge_point[1], edge_point[0]), membrane_radius, color='blue', linewidth=2, fill=True)
		ax.add_patch(c) 
	
	#plotting abs
	#for idx, (blobs, color, title) in enumerate(sequence):
	    #ax[idx].set_title(title)
	for blob in ab_blobs:
	    y, x, r = blob
	    if r <= 2.0:
	    	c = plt.Circle((x, y), r*4, color='green', linewidth=0, fill=True)
	    	ax.add_patch(c)

	#plotting seeds
	seed_blobs_list = [seed_blobs]

	for blob in seed_blobs:
	#print blob 
		y, x, r = blob
		if r <= 2.0:
			c = plt.Circle((x, y), r*8, color='red', linewidth=0, fill=True, alpha = .40)
			ax.add_patch(c)
	plt.savefig("plots/"+ab_filename+"processed_just_seeds.pdf")

	for ab in ab_blobs:
		if ab[2] <= 2.0:#only consider blobs with small radii, these are the seeds/abs
			ab_coordinates = ab[0:2]
			ab_coords_list.append(ab_coordinates)
	ab_coords_list = np.vstack(ab_coords_list)
	#print ab_coords_list

	for seed in seed_blobs:
		if seed[2] <= 2.0:#only consider blobs with small radii, these are the seeds/abs
			seed_coordinates = seed[0:2]
			seed_coords_list.append(seed_coordinates)
	if len(seed_coords_list) == 0:
		return "no seeds", "no seeds"
	seed_coords_array = np.vstack(seed_coords_list)


	#culling the antibody points for counting 
	#first we will remove any point that does not have at least 4 points close by 
	culled_ab_points = []
	for point1 in ab_coords_list:
		n_neighbors = 0
		include_point = False
		for point2 in ab_coords_list:
			if point1[0] == point2[0] and point1[1] == point2[1]:
				continue
			p2p_distance = calc_distance(point1, point2)
			if p2p_distance <= 150:
				n_neighbors+=1
			if n_neighbors>=3:
				include_point = True 
				break

		if include_point == True:
			culled_ab_points.append(point1)
		else:
			print "excluding point!"

	#print culled_ab_points

	culled_ab_points_array = np.vstack(culled_ab_points)
	#print culled_ab_points_array

	#now we will only count abs that are within a cutoff distance of one of the edge pixels
	membrane_radius = 10 
	culled_ab_points_on_edge_for_real_calculation = []
	for point1 in culled_ab_points_array:
		include_point = False
		for point2 in edge_indices:
			if point1[0] == point2[0] and point1[1] == point2[1]:
				culled_ab_points_on_edge_for_real_calculation.append(point1)
				continue
			p2p_distance = calc_distance(point1, point2)
			#print "distance is: ", p2p_distance
			if p2p_distance <= membrane_radius:
				print "close to edge, counting this one"
				include_point = True 
				culled_ab_points_on_edge_for_real_calculation.append(point1)
				break

	#now we will only count seeds that are within a cutoff distance of one of the edge pixels
	membrane_radius = 10 
	culled_seeds = []
	for point1 in seed_coords_array:
		include_point = False
		for point2 in edge_indices:
			if point1[0] == point2[0] and point1[1] == point2[1]:
				culled_seeds.append(point1)
				continue
			p2p_distance = calc_distance(point1, point2)
			#print "distance is: ", p2p_distance
			if p2p_distance <= membrane_radius:
				print "close to edge, counting this one"
				include_point = True 
				culled_seeds.append(point1)
				break
	if len(culled_seeds) == 0:
		return "no seeds", "no seeds"
	seed_coords_array = np.vstack(culled_seeds)

	culled_ab_points_array_on_edge_for_real_calculation = np.vstack(culled_ab_points_on_edge_for_real_calculation)
	culled_binding_fraction = report_binding_fraction_real_seeds(seed_coords_array, culled_ab_points_array_on_edge_for_real_calculation)

	membrane_radius = 1 
	culled_ab_points_on_edge = []
	for point1 in culled_ab_points_array:
		include_point = False
		for point2 in edge_indices:
			if point1[0] == point2[0] and point1[1] == point2[1]:
				culled_ab_points_on_edge.append(point1)
				continue
			p2p_distance = calc_distance(point1, point2)
			#print "distance is: ", p2p_distance
			if p2p_distance <= membrane_radius:
				print "close to edge, counting this one"
				include_point = True 
				culled_ab_points_on_edge.append(point1)
				break

	culled_ab_points_array_on_edge = np.vstack(culled_ab_points_on_edge)
	

	plt.clf()
	fig, ax = plt.subplots(figsize=(3, 3))
	ax.set_aspect('equal')
	ax.imshow(img_cell, interpolation='nearest')

	for blob in seed_blobs:
	#print blob 
		y, x, r = blob
		if r <= 2.0:
			c = plt.Circle((x, y), r*8, color='red', linewidth=1, fill=True, alpha = .5)
			ax.add_patch(c)
	for culled_ab in culled_ab_points_array_on_edge:
	#print blob 
		y, x = culled_ab[0], culled_ab[1]
		c = plt.Circle((x, y), 5, color='green', linewidth=1, fill=True)
		ax.add_patch(c)

	plt.savefig("plots/"+ab_filename+"post_culling_processed_just_seeds.pdf")

	#culled_ab_points_array_on_edge = np.vstack(culled_ab_points_on_edge)
	#culled_binding_fraction = report_binding_fraction_real_seeds(seed_coords_array, culled_ab_points_array_on_edge)

	pairwise_distances = []
	for i in range(len(culled_ab_points_array)):
		for j in range(i+1, len(culled_ab_points_array)):
			point1 = culled_ab_points_array[i]
			point2 = culled_ab_points_array[j]
			pairwise_distance = calc_distance(point1, point2)
			pairwise_distances.append(pairwise_distance)
	average_pairwise_distance = sum(pairwise_distances)/float(len(pairwise_distances))

		
	print culled_ab_points_array_on_edge
	print len(culled_ab_points_array_on_edge)


	#now for each edge pixel draw a circle of a cutoff radius
	#this should match the cutoff distance for considering an ab to be part of the membrane
	row_indices_mem, column_indices_mem = np.where(img_membrane)
	mem_indices = np.column_stack((row_indices_mem,column_indices_mem))
	#now we will randomly pick a number of points in mem_indices equal to the number of seeds
	#and calculate the binding fraction
	#we will do this 1000 times and report the average null binding
	null_fractions = []
	print "number of detected seeds: ", len(seed_coords_array)
	for i in range(1000):
		rand_indices = np.random.randint(len(mem_indices), size=(len(seed_coords_array)))
		simulated_seed_coords_array = []
		#print "simulated seed coords: ", simulated_seed_coords_array
		for index in rand_indices:
			simulated_seed_coords_array.append(mem_indices[index])
		#print "simulated seed coords: ", simulated_seed_coords_array
		null_fraction = report_binding_fraction(simulated_seed_coords_array, culled_ab_points_array_on_edge)
		#print "null fraction: ", null_fraction 
		null_fractions.append(null_fraction)
	#print null_fractions
	null_binding_fraction = sum(null_fractions)/len(null_fractions)


	return culled_binding_fraction, null_binding_fraction 
	

	


#time_dirs = ['3', '4', '5', '6', '7', '8', '9', '10', '11', '12']
#time_dirs = [ '1', '2', '3', '4']
time_dirs = ['test']
time_ab_membrane_counts = []
time_seed_pairwise_distance_avgs = []


for time_dir in time_dirs:	
	ab_membrane_counts = []
	culled_binding_fractions = []
	null_binding_fractions = []
	ab_file_list = glob.glob(time_dir+"/*AB*")
	seed_file_list = glob.glob(time_dir+"/*S*")
	print ab_file_list
	print seed_file_list

	for i in range(len(ab_file_list)):
		ab_full_file = ab_file_list[i]
		seed_full_file = seed_file_list[i]
		print 'ab file is: ', ab_full_file
		print 'seed file is: ', seed_full_file
		#try:
		culled_binding_fraction, null_binding_fraction = count_abs_in_membrane(ab_full_file, seed_full_file)
			#ab_area_fraction = detect_ab_null_membrane_fraction(ab_full_file, seed_full_file)
		#except ValueError:
		#	continue 
		#except IndexError:
		#	continue 
		if culled_binding_fraction == "no seeds":
			continue 

		culled_binding_fractions.append(culled_binding_fraction)
		null_binding_fractions.append(null_binding_fraction)
	
	
	total_binding_fraction = sum(culled_binding_fractions)/float(len(culled_binding_fractions))
	total_null_binding_fraction = sum(null_binding_fractions)/float(len(null_binding_fractions))
	print "printing culled_binding_fractions"
	f1=open(time_dir+'_culled_binding_fraction_optimal_membrane.dat','w+')
	#for fraction in culled_binding_fractions:
		#print >>f1, fraction
	print >>f1, total_binding_fraction
	f1.close()

	print "printing null_binding_fractions"
	f1=open(time_dir+'_null_binding_fraction_optimal_membrane.dat','w+')
	#for fraction in null_binding_fractions:
	#	print >>f1, fraction
	print >>f1, total_null_binding_fraction
	f1.close()

'''print "printing average pairwise distance for time points"
f1=open(time_dir+'_pairwise_distance.dat','w+')
for dist in seed_pairwise_distance_avg:
	print >>f1, dist  
f1.close()'''
#times = []
#for time in time_dirs:
#	times.append(float(time))

#print times
#print time_ab_membrane_counts

'''plt.clf()
plt.plot(times, time_ab_membrane_counts)
plt.xlabel('time (hours)')
plt.ylabel('antibody count on membrane')
plt.savefig("membrane_count_vs_time")'''	



