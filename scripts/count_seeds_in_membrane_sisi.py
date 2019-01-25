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
			if distance <= cutoff_distance:
				#print "attached seed detected"
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

def count_abs_in_membrane(ab_filename, seed_filename):
	#goal here is to calculate the fraction of the membrane area occupied by abs
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

	#plotting abs
	for idx, (blobs, color, title) in enumerate(sequence):
	    #ax[idx].set_title(title)
	    ax.imshow(image_ab, interpolation='nearest')
	    for blob in blobs:
	        y, x, r = blob
	        if r <= 2.0:
	        	c = plt.Circle((x, y), r, color='yellow', linewidth=2, fill=False)
	        	ax.add_patch(c)
	#plotting seeds 
	image_seed = io.imread(seed_filename)
	seed_blobs = blob_log(image_seed, max_sigma=30, num_sigma=10, threshold=.005)
	seed_coords_list = []
	for seed in seed_blobs:
		if seed[2] <= 2.0:#only consider blobs with small radii, these are the seeds/abs
			seed_coordinates = seed[0:2]
			seed_coords_list.append(seed_coordinates)
	seed_coords_array = np.vstack(seed_coords_list)

	'''for blob in seed_blobs:
		#print blob 
		y, x, r = blob
		if r <= 2.0:
			c = plt.Circle((x, y), r*8, color='green', linewidth=1, fill=False)
			ax.add_patch(c)'''
			#plt.plot(perimeter_points[0], 'r--', lw = 2)
	#plt.savefig("plots/"+ab_filename+"processed.pdf")

	#now we will perform edge detection on the seed image to find the perimeter of the cell

	#perfoming edge detection and morphological filling
	#hard coded threshold will be 4334
	mask = image_seed < 5000
	image_seed[mask] = 0
	edges_open = canny(image_seed, 5, 2, 25) #originally 2,1,25 last param can go up to 500 for improved performance, must lower for poorer images
	#edges_open = canny(image, 2) #originally 2,1,25
	selem = disk(10)#originally 5
	edges = closing(edges_open, selem)
	print edges 
	edges = img_as_uint(edges)
	print edges 
	#print "edge pixels are: ", edges 
	fill_tubes = ndi.binary_fill_holes(edges)
	edges_2 = canny(fill_tubes, 5, 2, 25)
	row_indices, column_indices = np.nonzero(edges_2)
	print "length of row indices: ", len(row_indices)
	edge_indices = []
	for i in range(len(row_indices)):
		edge_indices.append([row_indices[i], column_indices[i]])
	#edge_indices = np.column_stack((row_indices,column_indices))

	#print "edge indices are: ", edge_indices
	#fill_tubes = ndi.binary_fill_holes(edges)
		#io.imsave(str(i)+"_fill_tubes.png", img_as_uint(fill_tubes), cmap=cm.gray)
	io.imsave("edge_test.pdf", img_as_uint(edges))
	ax.imshow(img_as_uint(fill_tubes))
	plt.savefig("plots/"+ab_filename+"processed.pdf")
	for ab in ab_blobs:
		if ab[2] <= 2.0:#only consider blobs with small radii, these are the seeds/abs
			ab_coordinates = ab[0:2]
			ab_coords_list.append(ab_coordinates)
	ab_coords_list = np.vstack(ab_coords_list)
	#print ab_coords_list


	#culling the antibody points for null analysis 
	#first we will remove any point that does not have at least 4 points close by 
	culled_ab_points = []
	for point1 in ab_coords_list:
		n_neighbors = 0
		include_point = False
		for point2 in ab_coords_list:
			if point1[0] == point2[0] and point1[1] == point2[1]:
				continue
			p2p_distance = calc_distance(point1, point2)
			if p2p_distance <= 50:
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
	culled_ab_points_array_on_edge = []
	for point1 in culled_ab_points_array:
		include_point = False
		for point2 in edge_indices:
			if point1[0] == point2[0] and point1[1] == point2[1]:
				continue
			p2p_distance = calc_distance(point1, point2)
			#print "distance is: ", p2p_distance
			if p2p_distance <= 5.0:
				print "close to edge, counting this one"
				include_point = True 
				culled_ab_points_array_on_edge.append(point1)
				break
		
	print culled_ab_points_array_on_edge
	print len(culled_ab_points_array_on_edge)

	img_membrane = np.zeros((512,512), dtype=np.uint8)

	#now for each edge pixel draw a circle of a cutoff radius
	#this should match the cutoff distance for considering an ab to be part of the membrane
	membrane_radius = 5
	print "number of edge points: ", len(edge_indices)
	for edge_point in edge_indices:
		#print edge_point
		rr, cc = circle(edge_point[0], edge_point[1], membrane_radius)
		img_membrane[rr, cc] = 1
		c = plt.Circle((edge_point[1], edge_point[0]), membrane_radius, color='green', linewidth=2, fill=True)
		ax.add_patch(c) 
	plt.savefig("plots/"+ab_filename+"processed_with_membrane.pdf")
	return len(culled_ab_points_array_on_edge)

	


#time_dirs = [ '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18']
time_dirs = [ '1']
time_ab_membrane_counts = []

for time_dir in time_dirs:	
	ab_membrane_counts = []
	ab_file_list = glob.glob(time_dir+"/*1.tif")
	seed_file_list = glob.glob(time_dir+"/*2.tif")
	print ab_file_list
	print seed_file_list

	for i in range(len(ab_file_list)):
		ab_full_file = ab_file_list[i]
		seed_full_file = seed_file_list[i]
		print 'ab file is: ', ab_full_file
		print 'seed file is: ', seed_full_file
		try:
			ab_membrane_count = count_abs_in_membrane(ab_full_file, seed_full_file)
			#ab_area_fraction = detect_ab_null_membrane_fraction(ab_full_file, seed_full_file)
		except ValueError:
			continue 
		except IndexError:
			continue 
		

		ab_membrane_counts.append(ab_membrane_count)
	
	time_ab_membrane_counts.append(sum(ab_membrane_counts))

times = []
for time in time_dirs:
	times.append(float(time))

print times
print time_ab_membrane_counts

plt.clf()
plt.plot(times, time_ab_membrane_counts)
plt.xlabel('time (hours)')
plt.ylabel('antibody count on membrane')
plt.savefig("membrane_count_vs_time")	




