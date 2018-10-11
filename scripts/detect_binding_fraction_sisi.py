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
from scipy import ndimage

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

def calc_pereimeter_points(points):
	#given a list of points that we know are near the edge of the cell, return points
	#that can be used to calculate the perimeter of the cell
	return 0 

def concave(points,alpha_x=150,alpha_y=250):
    points = [(i[0],i[1]) if type(i) <> tuple else i for i in points]
    de = Delaunay(points)
    dec = []
    a = alpha_x
    b = alpha_y
    for i in de.simplices:
        tmp = []
        j = [points[c] for c in i]
        if abs(j[0][1] - j[1][1])>a or abs(j[1][1]-j[2][1])>a or abs(j[0][1]-j[2][1])>a or abs(j[0][0]-j[1][0])>b or abs(j[1][0]-j[2][0])>b or abs(j[0][0]-j[2][0])>b:
            continue
        for c in i:
            tmp.append(points[c])
        dec.append(tmp)
    G = nx.Graph()
    for i in dec:
            G.add_edge(i[0], i[1])
            G.add_edge(i[0], i[2])
            G.add_edge(i[1], i[2])
    ret = []
    for graph in nx.connected_component_subgraphs(G):
        ch = ConvexHull(graph.nodes())
        tmp = []
        for i in ch.simplices:
            tmp.append(graph.nodes()[i[0]])
            tmp.append(graph.nodes()[i[1]])
        ret.append(tmp)
    return ret 

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


def detect_ab_null_membrane_fraction(ab_filename, seed_filename):
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

	for idx, (blobs, color, title) in enumerate(sequence):
	    #ax[idx].set_title(title)
	    ax.imshow(image_ab, interpolation='nearest')
	    for blob in blobs:
	        y, x, r = blob
	        if r <= 2.0:
	        	c = plt.Circle((x, y), r, color=color, linewidth=2, fill=False)
	        	ax.add_patch(c)

	for ab in ab_blobs:
		if ab[2] <= 2.0:#only consider blobs with small radii, these are the seeds/abs
			ab_coordinates = ab[0:2]
			ab_coords_list.append(ab_coordinates)
	ab_coords_list = np.vstack(ab_coords_list)
	print ab_coords_list

	#convex hull does not see to be good enough here...
	'''hull = ConvexHull(ab_coords_list, True)
	print hull.vertices
	#fig, ax = plt.subplots(figsize=(9, 3))
	#ax.set_aspect('equal')
	#ax.imshow(image_ab, interpolation='nearest')
	#plt.imshow(image_ab)
	for simplex in hull.simplices:
		plt.plot(ab_coords_list[simplex, 1], ab_coords_list[simplex, 0], 'k-')
	for vertex in hull.vertices:
		c = plt.Circle((ab_coords_list[vertex][1], ab_coords_list[vertex][0]), radius = 5.0, color = 'red', fill=False)
		ax.add_patch(c)
	'''

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

	print culled_ab_points

	culled_ab_points_array = np.vstack(culled_ab_points)
	print culled_ab_points_array

	tri = Delaunay(culled_ab_points_array)
	plt.triplot(culled_ab_points_array[:,1], culled_ab_points_array[:,0], tri.simplices.copy(), color = 'black', alpha = .5)
	#com = ndimage.measurements.center_of_mass(culled_ab_points_array)
	#print com  
	culled_ab_points.append(np.array([250., 250.]))
	culled_ab_points.append(np.array([250., 255.]))
	culled_ab_points.append(np.array([250., 245.]))
	culled_ab_points.append(np.array([255., 250.]))
	culled_ab_points.append(np.array([245., 250.]))

	perimeter_points = concave(culled_ab_points, alpha_x = 75, alpha_y = 75)
	#can detect perimeter from low intensity pixels on cy3 channel!
	#edges_open = canny(image, 5, 2, 25)
	#selem = disk(10)
	print perimeter_points[0]
	for perimeter_point in perimeter_points[0]:
		y, x = perimeter_point[0], perimeter_point[1]
		c = plt.Circle((x, y), r*8, color='red', linewidth=1, fill=False)
		ax.add_patch(c)
	
	#plt.plot(perimeter_points[0], 'r--', lw=2)
	#plt.plot(ab_coords_list[hull.vertices[0],0], ab_coords_list[hull.vertices[0],1], 'ro')
	#ax.add_patch(hull_vertices)
	#ax.add_patch(hull_lines)

	image_seed = io.imread(seed_filename)
	seed_blobs = blob_log(image_seed, max_sigma=30, num_sigma=10, threshold=.005)
	for blob in seed_blobs:
		#print blob 
		y, x, r = blob
		c = plt.Circle((x, y), r*8, color='green', linewidth=1, fill=False)
		ax.add_patch(c)
	#plt.plot(perimeter_points[0], 'r--', lw = 2)
	plt.savefig("plots/"+ab_filename+"processed.pdf")

	#now we will compute the perimeter of the convex hull
	'''last_index = len(perimeter_points) - 1
	previous_point = perimeter_points[last_index]
	total_distance = 0.0
	for point in hull.vertices:
		point = ab_coords_list[vertex]
		distance = calc_distance(point, previous_point)
		total_distance+=distance
		previous_point = point

	print "total perimiter distance: ", total_distance

	#now compute the membrane "area" by multiplying the perimeter by the "width" which is the cutoff distance
	cutoff_distance = 2.0
	membrane_area = total_distance * cutoff_distance

	#now the total area taken by abs, the sum of all their circular areas
	ab_area = math.pi * cutoff_distance * cutoff_distance
	total_ab_area = ab_area * float(len(ab_coords_list))

	ab_area_fraction = total_ab_area/membrane_area
	print "ab area fraction is: ", ab_area_fraction

	return ab_area_fraction
	'''
	return 1.0 


#detect_ab_seed_binding_fraction("1_AB_00000001.tif", "1_S_00000001.tif")
#detect_ab_null_membrane_fraction("1_AB_00000001.tif")
#cell_dirs = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']


cell_dirs = ['1','2','3']
attached_seed_percentages = []
ab_area_fractions = []
for cell_dir in cell_dirs:	
	ab_file_list = glob.glob(cell_dir+"/*AB*")
	seed_file_list = glob.glob(cell_dir+"/*S*")
	print ab_file_list
	print seed_file_list

	for i in range(len(ab_file_list)):
		ab_full_file = ab_file_list[i]
		seed_full_file = seed_file_list[i]
		print 'ab file is: ', ab_full_file
		print 'seed file is: ', seed_full_file
	
		attached_seed_percentage = detect_ab_seed_binding_fraction(ab_full_file, seed_full_file)
		ab_area_fraction = detect_ab_null_membrane_fraction(ab_full_file, seed_full_file)
		print ab_area_fraction
		print attached_seed_percentage

		attached_seed_percentages.append(attached_seed_percentage)
		ab_area_fractions.append(ab_area_fraction)

	print "printing attached seed percentages"
	f1=open(cell_dir+'_attached_seed_percentages.dat','w+')
	for percentage in attached_seed_percentages:
		print >>f1, percentage
	f1.close()
	
	print "printing ab area fractions"
	f1=open(cell_dir+'_ab_area_fractions.dat','w+')
	for percentage in ab_area_fractions:
		print >>f1, percentage
	f1.close()


