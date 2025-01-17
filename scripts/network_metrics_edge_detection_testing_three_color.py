import numpy as np

from skimage.transform import (hough_line, hough_line_peaks,
                               probabilistic_hough_line)
from skimage.feature import canny
from skimage import data
from skimage import io 
from skimage.morphology import closing, disk
from skimage.morphology import skeletonize
from skimage.measure import label, regionprops, find_contours
from skimage.filters import threshold_otsu, threshold_local, rank
import sys
import os
import matplotlib.pyplot as plt
import math
from matplotlib import cm
from matplotlib import path
from skimage import img_as_uint, img_as_bool, img_as_float, img_as_ubyte
from scipy import ndimage
from scipy.spatial import distance
from scipy import ndimage as ndi
from numpy import unravel_index
import Tkinter, tkFileDialog
from skimage.external import tifffile

from math import sqrt
from skimage import data
from skimage.feature import blob_dog, blob_log, blob_doh
from skimage.color import rgb2gray
from skimage import io 
import sys
from skimage.data import camera
from skimage.filters import roberts, sobel, scharr, prewitt
from skimage.exposure import equalize_hist 
#from math import mean

def count_seeds(filename):
	"""count the number of seeds in a grayscale image using the laplacian
	of gaussians method"""

	image_gray = io.imread(filename)
	image = io.imread(filename)
	#image_gray = rgb2gray(image)

	blobs_log = blob_log(image_gray, max_sigma=30, num_sigma=10, threshold=.01)
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

	#plt.tight_layout()
	plt.show()



def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
  return math.acos(abs(dotproduct(v1, v2) )/ (length(v1) * length(v2)))

def line_length(line):
	p0, p1 = line
	a = np.array((p0[0],p0[1]))
	b = np.array((p1[0],p1[1]))
	dist = np.linalg.norm(a-b)
	#print dist
	return dist

def make_endpoints_mask(filled_binary_image):
	#function to determine the endpoints of a nanotube identified via edge detection/morphological filling
	#need to find all endpoint candidates and find the pair separated by the longest path

	#first skeletonize the filled binary image (must be a binary int image)
	filled_binary_image = filled_binary_image.astype(int)
	skeleton = skeletonize(filled_binary_image)
	skeleton = skeleton.astype(int)


	#now we make a kernel to compute the endpoints of the skeletonized image
	kernel = np.uint8([[1,  1, 1], [1, 10, 1], [1,  1, 1]])

	#now we convolve the kernel with the skeletonized image 
	convolved_skeleton = ndimage.convolve(skeleton, kernel, mode='constant', cval = 1)

	#now produce an output mask with only pixels with value 11, these are the endpoints
	endpoint_mask = np.zeros_like(convolved_skeleton)
	endpoint_mask[np.where(convolved_skeleton == 11)] = 1

	return endpoint_mask, skeleton

def endpoints(region_coords, endpoint_mask):
	#using a previously genereated endpoint mask to find the endpoints for a particular tube
	#this will return a pair of tubles with the x,y coordinates of the two endpoints 
    endpoints_labelled = label(endpoint_mask)
    potential_endpoints = []
    for endpoint in regionprops(endpoints_labelled):
    	if any(i in region_coords for i in endpoint.coords.tolist()):
    		potential_endpoints.append(endpoint.centroid)
    
    #now we will find the pair of potential endpoints with the maximal separation distance, those are the true endpoints
    if len(potential_endpoints) <= 1:
    	return None 

    pairwise_distances = distance.cdist(potential_endpoints, potential_endpoints, 'euclidean')
    indices_of_max_distance = unravel_index(pairwise_distances.argmax(), pairwise_distances.shape)

    endpoint1 = potential_endpoints[indices_of_max_distance[0]]
    endpoint2 = potential_endpoints[indices_of_max_distance[1]]
    #print endpoint1
    #print endpoint2
    endpoints = [endpoint1, endpoint2]
    return endpoints

def are_joined(endpoint1, endpoint2):
	#given two endpoints calculate the distance between them and return True or False for whether they meet the joining criteria
	cutoff = 5.0 
	distance = distance(endpoint1,endpoint2)
	if distance <= cutoff: 
		return True 

	else:
		return False 


def calc_distance(endpoint1, endpoint2):
	#simple distance calculation
	distance_squared = (endpoint1[0]-endpoint2[0]) * (endpoint1[0]-endpoint2[0]) + (endpoint1[1]-endpoint2[1]) * (endpoint1[1]-endpoint2[1])
	distance = math.sqrt(distance_squared)

	return distance

def generate_ordered_coords_from_skeleton(skeleton, endpoint1, endpoint2):
	#generate ordered coords (parametric curve) from a skeleton and its endpoints
	ordered_coords = [endpoint1]
	skeleton_point = endpoint1
	#neighbors = search_eight_neighbors(endpoint1, skeleton)
	while skeleton_point != endpoint2: 
		if skeleton_point == endpoint1: #if we are at the first endpoint
			print "at the first endpoint: ", skeleton_point
			positive_neighbors = search_eight_neighbors(skeleton_point, skeleton)
			print "these are the positive neighbors: ", positive_neighbors
			next_point = positive_neighbors[0]
			previous_point = skeleton_point
			skeleton_point = next_point
			ordered_coords.append(skeleton_point)

		else: #we are at subsequent points
			print "at a later point: ", skeleton_point
			print "these are the positive neighbors: "
			positive_neighbors = search_eight_neighbors(skeleton_point, skeleton)
			print positive_neighbors
			for point in positive_neighbors:
				if point == previous_point:
					continue
				else: 
					next_point = point
					previous_point = skeleton_point
					skeleton_point = next_point
					ordered_coords.append(skeleton_point)
		#print "skeleton point: ", skeleton_point
	

	#ordered_coords.append(endpoint2)
	print ordered_coords
	return ordered_coords 




		#start at endpoint 1 and search all 8 neighbors for a positive pixel

def search_eight_neighbors(pixel, skeleton):
	#search the surrounding 8 neighbors and return a list of positive valued pixels
	positive_neighbor_list = []
	print "searching for neighbors with positive pixels"
	for i in range(3):
		for j in range(3):
			print "pixel location: ", int(pixel[0]) - 1 + i, int(pixel[1]) - 1 + j
			print "value is: ", skeleton[int(pixel[0]) - 1 + i][int(pixel[1]) - 1 + j]
			if i == 1 and j == 1:
				continue
			elif skeleton[int(pixel[0]) - 1 + i][int(pixel[1]) - 1 + j] == 1:
				print "found a positive pixel!"
				positive_neighbor_list.append((int(pixel[0]) - 1 + i, int(pixel[1]) - 1 + j))
	#print "positive neighbor coords: "
	#print positive_neighbor_list
	return positive_neighbor_list

def calc_mean_tangent_from_ordered_coords(ordered_coords):
	print  "nothing"
	#given a list of ordered coordinates from a skeleton, calculat the mean tangent line 

def count_seeds_per_network(filename, filename_seeds_1, filename_seeds_2, plot_filename):
	#modifying this to account for both 488 and 647 seeds 
	number_of_seeds = []
	image= io.imread(filename)


			#cy3_file = cy3_file_list[i]

			#print "cy3 filename is "+str(cy3_file)
			#image_unthresholded = io.imread(cy3_file)

			#thresh = threshold_otsu(image_unthresholded)
			#image = image_unthresholded>thresh

		
			#image = threshold_local(image_unthresholded, block_size, offset=10)
			#image_647 = threshold_local(image_647_unthresholded, block_size, offset=10)

	radius = 1
	selem = disk(radius)

			#thresholding both files (getting rid of this because it should not be necessary!)
			#image = rank.otsu(image_unthresholded, selem)
			#image_647 = rank.otsu(image_647_unthresholded, selem)

			#image = image_unthresholded


			#perfoming edge detection and morphological filling
	edges_open = canny(image, 1, 15, 100) #originally 2,1,25 last param can go up to 500 for improved performance, must lower for poorer images
			#edges_open = canny(image, 2) #originally 2,1,25
	selem = disk(3)#originally 5
	edges = closing(edges_open, selem)
	fill_tubes = ndi.binary_fill_holes(edges)
	io.imsave("fill_tubes.png", img_as_uint(fill_tubes), cmap=cm.gray)
	#io.imshow(img_as_uint(fill_tubes))

	cy3_endpoint_mask, skeleton = make_endpoints_mask(fill_tubes)
			#io.imsave(str(i)+"_skeleton.png", img_as_uint(img_as_bool(skeleton)), cmap=cm.gray)
			#print fill_tubes
		#print skeleton
		#print cy3_endpoint_mask
			
			

			#io.imsave(str(i)+"_skeleton.png", skeleton, cmap=cm.gray)



			#label image 
	label_image = label(fill_tubes)


	print "number of detected networks: ", len(regionprops(label_image))

	##########################################################################
	#now we will count/identify the seeds
	image_seeds1 = io.imread(filename_seeds_1)
	image_seeds2 = io.imread(filename_seeds_2)
	
	print "now detecting seeds"
	blobs_log1 = blob_log(image_seeds1, max_sigma=30, num_sigma=10, threshold=.1)
	blobs_log2 = blob_log(image_seeds2, max_sigma=30, num_sigma=10, threshold=.1)
	#blobs_log = blob_log(image_gray, max_sigma=1, num_sigma=10, threshold=.1)

	print "number of seeds 1: "+str(len(blobs_log1))
	print "number of seeds 2: "+str(len(blobs_log2))

		# Compute radii in the 3rd column.
	blobs_log1[:, 2] = blobs_log1[:, 2] * sqrt(2)
	blobs_log2[:, 2] = blobs_log2[:, 2] * sqrt(2)


	blobs_list1 = blobs_log1
	blobs_list2 = blobs_log2

	fig, axes = plt.subplots(1, 4, figsize=(12, 4), sharex = True, sharey = True)
	#ax.set_aspect('equal')
	ax = axes.ravel()

	ax[0].set_title("original_tubes")
	ax[0].imshow(image, cmap = "gray")

	ax[1].set_title("original_488_seeds")
	ax[1].imshow(image_seeds1, cmap = 'gray')

	ax[2].set_title("original_647_seeds")
	ax[2].imshow(image_seeds2, cmap = 'gray')

	ax[3].set_title("processed")

	ax[3].imshow(fill_tubes, interpolation='nearest', cmap = 'Greens')
	print blobs_list1 
	for blob in blobs_list1:
		print blob 
	 	y, x, r = blob
		c = plt.Circle((x, y), r, color='blue', linewidth=1, fill=True)
		ax[3].add_patch(c)
		#ax[idx].set_axis_off()

	for blob in blobs_list2:
		print blob 
	 	y, x, r = blob
		c = plt.Circle((x, y), r, color='red', linewidth=1, fill=True)
		ax[3].add_patch(c)
		#ax[idx].set_axis_off()

	#plt.tight_layout()
	plt.savefig(plot_filename)

	plt.clf()
	fig, ax = plt.subplots()
	ax.imshow(fill_tubes, interpolation='nearest', cmap = 'Greens')
	for blob in blobs_list1:
		print blob
	 	y, x, r = blob
		c = plt.Circle((x, y), r, color='blue', linewidth=1, fill=True)
		ax.add_patch(c)
		#ax[idx].set_axis_off()

	for blob in blobs_list2:
		print blob 
	 	y, x, r = blob
		c = plt.Circle((x, y), r, color='red', linewidth=1, fill=True)
		ax.add_patch(c)
		#ax[idx].set_axis_off()

	plt.savefig(plot_filename+"individual.png")
	#now we will merge the two blob lists
	print "blobs_log_1: ", blobs_log1
	print "blobs_log_2: ", blobs_log2
	#blobs_log = np.append(blobs_log1, blobs_log2)
	#print "blobs_log: ", blobs_log 
		
	for region in regionprops(label_image):
		region_coords_list = region.coords.tolist()
		#print "region coords: ", region_coords_list 
		seed_count = 0
		for blob in blobs_log1:
			seed_coords = [blob[0], blob[1]]
			#print "seed coords: ", seed_coords
			for region_coords in region_coords_list:
				distance = calc_distance(seed_coords, region_coords)
				if distance <=3:
					seed_count += 1
					break

		for blob in blobs_log2:
			seed_coords = [blob[0], blob[1]]
			#print "seed coords: ", seed_coords
			for region_coords in region_coords_list:
				distance = calc_distance(seed_coords, region_coords)
				if distance <=3:
					seed_count += 1
					break
		if seed_count > 0:
			for i in range(seed_count):
				number_of_seeds.append(seed_count)

	'''print image.dtype
	image = (image/256).astype('uint8')
	print image.dtype 
	image_seeds = (image_seeds/256).astype('uint8')
	
	io.imsave(filename+"_converted.tiff", image)
	io.imsave(filename_seeds+"_converted.tiff", image_seeds)'''
	print number_of_seeds
	return number_of_seeds





def process_tiff_stacks(filename):
	# Line finding using the Probabilistic Hough Transform
	tube_lengths = []
	#tube_angles = []


	i=0
	#cy3_file_list = os.listdir('6_nt')
	'''root = Tkinter.Tk()
	root.withdraw()

	file_paths = tkFileDialog.askopenfilenames()
	cy3_file_list = list(file_paths)
	'''

	cy3_image_stack = tifffile.imread(filename)

	for image in cy3_image_stack:
		tube_angles = []
		if (i+1)%2 != 0:
			i+=1
			continue
		else:
			total_images = len(cy3_image_stack)
			current_frame = i
			print "processing frame " +str(i) + " of "+str(total_images)

			#cy3_file = cy3_file_list[i]

			#print "cy3 filename is "+str(cy3_file)
			#image_unthresholded = io.imread(cy3_file)

			#thresh = threshold_otsu(image_unthresholded)
			#image = image_unthresholded>thresh

			block_size = 15
			#image = threshold_local(image_unthresholded, block_size, offset=10)
			#image_647 = threshold_local(image_647_unthresholded, block_size, offset=10)

			radius = 5
			selem = disk(radius)

			#thresholding both files (getting rid of this because it should not be necessary!)
			#image = rank.otsu(image_unthresholded, selem)
			#image_647 = rank.otsu(image_647_unthresholded, selem)

			#image = image_unthresholded


			#perfoming edge detection and morphological filling
			edges_open = canny(image, .05, .5, 5) #originally 2,1,25 last param can go up to 500 for improved performance, must lower for poorer images
			#edges_open = canny(image, 2) #originally 2,1,25
			selem = disk(3)#originally 5
			edges = closing(edges_open, selem)
			fill_tubes = ndi.binary_fill_holes(edges)
			io.imsave(str(i)+"_fill_tubes.png", img_as_uint(fill_tubes), cmap=cm.gray)
			cy3_endpoint_mask, skeleton = make_endpoints_mask(fill_tubes)
			#io.imsave(str(i)+"_skeleton.png", img_as_uint(img_as_bool(skeleton)), cmap=cm.gray)
			#print fill_tubes
			print skeleton
			print cy3_endpoint_mask
			
			

			#io.imsave(str(i)+"_skeleton.png", skeleton, cmap=cm.gray)



			#label image 
			label_image = label(fill_tubes)


			print "detecting nanotube angles...."
			print len(regionprops(label_image))
			for region in regionprops(label_image):
				if region.area/tube_width >= length_cutoff and region.eccentricity >= eccentricity_cutoff:
					region_coords = region.coords.tolist()
					region_endpoints = endpoints(region_coords, cy3_endpoint_mask)

					if region_endpoints == None:
						continue

					print region_endpoints
					#print region_endpoints[0][0]
					#print region_endpoints[0][1]
					#print region_endpoints[1][0]
					#print region_endpoints[1][1]
					#cropped_skeleton = skeleton[int(region_endpoints[0][0]):int(region_endpoints[1][0]), int(region_endpoints[1][1]):int(region_endpoints[0][1])]
					#print "cropped skeleton"
					#print cropped_skeleton
					#contour = find_contours(cropped_skeleton, level = 1.0)
					#print contour 
					generate_ordered_coords_from_skeleton(skeleton, region_endpoints[0], region_endpoints[1])
					endpoint_to_endpoint_vector  = np.subtract(region_endpoints[0], region_endpoints[1])
					#x_axis_vector = np.array([0, 1])
					x_axis_vector = np.array([1,0])
					angle_with_x_axis = angle(endpoint_to_endpoint_vector, x_axis_vector)
					angle_with_x_axis *= 180.0/math.pi
					print 'angle with x axis is: ', angle_with_x_axis
					tube_angles.append(angle_with_x_axis)

						
			i+=1



		print "printing angles"
		f1=open('angles.dat','a')
		for angle_ in tube_angles:
			print >>f1, angle_
		f1.close()
#modifying the joining detection script to measure the angle of Sisi's nanotubes relative to the x-axis of her images 

#constants
#constants
tube_width = 5.0
length_cutoff = 3.0 
eccentricity_cutoff = 0.5
end_to_end_distance_cutoff = 10.0

#count_seeds_per_network("Lime 3.tif", "Blue 3.tif")

time_dirs = ['8']

network_sizes = []
for time_dir in time_dirs:
	network_size_time = []	
	tube_file_list = os.listdir(time_dir+'/tubes')
	seed1_file_list = os.listdir(time_dir+'/seeds1')
	seed2_file_list = os.listdir(time_dir+'/seeds2')

	for i in range(len(tube_file_list)):
		tube_file = tube_file_list[i]
		seed1_file = seed1_file_list[i]
		seed2_file = seed2_file_list[i]
		tube_full_file = time_dir+"/"+"tubes/"+tube_file
		seed_full_file_1 = time_dir+"/"+"seeds1/"+seed1_file
		seed_full_file_2 = time_dir+"/"+"seeds2/"+seed2_file
		print 'tube file is: ', tube_full_file
		print 'seed 1 file is: ', seed_full_file_1
		print 'seed 2 file is: ', seed_full_file_2
		plot_filename = time_dir+ '_' + tube_file + '.png'
		number_of_seeds_per_network = count_seeds_per_network(tube_full_file, seed_full_file_1, seed_full_file_2, plot_filename)
		print number_of_seeds_per_network
		

		network_size_time.extend(number_of_seeds_per_network)
		

	print "printing attached seed percentages"
	f1=open(time_dir+'_network_seeds.dat','w+')
	for network_seed in network_size_time:
		print >>f1, network_seed
	f1.close()




