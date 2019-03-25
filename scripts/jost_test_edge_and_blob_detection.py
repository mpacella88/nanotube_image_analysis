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
from skimage import img_as_uint, img_as_bool, img_as_float
from scipy import ndimage
from scipy.spatial import distance
from scipy import ndimage as ndi
from numpy import unravel_index
import Tkinter, tkFileDialog
from skimage.external import tifffile
from skimage.feature import blob_dog, blob_log, blob_doh

def count_josh_circles_edge(filename):
	image= io.imread(filename, as_grey=True)
	image_color = io.imread(filename)
	print image 
	print image_color
	selem = disk(0)
	edges_open = canny(image, .001) #originally 2,1,25 last param can go up to 500 for improved performance, must lower for poorer images
			#edges_open = canny(image, 2) #originally 2,1,25
	
	edges = closing(edges_open, selem)
	fill_tubes = ndi.binary_fill_holes(edges)
		#io.imsave(str(i)+"_fill_tubes.png", img_as_uint(fill_tubes), cmap=cm.gray)
	#io.imsave("cha_8_fill_tubes.png", img_as_uint(fill_tubes))
	#io.imsave("cha_8_greyscale.png", image)
	io.imsave("fill_tubes.png", img_as_uint(fill_tubes) )
#count_josh_circles("Cha 8.tif")

def count_josh_circles_blob(filename):
	image_gray = io.imread(filename, as_grey=True)

	#this is optimal for larger cell images
	#blobs_log = blob_log(image_gray, min_sigma=5, max_sigma= 6, num_sigma=10, threshold=.01)

	#this is optimal for flow cell images
	blobs_log = blob_log(image_gray, min_sigma=1, max_sigma=, num_sigma=10, threshold=.01)

	print "number of seeds: "+str(len(blobs_log))

	# Compute radii in the 3rd column.
	blobs_log[:, 2] = blobs_log[:, 2] * math.sqrt(2)


	blobs_list = [blobs_log]
	colors = ['yellow']
	titles = ['Laplacian of Gaussian']
	sequence = zip(blobs_list, colors, titles)

	fig, ax = plt.subplots(figsize=(9, 3))
	ax.set_aspect('equal')
	#ax = axes.ravel()

	for idx, (blobs, color, title) in enumerate(sequence):
	    #ax[idx].set_title(title)
	    ax.imshow(image_gray, interpolation='nearest')
	    for blob in blobs:
	        y, x, r = blob
	        c = plt.Circle((x, y), r, color=color, linewidth=2, fill=False)
	        ax.add_patch(c)
	    #ax[idx].set_axis_off()

	plt.tight_layout()
	plt.show()

#count_josh_circles_blob("Cha 8.tif")
count_josh_circles_blob("1hr.10ml.Channel_10.tif")