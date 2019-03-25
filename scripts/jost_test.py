import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import math
from matplotlib import cm
from matplotlib import path
from skimage import img_as_uint, img_as_bool, img_as_float
import Tkinter, tkFileDialog
from skimage.external import tifffile
from skimage import io 
from skimage.feature import blob_dog, blob_log, blob_doh



def count_josh_circles_blob(filename):
	image_gray = io.imread(filename, as_grey=True)

	#this is optimal for larger cell images
	blobs_log = blob_log(image_gray, min_sigma=5, max_sigma= 6, num_sigma=10, threshold=.01)

	#this is optimal for flow cell images
	#blobs_log = blob_log(image_gray, min_sigma=1, max_sigma=2, num_sigma=10, threshold=.01)

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
	#plt.show()
	plt.savefig("processed_image.tif")

count_josh_circles_blob("Cha 8.tif")
#count_josh_circles_blob("1hr.10ml.Channel_10.tif")