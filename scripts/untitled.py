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


#quick script to create sub-plots of histograms for each timepoint network growth time series data 
data_2 = open("2_network_seeds.dat")
data_2_network_sizes=data_2.readlines()
print data_2_network_sizes