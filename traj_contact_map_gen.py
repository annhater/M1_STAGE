import numpy
import numpy.linalg
import MDAnalysis
import MDAnalysis.core.distances
from numpy import *
#------------------------------------------------------------------------------------#
# Contact Analysis between two selections,
# 0 to 1 (0=never below threshold, 1=always above threshold)
# by M.E. Pitman
#------------------------------------------------------------------------------------#

#load trajectory
u = MDAnalysis.Universe("simulations_1HSI_APO_MP/V7/V7.gro", "simulations_1HSI_APO_MP/V7/md_V7.xtc")
timestep = 1
n_frames = len(u.trajectory)-1
start = 1

#select the regions to be measured against eachother from the trajecotry
chain1 = u.select_atoms('* and name CA')
chain2 = u.select_atoms('* and name CA')

#determine dimensions of matrix
n1 = len(chain1)
n2 = len(chain2)

#initialize array with zeros
contact_sum = numpy.zeros((n1, n2))

#define your max_distance (cutoff, example 5.0)
max_distance = 5.0

#sum 1 or 0 of each matrix position across timesteps
for ts in u.trajectory[start::timestep]:
	ch1 = chain1.positions
	ch2 = chain2.positions
	ts_dist = distance.cdist(ch1, ch2, 'euclidean')
	ts_dist[ts_dist < max_distance] = 1
	ts_dist[ts_dist > max_distance] = 0
	contact_sum = ts_dist + contact_sum	
						
contact_ratio = contact_sum/n_frames
#---------------------------------------------------------------------------#
#basic imaging
#---------------------------------------------------------------------------#
from pylab import imshow, xlabel, ylabel, xlim, ylim, colorbar, cm, clf
import matplotlib.pyplot as plt

class Formatter(object):
    def __init__(self, im):
        self.im = im
    def __call__(self, x, y):
        z = self.im.get_array()[int(y), int(x)]
        return 'x={:.01f}, y={:.01f}, z={:.01f}'.format(x, y, z)

#set x_min and y_min to the lowest residue index (example residue 50)
cr_shape = contact_ratio.shape
x_shift = cr_shape[1]
y_shift = cr_shape[0]
x_min = 1
y_min = 1
x_max = x_min + x_shift
y_max = y_min + y_shift

#had aspect= equal
im = plt.imshow(contact_ratio, vmin=0, vmax=1, aspect='equal', origin='lower', extent=[x_min,x_max, y_min, y_max] )

im.set_cmap('hot')
plt.grid(b=True, color='#737373')

im.set_interpolation('nearest')
plt.format_coord = Formatter(im)
delta = 1

plt.xticks( fontsize = 20)
plt.yticks( fontsize = 20)

colorbar()
plt.show()