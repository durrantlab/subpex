import h5py
import numpy
import numpy.random
import matplotlib.pyplot as plt

num_bins = 10

# Generate some test data
#x = numpy.random.randn(8873)
#y = numpy.random.randn(8873)
f = h5py.File("west.h5", "r+")

coors = numpy.vstack([f["iterations"][k]["pcoord"].value[:,1] for k in f["iterations"].keys()])

# Remove ones that are just 0s
coors = coors[numpy.nonzero(numpy.logical_and(coors[:,1] != 0, coors[:,0] != 0))[0]]

x = coors[:,0]
y = coors[:,1]

heatmap, xedges, yedges = numpy.histogram2d(x, y, bins=num_bins)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

plt.clf()
plt.imshow(heatmap.T, extent=extent, interpolation="gaussian", origin='lower')
plt.show()

