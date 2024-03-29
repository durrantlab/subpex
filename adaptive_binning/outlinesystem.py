from __future__ import division, print_function

import logging

import numpy
from west.propagators import WESTPropagator
from west.systems import WESTSystem
from westpa.binning import FuncBinMapper, RectilinearBinMapper, RecursiveBinMapper

log = logging.getLogger("westpa.rc")

PI = numpy.pi
from numpy import *

pcoord_dtype = numpy.float32


#########
def function_map(coords, mask, output):
    # The coords are your incoming progress coordinates in array form and the output for return is an array of the same length containg the bin number in the corresponding indicie of the walker with coordinate in the coords array
    return output


class System(WESTSystem):
    def initialize(self):
        # You will need to specify actual values here
        self.pcoord_ndim = numberofdim
        self.pcoord_len = pcoordlength
        self.pcoord_dtype = numpy.float32
        self.bin_mapper = FuncBinMapper(function_map, numberofbins)
        self.bin_target_counts = numpy.empty((self.bin_mapper.nbins,), numpy.int_)
        self.bin_target_counts[...] = bintargetcount
