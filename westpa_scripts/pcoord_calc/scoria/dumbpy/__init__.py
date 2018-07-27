# Copyright 2017 Jacob D. Durrant
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Loads external modules (numpy, scipy) if available. Otherwise, uses cheap
imitations."""

from __future__ import absolute_import
from __future__ import print_function
import textwrap
import sys

# A list of the available dependencies.
dependencies_available = []

# A flag for supressing errors
missing_dependency_throws_error = True

# Try to load numpy.
force_dumbpy = False  # True for debugging.

if '__pypy__' in sys.builtin_module_names:
    # It's pypy, so don't load numpy
    force_dumbpy = True
if len(sys.argv) > 1 and sys.argv[1].upper() == "NODEPENS":
    # The command-line parameter indicates you shouldn't use dependencies
    force_dumbpy = True

# Python3 requires soe minor tweaks
python_version = sys.version_info[0]

try:
    # Try to load traditional numpy
    if force_dumbpy: raise ValueError('Using dumbpy')

    from numpy import array
    from numpy.core.defchararray import strip as defchararray_strip 
    from numpy.core.defchararray import add as defchararray_add
    from numpy.core.defchararray import rjust as defchararray_rjust
    from numpy.core.defchararray import upper as defchararray_upper
    from numpy.core.defchararray import lstrip as defchararray_lstrip
    from numpy.core.defchararray import split as defchararray_split
    from numpy.lib.recfunctions import append_fields
    from numpy import genfromtxt
    from numpy import nonzero
    from numpy import logical_or
    from numpy import core
    from numpy import dtype
    from numpy import vstack
    from numpy import zeros
    from numpy import max
    from numpy import unique
    from numpy import insert
    from numpy import logical_not
    from numpy import logical_and
    from numpy import append
    from numpy import arange
    from numpy import savez
    from numpy import load
    from numpy import ones
    from numpy import empty
    from numpy import sum
    from numpy import min
    from numpy import delete
    from numpy import setdiff1d
    from numpy import hstack
    from numpy import intersect1d
    from numpy import setxor1d
    from numpy import power
    from numpy import cos
    from numpy import sin
    from numpy import sqrt
    from numpy import dot
    from numpy import linalg
    from numpy import amin
    from numpy import lib
    from numpy import identity
    from numpy import mean
    from numpy import transpose
    from numpy import argmax
    from numpy import ma
    from numpy import arccos
    from numpy import cross
    from numpy import arctan2
    from numpy import fabs
    from numpy import ones
    from numpy import empty
    from numpy import sum
    from numpy import delete
    from numpy import setdiff1d
    from numpy import power
    from numpy import cos
    from numpy import sin
    from numpy import sqrt
    from numpy import identity
    from numpy import mean
    from numpy.linalg import norm
    from numpy import fabs
    from numpy.lib.recfunctions import stack_arrays
    from numpy import std  # Note that there is no dumbpy equivalent yet.
    from numpy import ndarray
    from numpy import degrees
    from numpy import radians

    def get_col(arr, num):
        return arr[:, num]

    dependencies_available.append("NUMPY")
except:
    # Numpy not available, so load substitute (dumbpy)
    from .Array import array
    from .Utils import genfromtxt
    from .Utils import nonzero
    from .Utils import logical_or
    from .Utils import logical_not
    from .Utils import logical_and
    from .Utils import defchararray_strip
    from .Utils import defchararray_add
    from .Utils import defchararray_rjust
    from .Utils import defchararray_upper
    from .Utils import defchararray_lstrip
    from .Utils import defchararray_split
    from .DType import dtype
    from .Utils import vstack
    from .Utils import append_fields
    from .Utils import zeros
    from .Utils import _max as max
    from .Utils import _min as min
    from .Utils import unique
    from .Utils import insert
    from .Utils import append
    from .Utils import arange
    from .Utils import get_col
    from .Utils import ones
    from .Utils import empty
    from .Utils import sum
    from .Utils import delete
    from .Utils import setdiff1d
    from .Utils import power
    from .Utils import _cos as cos
    from .Utils import _sin as sin
    from .Utils import sqrt
    from .Utils import identity
    from .Utils import mean
    from .Utils import norm
    from .Utils import fabs
    from .Utils import stack_arrays
    from .Utils import transpose

    dependencies_available.append("DUMBPY")

try:
    # Try to load traditional scipy
    if force_dumbpy: raise ValueError('Using dumbpy')

    from scipy.spatial.distance import squareform
    from scipy.spatial.distance import pdist
    from scipy.spatial.distance import cdist
    dependencies_available.append("SCIPY")
except:
    # Use cheap replacement instead. None of these functions are available
    # through dumbpy.
    pass


def class_dependency(action, dependency):
    """Determines whether or not a given dependency is available.

        Args:
            action     -- A string describing the action you'd like to
                          perform.
            dependency -- A string, the dependency required for that action.

        Returns:
            A boolean, true if the dependency is available. Prints a message
                otherwise.
    """

    global dependencies_available
    global missing_dependency_throws_error

    if not dependency in dependencies_available:
        t = textwrap.TextWrapper(width = 60, subsequent_indent="       ")
        print()
        print(("\n".join(t.wrap("Error: Cannot " + action +
                               ". Install this recommended dependency: " +
                               dependency))))
        print()
        if missing_dependency_throws_error:
            raise ImportError("The " + dependency + " module is not available.")
            return False
        else:
            return False
    return True
