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

from __future__ import absolute_import
from .Array import RecArray
from .Array import array
from .DType import dtype as dtypeClass
from .Support import to_list
from .Support import var_type
import math
from ..six.moves import range
from ..six.moves import zip

def genfromtxt(fname, dtype = "", names = [], delimiter = []):
    """Generates an array from a text file.

        Args:
            fname     -- A string, the file name.
            dtype     -- The variable types.
            names     -- A list of string, the names of the variables.
            delimiter -- A list of numbers, the legths of each field.

        Returns:
            A RecArray object.
    """

    # fname here is a file object.

    # Load the data
    data = {}
    for n in names:
        data[n] = []

    for line in fname:
        parts = []
        for num in delimiter:
            parts.append(line[:num])
            line = line[num:]
        
        for i, name in enumerate(names):
            data[name].append(parts[i])

    # Fix any types
    dtype = [l.strip() for l in dtype.split(",")]
    for i, name in enumerate(names):
        if "int" in dtype[i] or dtype[i][:1] == "i":
            data[name] = [int(l) for l in data[name]]
            dtype[i] = "int"
        if "float" in dtype[i] or dtype[i][:1] == "f":
            data[name] = [float(l) for l in data[name]]
            dtype[i] = "float"
    
    dtypes = list(zip(names, dtype))

    return RecArray(data, dtypes)

def nonzero(arr):
    """Identifies which entries in an array are not zero.

        Args:
            arr  -- The array.

        Returns:
            A list of lists (the indecies along each dimension).
    """

    if len(arr.shape) == 1:
        indx_to_keep = []
        for i, val in enumerate(arr):
            if val != 0:
                indx_to_keep.append(i)
        indx_to_keep = (array(indx_to_keep),)
    elif len(arr.shape) == 2:
        indx1_to_keep = []
        indx2_to_keep = []

        for i, row in enumerate(arr):
            for j, val in enumerate(row):
                if val != 0:
                    indx1_to_keep.append(i)
                    indx2_to_keep.append(j)
        indx_to_keep = (array(indx1_to_keep), array(indx2_to_keep))

    return indx_to_keep

def logical_or(arr1, arr2):
    """Applies the logical or element wise.

        Args:
            arr1 -- The first array.
            arr2 -- The second array.

        Returns:
            An boolean array.
    """

    if var_type(arr1) == "list":
        arr1 = array(arr1)

    if len(arr1.shape) == 1:
        or_result = [x or y for x,y in zip(arr1, arr2)]
        return array(or_result)

def logical_and(arr1, arr2):
    """Applies the logical and element wise.

        Args:
            arr1 -- The first array.
            arr2 -- The second array.

        Returns:
            An boolean array.
    """

    if len(arr1.shape) == 1:
        and_result = [x and y for x,y in zip(arr1, arr2)]
        return array(and_result)

def logical_not(arr):
    """Applies the logical not element wise.

        Args:
            arr -- The first array.

        Returns:
            An boolean array.
    """

    # Only for 1D arrays
    return array([not i for i in arr])

def defchararray_strip(arr):
    """Strips spaces from the strings in a string array.

        Args:
            arr -- The string array.

        Returns:
            A string array.
    """

    if len(arr.shape) == 1:
        strip_result = [s.strip() for s in arr]
        return array(strip_result)

def defchararray_add(arr, addit):
    """Adds a string to each element in a strig array.

        Args:
            arr   -- The string array.
            addit -- The string to add.

        Returns:
            A string array.
    """

    # If addit is a string
    if var_type(addit) == "string":
        lst = [l + addit for l in arr.lst]
        return array(lst)
    else: # addit is another array of some sort
        added = [l[0] + l[1] for l in zip(arr, addit)]
        return array(added)
    raise ValueError('Could not add strings. Try installing NUMPY/SCIPY.')

def defchararray_rjust(arr, width):
    """Right justifies the strings in a string array.

        Args:
            arr   -- The string array.
            width -- The width of the new elements.

        Returns:
            A string array.
    """

    arr = to_list(arr)
    return array([a.rjust(width) for a in arr])

def defchararray_upper(arr):
    """Make the strings in a string array uppercase.

        Args:
            arr -- The string array.

        Returns:
            A string array, all upper case.
    """

    arr = to_list(arr)
    return array([a.upper() for a in arr])

def defchararray_lstrip(arr, chars = [" ", "\t"]):
    """Strips left spaces from the strings in a string array.

        Args:
            arr -- The string array.
            chars -- A list of the characters to consider white space.

        Returns:
            A string array.
    """

    arr = to_list(arr)
    return array([a.lstrip(chars) for a in arr])

def defchararray_split(arr, num = -1):
    """Splits the strings in a string array.

        Args:
            arr -- The string array.
            num -- The number of splits to make. Defaults to -1 (all).

        Returns:
            An array of string arrays.
    """

    arr = to_list(arr)
    return array([a.split(num) for a in arr])

def vstack(arrays):
    """Stacks arrays.

        Args:
            arrays -- A list of arrays.

        Returns:
            A string array.
    """

    if arrays[0].type == "1D":
        lsts = [a.lst for a in arrays]
        return array(lsts)
    elif arrays[0].type == "2D":
        new_array = to_list(arrays[0]) + to_list(arrays[1])
        return array(new_array)
    
    raise ValueError('Could not perform vstack. Try installing NUMPY/SCIPY.')

def append_fields(arr, field_name, data, usemask=False):
    """Append columns to a rec array.

        Args:
            arr        -- The rec array.
            field_name -- A string, the new field name.
            data       -- The data to add under that field name.

        Returns:
            A rec array.
    """

    arr[field_name] = data
    arr.dtype.names_ordered.append(field_name)
    arr.dtype.descr.append((field_name,))
    return arr

def extrema(func, arr, axis = 0):
    """Calculates the extrema (max or min) of an array.

        Args:
            func -- The function to apply, max or min.
            arr  -- The array.
            axis -- The axis (0 by default).

        Returns:
            The extrema value across the axis.
    """

    if axis < 0:
        axis = len(arr.shape) + axis

    var_typ = var_type(arr)

    if var_typ == "list":
        arr = array(arr)
        var_typ = var_type(arr)

    if var_typ == "1D":
        return func(arr.lst)
    elif axis == 0:
        new_list = []
        for line in arr.T:
            new_list.append(func(line))
        new_list = array(new_list)
        return new_list
    
    raise ValueError('Could not identify extrema (min or max). Try installing NUMPY/SCIPY.')  # An error here if it ever tries on axis 1!

def _max(arr, axis = 0):
    """The max function.

        Args:
            arr  -- The array.
            axis -- The axis (0 by default)

        Returns:
            The max value along the axis.
    """

    return extrema(max, arr, axis)

def _min(arr, axis = 0):
    """The min function.

        Args:
            arr  -- The array.
            axis -- The axis (0 by default)

        Returns:
            The min value along the axis.
    """

    return extrema(min, arr, axis)

def unique(arr):
    """Return the unique values in an array.

        Args:
            arr  -- The source array.

        Returns:
            An array, with unique values.
    """

    # This has only been implemented fo 1d arrays
    if var_type(arr) == "string":
        arr = array([arr])
    return array(list(set(arr.lst)))

def insert(arr, indx, val):
    """Insert a value into an array.

        Args:
            arr  -- The source array.
            indx -- The index at which to insert the value.
            val  -- The value to insert.

        Returns:
            An array, with the value inserted.
    """

    # A very limited implementation of insert.
    lst = arr.lst[:]
    lst.insert(indx, val)
    return array(lst)  

def append(arr1, to_append):
    """Add a value to the end of an array.

        Args:
            arr       -- The source array.
            to_append -- The value to append.

        Returns:
            An array, with the value appended.
    """

    lst = arr1.lst[:]

    # First, assume to_append is a list.
    if var_type(to_append) == "list":
        lst.extend(to_append)
    else:
        # It's not a list.
        lst.append(to_append)

    return array(lst)

def arange(start, stop, step, dtype = "f8"):  # float by default
    """An array with values that begin at start and end at stop, spaced step 
    apart.

        Args:
            start -- The starting value.
            stop  -- The stopping value.
            step  -- The distance between values.
            dtype -- The variable type (string).

        Returns:
            An array, with the specified equidistant values.
    """

    return [
        dtypeClass.convert(dtype, t * step) 
        for t in range(int(start / step), int(stop / step), 1)
    ]

def get_col(lst, num):
    """Return the column of a 2D array.

        Args:
            lst -- The array.
            num -- The column index.

        Returns:
            The specified column.
    """

    return array([a[num] for a in lst])

def all_same_num(dims, num, dtype="float"):
    """Make a square array filled with the same values.

        Args:
            dims  -- An int, the dimension of the square array.
            num   -- The number to fill.
            dtype -- A string, the variable type.

        Returns:
            An array.
    """

    arr = []
    val = dtypeClass.convert(dtype, num)

    # Start by assuming dims is a tuple or list
    try:
        for xi in range(dims[0]):
            arr.append([val] * dims[1])
        new_array = array(arr)
    except:
        # So it's an int.
        new_array = array([val] * dims)
    
    return new_array

def ones(dims, dtype="float"):
    """Make a square array filled with ones.

        Args:
            dims  -- An int, the dimension of the square array.
            dtype -- A string, the variable type.

        Returns:
            An array.
    """

    return all_same_num(dims, 1.0, dtype=dtype)

def zeros(dims, dtype="float"):
    """Make a square array filled with zeros.

        Args:
            dims  -- An int, the dimension of the square array.
            dtype -- A string, the variable type.

        Returns:
            An array.
    """

    return all_same_num(dims, 0.0, dtype=dtype)

def empty(dims, dtype="float"):
    """Make a square array filled with zeros.

        Args:
            dims  -- An int, the dimension of the square array.
            dtype -- A string, the variable type.

        Returns:
            An array.
    """

    return zeros(dims, dtype=dtype)

def sum(arr, axis = 0):
    """Sum the values of an array along an axis.

        Args:
            arr  -- The source array.
            axis -- The axis. Defaults to 0.

        Returns:
            The sum along the axis.
    """

    if axis < 0:
        axis = len(arr.shape) + axis

    if var_type(arr) == "1D":
        return sum(arr.lst)
    elif axis == 0:
        sumit = arr[0]
        
        try:
            remainder = arr.lst[1:]
        except:
            remainder = arr[1:] 

        for a in remainder:
            sumit = sumit + a

        return array(sumit)
    elif axis == 1:
        return sum(arr.T)

    raise ValueError('Could not perform sum. Try installing NUMPY/SCIPY.')

def delete(arr, indx_to_delete):
    """Delete the values in an array.

        Args:
            arr  -- The source array.
            axis -- The indecies to delete.

        Returns:
            The new array, with elements removed.
    """

    # make sure arr is an array
    try:
        arr.type
    except:
        arr = array(arr)

    if var_type(arr) == "1D":
        # This only works with 1D arrays
        arr_lst = to_list(arr)
        all_indx = list(range(len(arr_lst)))
        to_keep = array(list(set(all_indx) - set(indx_to_delete)))

        arr = array(arr)

        return arr[to_keep]

    raise ValueError('Could not delete array elements. Try installing NUMPY/SCIPY.')

def setdiff1d(arr1, arr2):
    """Get the difference between two arrays.

        Args:
            arr1 -- The first array.
            arr2 -- The second array.

        Returns:
            The difference between the arrays.
    """

    arr1_lst = set(to_list(arr1))
    arr2_lst = set(to_list(arr2))
    lst = list(arr1_lst - arr2_lst)
    return array(lst)

def power(num, p):
    """Calculate the power of a number.

        Args:
            num -- The number.
            p   -- The exponent.

        Returns:
            The number num ^ p.
    """

    if var_type(num) == "number":
        # This only works with 1D arrays
        # Just implemented for scalar at this point
        return math.pow(num, p)
    
    raise ValueError('Could not raise to the power. Try installing NUMPY/SCIPY.')

def _sin(num):
    """The sin function.

        Args:
            num -- A number.

        Returns:
            sin(num)
    """

    if var_type(num) == "number":
        # This only works with 1D arrays
        # Just implemented for scalar at this point
        return math.sin(num)
    
    raise ValueError('Could not calculate sin(). Try installing NUMPY/SCIPY.')

def _cos(num):
    """The cos function.

        Args:
            num -- A number.

        Returns:
            cos(num)
    """

    if var_type(num) == "number":
        # This only works with 1D arrays
        # Just implemented for scalar at this point
        return math.cos(num)
    
    raise ValueError('Could not calculate cos(). Try installing NUMPY/SCIPY.')

def sqrt(num):
    """The square root function.

        Args:
            num -- A number.

        Returns:
            sqrt(num)
    """

    if var_type(num) == "number":
        # This only works with 1D arrays
        # num must be a scalar for now.
        return math.sqrt(num)
    
    raise ValueError('Could not calculate square root. Try installing NUMPY/SCIPY.')

def identity(dimen):
    """A 2D identity array.

        Args:
            dimen -- The dimension of the square array.

        Returns:
            An identity array.
    """

    arr = zeros((dimen, dimen))
    for t in range(dimen):
        arr[t][t] = 1.0
    return arr

def mean(arr, axis = 0):
    """Calculates the mean of an array along an axis.

        Args:
            arr  -- The source array.
            axis -- The axis. Defaults to 0.

        Returns:
            The mean along the axis.
    """

    if axis < 0:
        axis = len(arr.shape) + axis

    s = sum(arr, axis=axis)
    mean = s * (1.0 / s.shape[axis])
    return mean

def norm(arr):
    """The length of a 1D vector.

        Args:
            arr -- The 1D vector.

        Returns:
            A number, the vector's length.
    """

    return sqrt(sum(arr * arr))

def fabs(num):
    """The aboslute value function.

        Args:
            num -- A number.

        Returns:
            fabs(num)
    """

    if var_type(num) == "number":
        return math.fabs(num)
    
    raise ValueError('Could not determine absolute value. Try installing NUMPY/SCIPY.')

def stack_arrays(arr_list, usemask = False):
    """Like vstack, but for rec arrays.

        Args:
            arr_list -- A list of rec arrays.
            usermask -- Not sure. Always false, not used.

        Returns:
            A rec array.
    """

    # arr1 and arr2 are Rec arrays! Known as recarrays in numpy.
    # Separate out the arrays
    arr1 = arr_list[0]
    arr2 = arr_list[1]
    
    # Identify the keys they have in common.
    keys1 = set(arr1.dict.keys())
    keys2 = set(arr2.dict.keys())
    common_keys = keys1 & keys2

    dict = {}
    for key in common_keys:
        dict[key] = to_list(arr1[key]) + to_list(arr2[key])

    return RecArray(dict)

def transpose(arr):
    """Transpose an array

        Args:
            arr -- The array to transpose.

        Returns:
            The transposed array.
    """
    
    if var_type(arr) == "list":
        arr = array(arr)

    return arr.T