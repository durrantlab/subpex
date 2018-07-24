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

"""Note that this code is NOT a replacement for numpy. It imitates numpy just
well enough to run some of the scoria functions that couldn't run
otherwise. Installing numpy/scipy is strongly recommended.
"""

from __future__ import absolute_import
import copy
from .DType import dtype as dtypeClass
from .Support import to_list
from .Support import var_type
from scoria.six.moves import range
from scoria.six.moves import zip
import scoria.six as six

#from .six import string_types as six_string_types 

def array(lst, dtype=""):
    """Determines whether or not a 1D or 2D array should be used.

        Args:
            lst -- A list to convert to an array.

        Returns:
            An Array2D or Array1D object, as required.
    """

    # If it's a 1D array with non-list elements, just return it, but
    # converting the types.
    try:
        if var_type(lst) == "1D":
            if var_type(lst[0]) in ["number", "string"]:
                if dtype == "":
                    return lst
                else:
                    return lst_cp.astype(dtype)
    except: pass

    # check if it's a list of 1D arrays. If so, convert to list of lists.
    try: 
        if var_type(lst[0]) == "1D":
            lst = [l.lst for l in lst]
    except: pass

    if var_type(lst) == "number":
        return lst
    elif len(lst) == 0:
        # It's a 1D array
        return Array1D(lst, dtype)
    elif type(lst[0]) is list or type(lst[0]) is tuple:
        # It's a list of lists (2D arraay)
        return Array2D(lst, dtype)
    else:
        return Array1D(lst, dtype)
    
    raise ValueError('Could not create array. Try installing NUMPY/SCIPY.')
    

class ArrayParent:
    """The parent of all Array classes."""
    
    lst = []
    shape = ()
    type = ""

    def __iter__(self):
        # So you could use for... in...
        for x in self.lst:
            yield x

    def __str__(self):
        return str(self.lst)

    def __repr__(self):
        return "array(" + self.__str__() + ")"

    def __getitem__(self, selection):
        # []
        new_lst = []

        # If selection is itself an array, convert it to a list.
        selection = to_list(selection)

        try:
            # Assume the selection is iterable.
            for i in selection:
                new_lst.append(self.lst[i])
        except:
            # selection must not be iterable. Just a number.
            new_lst.append(self.lst[selection])
    
        # If it's just one item, then just return the value, not the value in
        # an array
        if len(new_lst) == 1:
            return new_lst[0]

        new_arr = array(new_lst)
        new_arr.set_shape()

        return new_arr
    
    def __setitem__(self, key, item):       
        if var_type(key) in ["string", "number"]:
            # The key is not iterable. Probably just a number.
            self.lst[key] = item
            return
        else:
            # The key is iterable
            if var_type(item) in ["string", "number"]:
                # But the item is not iterable
                for k in key:
                    self.lst[k] = item
                return
            else:
                # So the item is also iterable.
                for k, i in zip(key, item):
                    self.lst[k] = i
                return
        
        raise ValueError('Could not set array item. Try installing NUMPY/SCIPY.')

    def __len__(self):
        return len(self.lst)

    def copy(self):
        return copy.deepcopy(self)

class Array1D(ArrayParent):
    """A 1D Array."""

    type = "1D"

    def __init__(self, lst, dtype=""):
        if dtype != "":
            for i, val in enumerate(lst):
                lst[i] = dtypeClass.convert(dtype, val)

        self.lst = lst

        self.set_shape()
    
    def set_shape(self):
        """Sets the shape of this array."""

        self.shape = (len(self.lst),)

    def __eq__(self, other):
        bools = copy.deepcopy(self.lst)  # Just to hold the booleans.
        
        for x in range(self.shape[0]):
            if var_type(other) in ["string", "number"]:
                bools[x] = (bools[x] == other)
            else:
                bools[x] = (bools[x] == other[x])
            
        return array(bools)

    def astype(self, dtype):
        """Casts this array as a given type.

            Args:
                dtype -- The type to cast.

            Returns:
                The new array.
        """
    
        cp = copy.deepcopy(self)
        for i, val in enumerate(cp.lst):
            cp.lst[i] = dtypeClass.convert(dtype, val)
        return cp

    def __mul__(self, other):
        # pairwise multiplication

        # What is other?
        typ = var_type(other)

        if typ == "1D":
            new_lst = []

            return array([l[0] * l[1] for l in zip(self.lst, other.lst)])

        elif typ == "number":
            new_lst = [other * l for l in self.lst]
            return array(new_lst)
        
        raise ValueError('Could not perform multiplication. Try installing NUMPY/SCIPY.')
    
    def __rmul__(self, other):
        return self.__mul__(other)

    def __div__(self, other):
        typ = var_type(other)

        if typ == "number":
            new_lst = [l / other for l in self.lst]
            return array(new_lst)
        
        raise ValueError('Could not perform division. Try installing NUMPY/SCIPY.')

    def __add__(self, other):
        typ = var_type(other)
        if typ == "1D":
            new_lst = []
            for t in zip(self.lst, other.lst):
                new_lst.append(sum(t))
            
            return array(new_lst)
        elif typ == "number":
            new_lst = [t + other for t in self.lst]
            return array(new_lst)     
        
        raise ValueError('Could not perform addition. Try installing NUMPY/SCIPY.')

    def __sub__(self, other):
        other = -1 * other
        return self.__add__(other)   
    
class Array2D(ArrayParent):
    """A 2D Array."""

    type = "2D"

    def __init__(self, lst, dtype="float"):
        self.lst = []
        for row in lst:
            self.lst.append(array(row))
        
        self.set_shape()
    
    def set_shape(self):
        """Sets the shape of this array."""

        self.shape = (len(self.lst), len(self.lst[0]))

    def __eq__(self, other):
        bools = copy.deepcopy(self.lst)
        for x in range(self.shape[0]):
            for y in range(self.shape[1]):
                bools[x][y] = (bools[x][y] == other)
        return array(bools)
    
    @property
    def T(self):
        """Returns the transpose of this array.

            Returns:
                An Array2D, the transpose.
        """
        
        coors = []
        for coor in zip(*self.lst):
            coors.append(array(list(coor)))
        return array(coors)
    
    def __mul__(self, other):
        # pairwise multiplication

        # What is other?
        typ = var_type(other)

        if typ == "2D":
            new_lst = []
            for t in zip(self.lst, other.lst):
                new_ln = []
                for i in range(len(t[0])):
                    new_ln.append(t[0][i] * t[1][i]) 
                new_ln = array(new_ln)
                new_lst.append(new_ln)
            
            return array(new_lst)
        elif typ == "1D":
            new_lst = []
            # So assume the dimensions of 1D are the same as the elements of this 2D array.
            for l in self.lst:
                l = array(l)
                other = array(other)
                new_lst.append(l * other)
            return array(new_lst)
        elif typ == "number":
            new_lst = []
            for row in self.lst:
                new_lst.append(copy.deepcopy(row) * other)
            return array(new_lst)
        
        raise ValueError('Could not perform multiplication. Try installing NUMPY/SCIPY.')
    
    def __add__(self, other):
        typ = var_type(other)
        if typ == "number":
            new_lst = []
            for i in range(len(self.lst)):
                new_lst.append([])
                for j in range(len(self.lst[i])):
                    new_lst[i].append(self.lst[i][j] + other)
            return array(new_lst)
        elif typ == "2D":
            new_lst = []
            for i in range(len(self.lst)):

                if i >= len(other.lst):
                    toadd = other.lst[-1]
                else:
                    toadd = other.lst[i]

                new_lst.append(self.lst[i] + toadd)
            return array(new_lst)
        elif typ == "1D":
            # assume 1d has same dimensions as element of 2D
            new_lst = [l + other for l in self.lst]
            return array(new_lst)
        
        raise ValueError('Could not perform addition. Try installing NUMPY/SCIPY.')

    def __sub__(self, other):
        other = other * -1
        return self.__add__(other)
    
class RecArray:
    """A record array."""

    type = "Rec"
    dtype = dtypeClass([])

    # Collection of arrays accessible through dictionary keys.
    dict = {}
    ndim = 1  # I think for a RecArray this is always 1?

    def __init__(self, dict, dtypes = None):
        # dict maps keys to arrays.
        for key in dict:
            self.dict[key] = array(dict[key])
        
        if dtypes is not None:
            dtypes = dtypeClass(dtypes)
            self.dtype = None
            self.dtype = dtypes
            self.astype(dtypes)
        
    def __len__(self):
        key = list(self.dict.keys())[0]
        return len(self.dict[key])

    def __repr__(self):
        return "recarray(" + self.__str__() + ")"
    
    def __str__(self):
        return str(self.dict)

    def __getitem__(self, lookup_key):
        var_typ = var_type(lookup_key)

        # If it's just an int, convert it to an array
        if var_typ == "number":
            lookup_key = [lookup_key]

        if var_typ == "string":
            # It's a string, so just lookup the corresponding column
            return self.dict[lookup_key]
        elif isinstance(lookup_key[0], six.string_types):
            # So it's a list of strings. Return the columns.

            # So key must be something that should act on the individual
            # Arrays, like ["key1", "key2", "key3", "key4""]
            new_dict = {}
            for str_key in lookup_key:
                if str_key in list(self.dict.keys()):
                    new_dict[str_key] = self.dict[str_key]
            
            #print "DDD", new_dict.keys()

            updated_dict = RecArray({})
            updated_dict.dict = new_dict

            # get the new descr
            descr = [d for d in self.dtype.descr if d[0] in list(new_dict.keys())]
            names_ordered = [d for d in self.dtype.names_ordered if d in list(new_dict.keys())]

            updated_dict.dtype = dtypeClass(descr)
            updated_dict.dtype.names_ordered = names_ordered

            return updated_dict

        else: # So it's a list of integers
            # So key must be something that should act on the individual
            # Arrays, like [1, 2, 3, 4]
            cp_dict = self.dict.copy()
            new_dict = {}
            for str_key in cp_dict.keys():
                new_col = cp_dict[str_key][lookup_key]
                new_dict[str_key] = new_col

            new_array = self.clone()
            new_array.dict = new_dict
            new_array.dtype = self.dtype

            return new_array
    
    def clone(self):
        """Makes a clone of this array.

            Returns:
                A RecArray, the clone.
        """
        
        cp = copy.deepcopy(self)
        cp.dict = self.dict.copy()
        cp.dtype = self.dtype.clone()
        cp.ndim = self.ndim
        cp.type = self.type
        return cp
    
    def __setitem__(self, key, vals):
        # In this implementation, key must be a string.
        var_typ = var_type(vals)

        if var_typ in ["list", "1D", "2D", "Rec"]:
            self.dict[key] = vals
        else:  # All vals in this col being set to same value.
            for i in range(len(self.dict[key])):
                self.dict[key][i] = vals

    def astype(self, dtype):
        """Cast this array as a certain type.

            Returns:
                The new array.
        """

        cp = copy.deepcopy(self)

        cp.dtype = dtype
        for key, tp in dtype.descr:
            for i, val in enumerate(cp.dict[key]):
                cp.dict[key][i] = dtypeClass.convert(tp, val)
        return cp
    
    def copy(self):
        """Makes a copy of this array.

            Returns:
                A RecArray, the copy.
        """

        return copy.deepcopy(self)

