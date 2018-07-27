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
import copy

class dtype():
    """A class that stores variable-type information."""
    
    descr = []
    names_ordered = []

    def clone(self):
        """Make a clone of this dtype.

            Returns:
                A DType, the clone.
        """

        cp = copy.deepcopy(self)
        cp.descr = copy.deepcopy(self.descr)
        cp.names_ordered = copy.deepcopy(self.names_ordered)
        return cp       

    def __init__(self, descr):
        # descr comes in as a list of (var, type) tuples

        # Temporarily convert descr to a dictionary
        sig_dict = {}
        for var, typ in self.descr:
            sig_dict[var] = typ
        
        # Add new elements
        self.names_ordered = []
        for var, typ in descr:
            typ = typ.replace("|", ""). replace("<", "")
            sig_dict[var] = typ
            self.names_ordered.append(var)
        
        # Convert back to list and store
        self.descr = [(k, sig_dict[k]) for k in self.names_ordered]

        #self.names = sig_dict.keys() #[i[0] for i in descr]
    
    def __str__(self):
        return str(self.descr)

    @staticmethod
    def convert(tp, val):
        """Convert a variable given a variable type.

            Args:
                tp  -- The variable type.
                val -- The variable to convert.

            Returns:
                The converted variable.
        """
        
        if tp == bool:
            return bool(val)

        if tp == int:
            return int(val)

        tp = tp.replace("|", "").replace("<", "")
        
        if tp[:1].upper() == "F":
            return float(val)
        
        if tp[:1].upper() == "I":
            return int(val)
        
        if tp[:1].upper() == "S" or tp[:1].upper() == "U":
            str_len = int(tp[1:])  # Get the number part
            val = str(val)[:str_len] # Keep only the correct few letters
            return val

    @property
    def names(self):
        """Get the names of the variable types.

            Returns:
                A list of strings, the variable types.
        """

        self.names_ordered =self.names_ordered[-len(self.descr):]
        return self.names_ordered
