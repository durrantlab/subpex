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
import scoria.six as six
#from .six import string_types as six_string_types


def var_type(var):
    """A helper function to identify a variable's type.

        Args:
            var  -- The variable.

        Returns:
            A string, the variable type.
    """

    if isinstance(var, six.string_types):
        return "string"
    
    # See if it's one of your dumbpy array objects
    try: return var.type
    except: pass

    # see if it's a list
    if type(var) is list:
        return "list"
    else:
        # a number perhaps?
        return "number"

def to_list(arr):
    """Convert an array to a list.

        Args:
            arr  -- A 1D or 2D array.

        Returns:
            The list.
    """

    # return a list regardless of whether arr is a list or array.
    typ = var_type(arr) 
    if typ in ["1D", "2D"]:
        return arr.lst[:]
    elif typ == "list":
        # already a python list
        return arr[:]
    else:
        # a number or string
        return arr
