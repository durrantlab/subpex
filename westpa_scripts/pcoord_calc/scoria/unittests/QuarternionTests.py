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
import unittest
import os
import sys

#import numpy as np
from scoria import dumbpy as np

import scoria


class OtherMoleculeTests(unittest.TestCase):
    """
    Base Test Suite
    """

    # Initialization and destruction for each test.

    def setUp(self):
        """
        Setting up the test molecule.
        """
        self.mol = scoria.Molecule()

    def tearDown(self):
        """
        Cleans up variables for the next test.
        """
        self.mol = None

        # Remove all files from the tmp folder?

    ### Tests

    @unittest.skip("Needs test written")
    def test_realignment(self):
        """
        Empty test.
        """
        
