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
import warnings
import unittest
from . import InformationTests as IT
from . import FileIOTests as FIOT
from . import GeometryTests as GT
from . import ManipulationTests as MT
from . import OtherMoleculeTests as OMT
from . import SelectionTests as ST

warnings.filterwarnings("ignore", category=DeprecationWarning)

class UnitTests(object):
    """
    Unit testing object for scoria.
    """
    def __init__(self):
        """
        Initalizes the unit tests.
        """
        self._suite = unittest.TestSuite()
        self._runner = unittest.TextTestRunner()

    # Running Suite

    def run(self):
        """
        Runs the currently queued suite of tests.
        """
        self._runner.run(self._suite)

    def run_all(self):
        """
        Quickly runs all unit tests.
        """
        self.add_all_tests()
        self.run()

    # Add specific module tests

    def add_all_tests(self):
        """
        Adds all available tests to the suite.
        """
        self.add_information_tests()
        self.add_fileio_tests()
        self.add_geometry_tests()
        self.add_selection_tests()
        self.add_manipulation_tests()
        self.add_other_molecule_tests()

    def add_information_tests(self):
        """
        Adds the information tests.
        """
        information_tests = unittest.makeSuite(IT.InformationTests)
        self._suite.addTests(information_tests)

    def add_fileio_tests(self):
        """
        Adds the information tests.
        """
        fileio_tests = unittest.makeSuite(FIOT.FileIOTests)
        self._suite.addTests(fileio_tests)

    def add_geometry_tests(self):
        """
        Adds the information tests.
        """
        tests = unittest.makeSuite(GT.GeometryTests)
        self._suite.addTests(tests)
        
    def add_manipulation_tests(self):
        """
        Adds the information tests.
        """
        tests = unittest.makeSuite(MT.ManipulationTests)
        self._suite.addTests(tests)

    def add_other_molecule_tests(self):
        """
        Adds the information tests.
        """
        tests = unittest.makeSuite(OMT.OtherMoleculeTests)
        self._suite.addTests(tests)

    def add_selection_tests(self):
        """
        Adds the information tests.
        """
        tests = unittest.makeSuite(ST.SelectionsTests)
        self._suite.addTests(tests)
        
        