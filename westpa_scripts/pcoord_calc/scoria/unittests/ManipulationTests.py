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
import copy

#import numpy as np
from scoria import dumbpy as np

import scoria
from ..six.moves import range


class ManipulationTests(unittest.TestCase):
    """
    Base Test Suite
    """

    # Initialization and destruction for each test.

    def setUp(self):
        """
        Setting up the test molecule.
        """
        info_path = os.path.dirname(os.path.abspath(__file__)) + '/../sample-files/'
        self.mol = scoria.Molecule(info_path + '3_mol_test.pdb')

        self.accuracy = 4

    def tearDown(self):
        """
        Cleans up variables for the next test.
        """
        self.mol = None

        # Remove all files from the tmp folder?

    ### Tests

    def test_set_coordinate_undo_point(self):
        """
        Empty test.
        """
        original = self.mol.get_coordinates_undo_point()
        self.assertIsNone(original)

        expected = self.mol.get_trajectory_coordinates()
        self.mol.set_coordinates_undo_point(expected)
        undo_point = self.mol.get_coordinates_undo_point()

        for h in range(self.mol.get_trajectory_frame_count()):
            for i in range(self.mol.get_total_number_of_atoms()):
                for j in range(3):
                    self.assertAlmostEqual(expected[h][i][j], undo_point[h][i][j], self.accuracy)

    def test_coordinate_undo(self):
        """
        Empty test.
        """
        shifted = [[None]]
        undone = [[None]]
        expected = [[None]]

        expected = self.mol.get_trajectory_coordinates()
        self.mol.set_coordinates_undo_point(expected)

        self.mol.set_atom_location(0, np.array([100.0, 100.0, 100.0]))

        # print(shifted[0][0], undone[0][0], expected[0][0])

        shifted = self.mol.get_trajectory_coordinates()

        # print(shifted[0][0], undone[0][0], expected[0][0])


        self.mol.coordinate_undo()
        undone = self.mol.get_trajectory_coordinates()

        # print(shifted[0][0], undone[0][0], expected[0][0])


        for h in range(self.mol.get_trajectory_frame_count()):
            for i in range(self.mol.get_total_number_of_atoms()):
                for j in range(3):
                    self.assertAlmostEqual(expected[h][i][j],
                        undone[h][i][j], self.accuracy, msg= str(i) + ' ' + str(j))
                    self.assertNotAlmostEqual(expected[h][i][j],
                        shifted[h][i][j], self.accuracy)

    def test_set_atom_location(self):
        """
        Empty test.
        """
        original = self.mol.get_coordinates()[1]
        delta = self.mol.set_atom_location(0, np.array([[20.0, 20.0, 20.0]]))

        self.assertEqual(len(delta[0]), 3)

        distance = (original + delta[0])
        moved = self.mol.get_coordinates()[1]

        self.assertEqual(list(distance), list(moved))

    def test_translate_molecule(self):
        """
        Empty test.
        """
        original = self.mol.get_coordinates()[0]
        delta = np.array([10.0, 10.0, 10.0])

        self.mol.translate_molecule(delta)

        moved = self.mol.get_coordinates()[0]
        distance = original + delta

        self.assertEqual(list(moved), list(distance))

    def test_rotate_molecule_around_a_line_between_points(self):
        """
        Empty test.
        """ 
        one = np.array([0.0, 0.0, 0.0])
        two = np.array([1.0, 0.0, 0.0])
        radians = np.radians(180.0)

        original = self.mol.get_coordinates()[0]

        self.mol.rotate_molecule_around_a_line_between_points(one, two, radians)

        new = self.mol.get_coordinates()[0]
        expected = [10.0, -10.0, -10.0]
        for i in range(0,3):
            self.assertAlmostEqual(list(new)[i], expected[i], self.accuracy)

    @unittest.skip("Needs test written")
    def test_rotate_molecule_around_a_line_between_atoms(self):
        """
        Empty test.
        """
        pass

    @unittest.skip("Needs test written")
    def test_rotate_molecule_around_pivot_point(self):
        """
        Empty test.
        """
        self.mol.rotate_molecule_around_pivot_point([0,0,0], 0.0, 180.0, 180.0)

    @unittest.skip("Needs test written")
    def test_rotate_molecule_around_pivot_atom(self):
        """
        Empty test.
        """
        pass
