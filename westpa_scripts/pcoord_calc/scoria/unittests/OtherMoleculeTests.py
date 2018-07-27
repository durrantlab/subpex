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
from ..six.moves import range


class OtherMoleculeTests(unittest.TestCase):
    """
    Base Test Suite for testing OtherMolecule functions
    """

    # Initialization and destruction for each test.

    def setUp(self):
        """
        Setting up the test molecules.
        """
        info_path = os.path.dirname(os.path.abspath(__file__)) + '/../sample-files/'
        self.mol = scoria.Molecule(info_path + '3_mol_test.pdb')
        self.other_mol = scoria.Molecule(info_path + 'other_mol_test.pdb')

        self.accuracy = 4

        self.tethers = []
        for i in range(12):
            self.tethers.append([i, i]) 

    def tearDown(self):
        """
        Cleans up variables for the next test.
        """
        self.mol = None

        # Remove all files from the tmp folder?

    ### Tests

    def test_get_other_molecules_aligned_to_this(self):
        """
        Empty test.
        """


        rotated = self.mol.get_other_molecules_aligned_to_this(self.other_mol, self.tethers)

        origianl_coords = self.mol.get_coordinates()
        aligned_coords = rotated.get_coordinates()

        for i in range(0, len(origianl_coords)):
            for j in range(0, 3):
                self.assertAlmostEqual(origianl_coords[i][j], aligned_coords[i][j], self.accuracy)


    def test_steric_clash_with_another_molecules(self):
        """
        Empty test.
        """
        should_clash = self.mol.steric_clash_with_another_molecules(self.mol, 5)
        self.assertTrue(should_clash)

        shouldnt_clash = self.mol.steric_clash_with_another_molecules(self.other_mol, 5)
        self.assertFalse(shouldnt_clash)


    def test_merge_with_another_molecules(self):
        """
        Empty test.
        """
        expected_total = self.mol.get_total_number_of_atoms() +\
                         self.other_mol.get_total_number_of_atoms()

        new_mol = self.mol.merge_with_another_molecules(self.other_mol)

        self.assertEqual(new_mol.get_total_number_of_atoms(), expected_total)


    def test_get_distance_to_another_molecule(self):
        """
        Empty test.
        """
        expected_distance = 20.0
        distance = self.mol.get_distance_to_another_molecule(self.other_mol)

        self.assertAlmostEqual(expected_distance, distance, self.accuracy)

        distance = self.mol.get_distance_to_another_molecule(self.mol)

        self.assertAlmostEqual(0.0, distance, self.accuracy)


    def test_get_rmsd_equivalent_atoms_specified(self):
        """
        Empty test.
        """
        rmsd = self.mol.get_rmsd_equivalent_atoms_specified(self.mol, self.tethers)

        self.assertAlmostEqual(rmsd, 0.0, self.accuracy)


    def test_get_rmsd_order_dependent(self):
        """
        Empty test.
        """
        rmsd = self.mol.get_rmsd_order_dependent(self.mol)

        self.assertAlmostEqual(rmsd, 0.0, self.accuracy)

        expected_rmsd = 86.019582323832893
        rmsd = self.mol.get_rmsd_order_dependent(self.other_mol)

        self.assertAlmostEqual(rmsd, expected_rmsd, self.accuracy)

    def test_get_rmsd_heuristic(self):
        """
        Empty test.
        """
        rmsd = self.mol.get_rmsd_heuristic(self.mol)

        self.assertAlmostEqual(rmsd, 0.0, self.accuracy)

        expected_rmsd = 86.019582323832893
        rmsd = self.mol.get_rmsd_order_dependent(self.other_mol)

        self.assertAlmostEqual(rmsd, expected_rmsd, self.accuracy)
        # We need to find out what the reduced RMSD is from, and if that's accurate
        