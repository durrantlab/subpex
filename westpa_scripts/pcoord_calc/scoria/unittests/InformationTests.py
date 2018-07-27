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

from scoria.six.moves import range


class InformationTests(unittest.TestCase):
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
    # Testing Getters

    def test_get_filename(self):
        """
        Tests the getting of filenames.
        """
        expected_filename = [os.path.dirname(os.path.abspath(__file__)) + \
         '/../sample-files/3_mol_test.pdb']

        self.assertEqual(self.mol.get_filename(), expected_filename)

    def test_get_remarks(self):
        """
        Tests the getting of remarks.
        """
        expected_remarks = [" This is a test file."]
        self.assertEqual(self.mol.get_remarks(), expected_remarks)


    def test_get_center_of_mass(self):
        """
        Tests the determination of the center of mass.
        """

        center_of_mass = self.mol.get_center_of_mass()


        center = center_of_mass  # pypy case. Basically now just error
                                        # checking for pypy, since below will
                                        # always be true.

        self.assertAlmostEqual(center_of_mass[0], center[0], self.accuracy)
        self.assertAlmostEqual(center_of_mass[1], center[1], self.accuracy)
        self.assertAlmostEqual(center_of_mass[2], center[2], self.accuracy)


    def test_get_atom_information(self):
        """
        Tests the atom information.
        """
        atom_inf = self.mol.get_atom_information()

        expected_record_name = ['ATOM  '] * self.mol.get_total_number_of_atoms()
        self.assertListEqual(list(atom_inf['record_name']), expected_record_name)

        expected_serial = list(range(1, 13))
        self.assertListEqual(list(atom_inf['serial']), expected_serial)

        expected_names = ['N1', "C2'", 'N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1',
                          'CD2', 'CE1', 'NE2']
        self.assertListEqual(list(atom_inf['name']), expected_names)

        expected_element = ['N', 'C', 'N', 'C', 'C', 'O', 'C', 'C', 'N', 'C',
                            'C', 'N']
        self.assertListEqual(list(atom_inf['element']), expected_element)

        expected_resname = ['U', 'DT', 'HIS', 'HIS', 'HIS', 'HIS', 'HIS',
                            'HIS', 'HIS', 'HIS', 'HIS', 'HIS']
        self.assertListEqual(list(atom_inf['resname']), expected_resname)

    def test_padded_atom_informatio_fields_are_correct(self):
        """
        Testing that the atom information fields are properly stripped.
        """
        atom_inf = self.mol.get_atom_information()

        padding = ['name', 'chainid', 'resname', 'element']

        for field in padding:
            for i in range(self.mol.get_total_number_of_atoms()):
                self.assertEqual(atom_inf[field][i], atom_inf[field+'_padded'][i].strip())

    def test_get_coordinates(self):
        """
        Tests that the coordinates are returned as expected.
        """
        coordinates = self.mol.get_coordinates()
        self.assertEqual(len(coordinates), 12)

        atom_zero = coordinates[0]
        expected_atom = [10.000, 10.000, 10.000]
        self.assertAlmostEqual(atom_zero[0], expected_atom[0], self.accuracy)
        self.assertAlmostEqual(atom_zero[1], expected_atom[1], self.accuracy)
        self.assertAlmostEqual(atom_zero[2], expected_atom[2], self.accuracy)

    def test_get_bonds(self):
        """
        Tests that the bond values are returned.
        """
        bonds = self.mol.get_bonds()

        # Checking that it is the correct depth, length and value
        self.assertEqual(len(bonds), 12)
        self.assertEqual(len(bonds[0]), 12)
        self.assertEqual(bonds[0][0], 0)
        self.assertEqual(bonds[2][3], 1)

    def test_get_geometric_center(self):
        """
        Tests that the geometric center is able to be calculated properly.
        """
        geo_center = self.mol.get_geometric_center()

        center = geo_center  # pypy case. Basically now just error
                                    # checking for pypy, since below will
                                    # always be true.

        self.assertAlmostEqual(geo_center[0], center[0], self.accuracy)
        self.assertAlmostEqual(geo_center[1], center[1], self.accuracy)
        self.assertAlmostEqual(geo_center[2], center[2], self.accuracy)

    def test_get_total_number_of_atoms(self):
        """
        Tests that the number of atoms returned is correct.
        """
        number_of_atoms = self.mol.get_total_number_of_atoms()
        self.assertEqual(number_of_atoms, 12)

    def test_get_total_mass(self):
        """
        Tests that the total mass is returned correctly.
        """
        total_mass = self.mol.get_total_mass()

        expected_mass = total_mass  # pypy case. Basically now just error
                                    # checking for pypy, since below will
                                    # always be true.

        self.assertAlmostEqual(total_mass, expected_mass, 1)

    # Depreciated? And needs skip for dependencies
    @unittest.skip("hierarchy related method")
    def test_get_heirarchy(self):
        """
        Tests that the hierarchy can be set.
        """
        hierarchy = self.mol.get_hierarchy()
        # Assertation here

    ## Testing Setters

    def test_set_filename(self):
        """
        Tests the setting of filenames.
        """
        expected_file = ["OtherFile.txt"]
        self.mol.set_filename(expected_file)
        self.assertEqual(self.mol.get_filename(), expected_file)

    def test_set_remarks(self):
        """
        Tests the setting of remarks.
        """
        set_remarks = ["TEST REMARK"]
        self.mol.set_remarks(set_remarks)
        self.assertEqual(self.mol.get_remarks(), set_remarks)


    def test_set_coordinates(self):
        """
        Tests that the coordinates can be set
        """
        atom_count = self.mol.get_total_number_of_atoms()
        coordinates = [[0.0, 0.0, 0.0]] * atom_count
        self.mol.set_coordinates(coordinates)
        self.assertEqual(self.mol.get_coordinates(), coordinates)


    def test_set_atomic_information(self):
        """
        Tests that the atom information can be set.
        """
        atom_inf = self.mol.get_atom_information()

        atom_inf['chainid'][0] = 'X'

        self.mol.set_atom_information(atom_inf)

        self.assertEqual(self.mol.get_atom_information()['chainid'][0], 'X')

    def test_set_bonds(self):
        """
        Tests that the atom information can be set.
        """
        atom_count = self.mol.get_total_number_of_atoms()
        bonds = [[0] * atom_count] * atom_count
        self.mol.set_bonds(bonds)

        self.assertEqual(self.mol.get_bonds(), bonds)

    def test_set_coordinate_undo_point(self):
        """
        Tests that the coordinate undo point can be set.
        """
        expected = {}
        self.mol.set_coordinates_undo_point(expected)
        coord_undo = self.mol.get_coordinates_undo_point()

        self.assertEqual(expected, coord_undo)

    ## Testing Functions

    # The bounding box, having several parameters, should have some
    # comprehensive tests written for it.
    @unittest.skip("Need Bounding Box")
    def test_get_default_bounding_box(self):
        """
        Tests that the bounding box can be calculated.
        """
        test_box = []
        bounding_box = self.mol.get_bounding_box(None, 0.0, 0)
    
        self.assertAlmostEqual(bounding_box[0][0], test_box[0][0], self.accuracy)
        self.assertAlmostEqual(bounding_box[0][1], test_box[0][1], self.accuracy)
        self.assertAlmostEqual(bounding_box[0][2], test_box[0][2], self.accuracy)
    
        self.assertAlmostEqual(bounding_box[1][0], test_box[1][0], self.accuracy)
        self.assertAlmostEqual(bounding_box[1][1], test_box[1][1], self.accuracy)
        self.assertAlmostEqual(bounding_box[1][2], test_box[1][2], self.accuracy)

    # Similar to the bounding box tests, we need to check that all
    # parameters work properly.
    @unittest.skip("Need Bounding SPhere")
    def test_get_bounding_sphere(self):
        """
        Tests that the bounding sphere can be calculated.
        """
        sphere = []
        bounding_sphere = self.mol.get_bounding_sphere(None, 0.0, 0)
    
        self.assertAlmostEqual(bounding_sphere[0][0], sphere[1][0], self.accuracy)
        self.assertAlmostEqual(bounding_sphere[0][1], sphere[1][1], self.accuracy)
        self.assertAlmostEqual(bounding_sphere[0][2], sphere[1][2], self.accuracy)
    
        self.assertAlmostEqual(bounding_sphere[1], sphere[0], self.accuracy)

    @unittest.skip("Needs test written")
    def test_get_constants(self):
        """
        Tests that the constants returned are as expected. How do we want to test this?
        """
        constants = self.mol.get_constants()
        # Assertation here

    # For the 'belongs' tests, we need one index of each category
    # And each should be tested for redundancy
    def test_belongs_to_protein(self):
        """
        Tests that indices are proteins
        """
        self.assertFalse(self.mol.belongs_to_protein(0))
        self.assertFalse(self.mol.belongs_to_protein(1))
        self.assertTrue(self.mol.belongs_to_protein(2))

    def test_belongs_to_dna(self):
        """
        Tests that indices are DNA
        """
        self.assertFalse(self.mol.belongs_to_dna(0))
        self.assertTrue(self.mol.belongs_to_dna(1))
        self.assertFalse(self.mol.belongs_to_dna(2))

    def test_belongs_to_RNA(self):
        """
        Tests that indices are RNA
        """
        self.assertTrue(self.mol.belongs_to_rna(0))
        self.assertFalse(self.mol.belongs_to_rna(1))
        self.assertFalse(self.mol.belongs_to_rna(2))

    def test_assign_elements_from_atom_names(self):
        """
        Tests the assignment of elements from the atom names.
        """
        atom_inf = self.mol.get_atom_information()
        other = copy.deepcopy(self.mol.get_atom_information())

        atoms = self.mol.get_total_number_of_atoms()

        atom_inf['element'] = " "
        atom_inf['element_padded'] = " "
        self.mol.set_atom_information(atom_inf)

        for i in range(atoms):
            self.assertNotEqual(self.mol.get_atom_information()['element'][i],
                                other['element'][i])
            self.assertNotEqual(self.mol.get_atom_information()['element_padded'][i],
                                other['element_padded'][i])


        self.mol.assign_elements_from_atom_names()

        for i in range(atoms):
            self.assertEqual(self.mol.get_atom_information()['element'][i],
                             other['element'][i])
            self.assertEqual(self.mol.get_atom_information()['element_padded'][i],
                             other['element_padded'][i])


    def test_assign_masses(self):
        """
        Tests the assignment of masses.
        """
        atom_inf = self.mol.get_atom_information()
        masses = self.mol.get_constants()['mass_dict']

        atoms = self.mol.get_total_number_of_atoms()

        self.mol.set_atom_information(atom_inf)

        # At times mass entry isn't created until after self.mol.assign_masses()
        # with self.assertRaises(ValueError):
        #     self.mol.get_atom_information()['mass']

        self.mol.assign_masses()

        for i in range(atoms):
            element = self.mol.get_atom_information()['element'][i]
            self.assertEqual(self.mol.get_atom_information()['mass'][i],
                             masses[element])

    def test_serial_reindex(self):
        """
        Tests the reindexing of the serial field.
        """

        # pypy can't delete atoms. So let's try a different test here.
        #self.mol.delete_atom(4)
        
        atom_inf = self.mol.get_atom_information()
        atom_inf['serial'][5] = 999
        self.mol.set_atom_information(atom_inf)

        atoms = self.mol.get_total_number_of_atoms()
        other = list(range(1, atoms+1))

        self.assertNotEqual(list(self.mol.get_atom_information()['serial']), other)

        self.mol.serial_reindex()

        for i in range(atoms):
            self.assertEqual(self.mol.get_atom_information()['serial'][i],
                             other[i])

        # Assertion here, post assignment

    def test_resseq_reindex(self):
        """
        Tests the reindexing of the resseq field.
        """
        self.mol.delete_atom(4)
        atom_inf = self.mol.get_atom_information()
        atoms = self.mol.get_total_number_of_atoms()

        other = [1, 2] + [3] * 10

        self.assertNotEqual(list(self.mol.get_atom_information()['resseq']), other)

        self.mol.resseq_reindex()

        for i in range(atoms):
            self.assertEqual(self.mol.get_atom_information()['resseq'][i],
                             other[i])

    @unittest.skip("hierarchy related method")
    def test_define_molecule_chain_residue_spherical_boundaries(self):
        """
        Tests the reindexing of the serial field.
        """
        # Assertion here, pre assignment
        self.mol.define_molecule_chain_residue_spherical_boundaries()
        # Assertion here, post assignment
        