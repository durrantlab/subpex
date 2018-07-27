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
import scoria

class InformationTests(unittest.TestCase):
    """
    Base Test Suite
    """

    # Initialization and destruction for each test.

    def setUp(self):
        """
        Setting up the test molecule.
        """
        
        self.mol = scoria.Molecule("PDB", "TEST.pdb")

    def tearDown(self):
        """
        Cleans up variables for the next test.
        """
        self.mol = None

        # Remove all files from the tmp folder?

    ### Tests
    # Testing Getters

    @unittest.skip("Needs final value")
    def test_get_filename(self):
        """
        Tests the getting of filenames.
        """
        # Fill out filename.
        self.assertEqual(self.mol.get_filename(), "TEST.pdb")

    @unittest.skip("Needs final value")
    def test_get_remarks(self):
        """
        Tests the getting of remarks.
        """
        self.assertEqual(self.mol.get_remarks(), "REMARKS")

    @unittest.skip("Need final value")
    def test_get_center_of_mass(self):
        """
        Tests the determination of the center of mass.
        """
        center_of_mass = self.mol.get_center_of_mass()
        self.assertEqual(center_of_mass, ["Expected"])

    @unittest.skip("Needs test written")
    def test_get_atom_information(self):
        """
        Tests the atom information.
        """
        atom_inf = self.mol.get_atom_information()
        # Assertation here, post

    @unittest.skip("Needs test written")
    def test_get_coordinates(self):
        """
        Tests that the coordinates are returned as expected.
        """
        coordinates = self.mol.get_coordinates()
        # Assertation here, post

    @unittest.skip("Needs test written")
    def test_get_bonds(self):
        """
        Tests that the bond values are returned.
        """
        bonds = self.mol.get_bonds()
        # Assertation here

    @unittest.skip("Needs test written")
    def test_get_geometric_center(self):
        """
        Tests that the geometric center is able to be calculated properly.
        """
        geo_center = self.mol.get_geometric_center()
        # Assertation here

    @unittest.skip("Needs test written")
    def test_get_total_number_of_atoms(self):
        """
        Tests that the number of atoms returned is correct.
        """
        number_of_atoms = self.mol.get_total_number_of_atoms()
        # Assertation here

    @unittest.skip("Needs test written")
    def test_get_total_mass(self):
        """
        Tests that the total mass is returned correctly.
        """
        total_mass = self.mol.get_total_mass()
        # Assertation here

    # Depreciated? And needs skip for dependencies
    @unittest.skip("Needs test written")
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
        self.mol.set_filename("OtherFile.txt")
        self.assertEqual(self.mol.get_filename(), "OtherFile.txt")

    def test_set_remarks(self):
        """
        Tests the setting of remarks.
        """
        self.mol.set_remarks("TEST REMARK")
        self.assertEqual(self.mol.get_remarks(), "TEST REMARK")

    @unittest.skip("Needs test written")
    def test_set_coordinates(self):
        """
        Tests that the coordinates can be set
        """
        coordinates = []
        self.mol.set_coordinates(coordinates)
        # Assertation here

    @unittest.skip("Needs test written")
    def test_set_atomic_information(self):
        """
        Tests that the atom information can be set.
        """
        atom_inf = {}
        self.mol.set_atom_information(atom_inf)
        # Assertation here

    # Add skip for dependencies
    @unittest.skip("Needs test written")
    def test_set_bonds(self):
        """
        Tests that the atom information can be set.
        """
        bonds = {}
        self.mol.set_bonds(bonds)
        # Assertation here

    @unittest.skip("Needs test written")
    def test_set_coordinate_undo_point(self):
        """
        Tests that the coordinate undo point can be set.
        """
        coord_undo = {}
        self.mol.set_coordinates_undo_point(coord_undo)
        # Assertation here

    # Depreciated? And needs skip for dependencies
    @unittest.skip("Needs test written")
    def test_set_heirarchy(self):
        """
        Tests that the hierarchy can be set.
        """
        hierarchy = {}
        self.mol.set_hierarchy(hierarchy)
        # Assertation here


    ## Testing Functions

    # The bounding box, having several parameters, should have some
    # comprehensive tests written for it.
    @unittest.skip("Needs test written")
    def test_get_default_bounding_box(self):
        """
        Tests that the bounding box can be calculated.
        """
        bounding_box = self.mol.get_bounding_box(None, 0.0, 0)
        # Assertation here

    # Similar to the bounding box tests, we need to check that all
    # parameters work properly.
    @unittest.skip("Needs test written")
    def test_get_bounding_sphere(self):
        """
        Tests that the bounding sphere can be calculated.
        """
        bounding_sphere = self.mol.get_bounding_sphere(None, 0.0, 0)
        # Assertation here

    @unittest.skip("Needs test written")
    def test_get_constants(self):
        """
        Tests that the constants returned are as expected
        """
        constants = self.mol.get_constants(self)
        # Assertation here

    # For the 'belongs' tests, we need one index of each category
    # And each should be tested for redundancy
    @unittest.skip("Needs correct test indexes")
    def test_belongs_to_protein(self):
        """
        Tests that indices are proteins
        """
        self.assertTrue(self.mol.belongs_to_protein(0))
        self.assertFalse(self.mol.belongs_to_protein(1))
        self.assertFalse(self.mol.belongs_to_protein(2))

    @unittest.skip("Needs correct test indexes")
    def test_belongs_to_dna(self):
        """
        Tests that indices are DNA
        """
        self.assertTrue(self.mol.belongs_to_dna(0))
        self.assertFalse(self.mol.belongs_to_dna(1))
        self.assertFalse(self.mol.belongs_to_dna(2))

    @unittest.skip("Needs correct test indexes")
    def test_belongs_to_RNA(self):
        """
        Tests that indices are RNA
        """
        self.assertTrue(self.mol.belongs_to_rna(0))
        self.assertFalse(self.mol.belongs_to_rna(1))
        self.assertFalse(self.mol.belongs_to_rna(2))

    @unittest.skip("Needs test written")
    def test_assign_elements_from_atom_names(self):
        """
        Tests the assignment of elements from the atom names.
        """
        # Assertion here, pre assignment
        self.mol.assign_elements_from_atom_names()
        # Assertion here, post assignment

    @unittest.skip("Needs test written")
    def test_assign_masses(self):
        """
        Tests the assignment of masses.
        """
        # Assertion here, pre assignment
        self.mol.assign_masses()
        # Assertion here, post assignment

    @unittest.skip("Needs test written")
    def test_serial_reindex(self):
        """
        Tests the reindexing of the serial field.
        """
        # Assertion here, pre assignment
        self.mol.serial_reindex()
        # Assertion here, post assignment

    @unittest.skip("Needs test written")
    def test_resseq_reindex(self):
        """
        Tests the reindexing of the resseq field.
        """
        # Assertion here, pre assignment
        self.mol.resseq_reindex()
        # Assertion here, post assignment

    @unittest.skip("Needs test written")
    def test_define_molecule_chain_residue_spherical_boundaries(self):
        """
        Tests the reindexing of the serial field.
        """
        # Assertion here, pre assignment
        self.mol.define_molecule_chain_residue_spherical_boundaries()
        # Assertion here, post assignment
