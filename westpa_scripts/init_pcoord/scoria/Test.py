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
from __future__ import print_function
from scoria.Molecule import Molecule
from .six.moves import range

try: import cStringIO as StringIO  # python2
except: from io import StringIO  # python3

import os
from scoria import dumbpy as numpy
import inspect
from time import time
import shutil
import math


class Test:
    """A class for testing all scoria functions."""

    mol = None
    numpy.missing_dependency_throws_error = False

    def test_all(self):
        """Test all scoria functions."""

        print("Creating directory to store temporary files.")

        if not os.path.exists("./scoria_tests_tmp"):
            os.mkdir("./scoria_tests_tmp")

        self.test_file_io()
        self.test_information()
        #self.test_selection()
        self.test_manipulation()
        #self.test_other_moleculess()
        self.test_atoms_and_bonds()
        self.test_geometry()

    def test_file_io(self):
        """Test the functions in FileIO."""

        file_io_filename = "./scoria_tests_tmp/file_io_test"
        sample_structures_dir = os.path.dirname(inspect.stack()[0][1]) + os.sep + "sample-files" + os.sep

        print("FileIO Functions")
        print("    load_pdb_into()")
        self.mol = Molecule()
        self.mol.load_pdb_into(sample_structures_dir + "single_frame.pdb", True, True, True)

        pdb_str = self.mol.save_pdb("", False, False, True)
        print("        Initial part of PDB string:")
        print((pdb_str[:500]))

        print("    load_pdb_into_using_file_object()")
        self.mol = Molecule()
        self.mol.load_pdb_into_using_file_object(
            StringIO.StringIO(pdb_str),
            True,
            True,
            True
        )

        print("    save_pdb()")
        self.mol.save_pdb(file_io_filename + ".pdb", True, True, False)

        print("    save_pym()")
        self.mol.save_pym(file_io_filename + ".pym", True, True, True, True, True)

        print("    load_pdbqt_into()")
        self.mol = Molecule()
        self.mol.load_pdbqt_into(sample_structures_dir + "single_frame.pdbqt", False, True, True)

        pdbqt_str = self.mol.save_pdb("", False, False, True)
        print("        Initial part of PDBQT string:")
        print((pdbqt_str[:500]))

        print("    load_pdbqt_into_using_file_object()")
        self.mol = Molecule()
        self.mol.load_pdb_into_using_file_object(
            StringIO.StringIO(pdbqt_str),
            False,  # no bonds by distance, because pdbqt file has unrecognized atom types.
            True,
            True
        )

        print("    load_MDAnalysis_into()")
        self.mol = Molecule()
        self.mol.load_MDAnalysis_into(sample_structures_dir + "M2_traj.psf",
                                     sample_structures_dir + "M2_traj.dcd")
    
        # Temporarily commented out because no dumbpy implementation needed.
        print("    load_pym_into()")
        self.mol = Molecule()
        self.mol.load_pym_into(file_io_filename + ".pym")

    def test_information(self):
        """Test the functions in Information"""

        print("Information Functions")
        print("    get_filename()")
        print(("        Filename: " + self.mol.get_filename()[0]))

        print("    get_remarks()")
        print(("        " + str(self.mol.get_remarks())))

        print("    assign_elements_from_atom_names()")
        self.mol.assign_elements_from_atom_names()

        print("    assign_masses()")
        self.mol.assign_masses()

        print("    get_center_of_mass()")
        print(("        " + str(self.mol.get_center_of_mass())))

        print("    get_atom_information()")
        print(("        " + str(self.mol.get_atom_information()[0])))

        print("    get_coordinates()")
        print(("        " + str(self.mol.get_coordinates()[0])))

        print("    get_bonds()")
        if numpy.class_dependency("test get_bonds()", "NUMPY"):
            if numpy.class_dependency("test set_bonds()", "SCIPY"):
                print(("        " + str(self.mol.get_bonds()[0][:10])))

        print("    get_geometric_center()")
        print(("        " + str(self.mol.get_geometric_center())))

        print("    get_total_number_of_atoms()")
        print(("        " + str(self.mol.get_total_number_of_atoms())))

        print("    get_total_number_of_heavy_atoms()")
        print(("        " + str(self.mol.get_total_number_of_heavy_atoms())))

        print("    get_total_mass()")
        print(("        " + str(self.mol.get_total_mass())))

        print("    get_bounding_box()")
        print(("        " + str(self.mol.get_bounding_box())))

        # Temporarily commented out because no dumbpy implementation needed.
        print("    get_bounding_sphere()")
        print(("        " + str(self.mol.get_bounding_sphere())))

        print("    get_constants()")
        print(("        " + str(self.mol.get_constants())))

        print("    belongs_to_protein()")
        print(("        " + str(self.mol.belongs_to_protein(1))))

        print("    belongs_to_dna()")
        print(("        " + str(self.mol.belongs_to_dna(1))))

        print("    belongs_to_rna()")
        print(("        " + str(self.mol.belongs_to_rna(1))))

        print("    set_filename()")
        self.mol.set_filename("test_name.pdb")
        print(("        " + str(self.mol.get_filename())))

        print("    set_remarks()")
        self.mol.set_remarks("Test remark.")
        print(("        " + str(self.mol.get_remarks())))

        print("    set_coordinates()")
        coors = self.mol.get_coordinates()
        
        for i in range(len(coors)):  # Doing it this way instead of coors[:,0] = 999.999
            coors[i][0] = 999.999    # so dumbpy will pass
        
        self.mol.set_coordinates(coors)
        print(("        " + str(self.mol.get_coordinates()[0])))

        print("    set_atom_information()")
        inf = self.mol.get_atom_information()
        inf["resname_padded"] = " TST "
        self.mol.set_atom_information(inf)
        print(("        " + str(self.mol.get_atom_information()[0])))

        print("    set_bonds()")
        if numpy.class_dependency("test set_bonds()", "NUMPY"):
            if numpy.class_dependency("test set_bonds()", "SCIPY"):
                bnds = self.mol.get_bonds()
                bnds[:,:10] = 1
                self.mol.set_bonds(bnds)
                print(("        " + str(self.mol.get_bonds()[0][:10])))

        print("    serial_reindex()")
        self.mol.serial_reindex()
        
        print("    resseq_reindex()")
        self.mol.resseq_reindex()
        
        print("    set_coordinates_undo_point()")
        self.mol.set_coordinates_undo_point(
            self.mol.get_coordinates() + 10
        )
        
        print("    define_molecule_chain_residue_spherical_boundaries()")
        self.mol.define_molecule_chain_residue_spherical_boundaries()
        
        print("    get_hierarchy()")
        if numpy.class_dependency("test the get_hierarchy() function", "NUMPY"):
            if numpy.class_dependency("test the get_hierarchy() function", "SCIPY"):
                print(("        " + str(self.mol.get_hierarchy()['residues']['indices']['VAL-32-A'])))

        print("    set_hierarchy()")
        if numpy.class_dependency("test the set_hierarchy() function", "NUMPY"):
            self.mol.set_hierarchy(
                self.mol.get_hierarchy()
            )

    def test_selection(self):
        """Test the functions in Selections."""

        print("Selection Functions")

        # Get new molecule... clean slate.
        self.mol = Molecule()
        self.mol.load_pdb_into_using_file_object(
            StringIO.StringIO(sample_structures_dir + "single_frame.pdb"),
            True,
            True,
            True
        )

        print("    select_all()")
        print(("        Atoms in selection: " + str(
            len(
                self.mol.select_all()
            )
        )))

        print("    select_atoms()")

        sel = self.mol.select_atoms({"resname": "TRP"})
        print("        Atoms in selection: " + str(len(sel)))

        print("    invert_selection()")
        print(("        Atoms in selection: " + str(len(self.mol.invert_selection(sel)))))
        
        print("    select_all_atoms_bound_to_selection()")
        if numpy.class_dependency("test select_all_atoms_bound_to_selection()", "NUMPY"):
            if numpy.class_dependency("test select_all_atoms_bound_to_selection()", "SCIPY"):
                print(("        Atoms in selection: " + str(len(self.mol.select_all_atoms_bound_to_selection(sel)))))
        
        # Temporarily commented out because no dumbpy implementation needed.
        print("    select_atoms_near_other_selection()")
        if numpy.class_dependency("test select_atoms_near_other_selection()", "NUMPY"):
            if numpy.class_dependency("test select_atoms_near_other_selection()", "SCIPY"):
                print(("        Atoms in selection: " + str(len(self.mol.select_atoms_near_other_selection(sel, 8.0)))))
        
        print("    select_atoms_in_same_residue()")
        print(("        Atoms in selection: " + str(len(self.mol.select_atoms_in_same_residue([1])))))
        
        print("    select_atoms_from_same_molecule()")
        if numpy.class_dependency("test select_atoms_from_same_molecule()", "NUMPY"):        
            if numpy.class_dependency("test select_atoms_from_same_molecule()", "SCIPY"):        
                print(("        Atoms in selection: " + str(len(self.mol.select_atoms_from_same_molecule([1])))))
        
        print("    select_atoms_in_bounding_box()")
        if numpy.class_dependency("test select_atoms_in_bounding_box()", "NUMPY"):        
            print(("        Atoms in selection: " + str(
                len(
                    self.mol.select_atoms_in_bounding_box(
                        self.mol.get_bounding_box(list(range(15)), 5.0)
                    )
                )
            )))

        print("    select_branch()")  # Get a side chain.
        if numpy.class_dependency("test select_branch()", "NUMPY"):        
            if numpy.class_dependency("test select_branch()", "SCIPY"):        
                print(("        Atoms in selection: " + str(len(self.mol.select_branch(1, 4)))))

        print("    selections_of_constituent_molecules()" )
        if numpy.class_dependency("test select_branch()", "NUMPY"):
            if numpy.class_dependency("test select_branch()", "SCIPY"):
                print(("        Atoms in selection: " + str(self.mol.selections_of_constituent_molecules())))

        print("    selections_of_chains()")
        if numpy.class_dependency("test selections_of_chains()", "NUMPY"):
            print(("        Chains mapped to indices: " + str(list(self.mol.selections_of_chains().keys()))))

        print("    selections_of_residues()")
        if numpy.class_dependency("test selections_of_residues()", "NUMPY"):
            print(("        Residues mapped to indices: " + str(list(self.mol.selections_of_residues().keys())[:5])))

        # Temporarily commented out because no dumbpy implementation needed.
        print("    select_close_atoms_from_different_molecules()")
        if numpy.class_dependency("test select_close_atoms_from_different_molecules()", "NUMPY"):
            if numpy.class_dependency("test select_close_atoms_from_different_molecules()", "SCIPY"):
                sels = self.mol.select_close_atoms_from_different_molecules(self.mol, 1.0)
                print(("        Atoms in mol1 selection: " + str(len(sels[0]))))
                print(("        Atoms in mol2 selection: " + str(len(sels[1]))))

        print("    get_molecule_from_selection()")
        sel = [1,2,3,4,5]
        print(("        Atoms in selection: " + str(len(sel))))
        mol2 = self.mol.get_molecule_from_selection(sel)
        mol2.save_pdb("./scoria_tests_tmp/save_selection.pdb", True, True, False)

    def test_manipulation(self):
        """Test the functions in Manipulation."""

        manip_filename = "./scoria_tests_tmp/manipulation_test"

        print("Manipulation Functions")
        
        print("    set_coordinates_undo_point()")
        self.mol.set_coordinates_undo_point(
            self.mol.get_coordinates()
        )
        
        print("    translate_molecule()")
        self.mol.translate_molecule(numpy.array([10.0, 10.0, 10.0]))
        self.mol.save_pdb(manip_filename + "1.pdb", False, False, False)
        
        print("    set_atom_location()")
        self.mol.set_atom_location(1, numpy.array([10.0, 10.0, 10.0]))
        self.mol.save_pdb(manip_filename + "2.pdb", False, False, False)
        
        print("    rotate_molecule_around_a_line_between_points()")
        self.mol.rotate_molecule_around_a_line_between_points(
            numpy.array([0.0, 0.0, 0.0]),
            numpy.array([10.0, 0.0, 0.0]),
            45.0
        )
        self.mol.save_pdb(manip_filename + "3.pdb", False, False, False)
        
        print("    rotate_molecule_around_a_line_between_atoms()")
        self.mol.rotate_molecule_around_a_line_between_atoms(1, 2, 45.0)
        self.mol.save_pdb(manip_filename + "4.pdb", False, False, False)
        
        print("    rotate_molecule_around_pivot_point()")
        if numpy.class_dependency("test rotate_molecule_around_pivot_point(). Missing the dot-product function", "NUMPY"):
            self.mol.rotate_molecule_around_pivot_point(numpy.array([0.0, 0.0, 0.0]), 45, 45, 45)
            self.mol.save_pdb(manip_filename + "5.pdb", False, False, False)
        
        print("    rotate_molecule_around_pivot_atom()")
        if numpy.class_dependency("test rotate_molecule_around_pivot_atom(). Missing the dot-product function", "NUMPY"):
            self.mol.rotate_molecule_around_pivot_atom(5, 45, 45, 45)
            self.mol.save_pdb(manip_filename + "6.pdb", False, False, False)
        
        print("    coordinate_undo()")
        self.mol.coordinate_undo()
        self.mol.save_pdb(manip_filename + "7.pdb", False, False, False)
        
    def test_other_moleculess(self):
        """Test the functions in OtherMolecules"""

        print("OtherMolecules Functions")

        # Make two molecules to play with.
        mol1 = Molecule()
        mol1.load_pdb_into_using_file_object(
            StringIO.StringIO(sample_structures_dir + "single_frame.pdb"),
            True,
            True,
            True
        )

        mol2 = Molecule()
        mol2.load_pdb_into_using_file_object(
            StringIO.StringIO(sample_structures_dir + "single_frame.pdb"),
            True,
            True,
            True
        )

        mol2.translate_molecule(numpy.array([10.0, 10.0, 10.0]))

        # Temporarily commented out because no dumbpy implementation needed.
        print("    steric_clash_with_another_molecule()")
        print(("        " + str(mol1.steric_clash_with_another_molecules(mol2, 5.0, False))))
        print(("        " + str(mol1.steric_clash_with_another_molecules(mol2, 5.0, True))))

        # Temporarily commented out because no dumbpy implementation needed.
        print("    get_distance_to_another_molecule()")
        print(("        " + str(mol1.get_distance_to_another_molecule(mol2, False))))
        print(("        " + str(mol1.get_distance_to_another_molecule(mol2, True))))

        print("    get_rmsd_order_dependent()")
        print(("        " + str(mol1.get_rmsd_order_dependent(mol2))))

        # Temporarily commented out because no dumbpy implementation needed.
        print("    get_rmsd_heuristic()")
        print(("        " + str(mol1.get_rmsd_heuristic(mol2))))

        print("    get_rmsd_equivalent_atoms_specified()")
        tethers = [numpy.array([1, 1]), numpy.array([2, 2]), numpy.array([3, 3])]
        print(("        " + str(mol1.get_rmsd_equivalent_atoms_specified(mol2, tethers))))

        print("    merge_with_another_molecules()")
        merged_mol = mol1.merge_with_another_molecules(mol2)
        merged_mol.save_pdb("./scoria_tests_tmp/merged.pdb", False, False, False)

        print("    get_other_molecules_aligned_to_this()")
        if numpy.class_dependency("test get_other_molecules_aligned_to_this(). Missing the dot-product function", "NUMPY"):
            aligned_mol = mol1.get_other_molecules_aligned_to_this(mol2, tethers)
            print(("        New RMSD: " + str(mol1.get_rmsd_order_dependent(aligned_mol))))

    def test_atoms_and_bonds(self):
        """Test the functions in AtomsAndBonds."""

        self.mol = Molecule()
        
        print("AtomsAndBonds Functions")

        print("    add_atom()")
        self.mol.add_atom()
        self.mol.add_atom()
        self.mol.add_atom()
        self.mol.add_atom()
        self.mol.add_atom()

        self.mol.set_coordinates(numpy.array(
            [[0.0, 0.0, 0.0],
             [1.5, 0.0, 0.0],
             [1.5, 1.5, 0.0],
             [3.0, 1.5, 0.0],
             [999.99, 2.0, 0.0]]
        ))

        print("    delete_atom()")
        self.mol.delete_atom(4)

        print("    add_bond()")
        self.mol.add_bond(0, 1, 2)

        print("    delete_bond()")
        self.mol.delete_bond(0, 1)

        # Temporarily commented out because no dumbpy implementation needed.
        print("    create_bonds_by_distance()")
        self.mol.create_bonds_by_distance()

        print("    get_number_of_bond_partners_of_element()")
        print(("        " + str(self.mol.get_number_of_bond_partners_of_element(0, "X"))))

        print("    get_index_of_first_bond_partner_of_element()")
        print(("        " + str(self.mol.get_index_of_first_bond_partner_of_element(1, "X"))))

    def test_geometry(self):
        """Test the functions in Geometry."""

        print("Geometry Functions")

        coors = self.mol.get_coordinates()

        print("    get_angle_between_three_points()")
        print(("        " + str(self.mol.get_angle_between_three_points(coors[0], coors[1], coors[2]))))

        print("    get_dihedral_angle()")
        print(("        " + str(self.mol.get_dihedral_angle(coors[0], coors[1], coors[2], coors[3]))))

        print("    is_planar()")
        print(("        " + str(self.mol.is_planar(coors[0], coors[1], coors[2], coors[3]))))

        print("    get_planarity_deviation()")
        print(("        " + str(self.mol.get_planarity_deviation(coors[0], coors[1], coors[2], coors[3]))))

class FileIOBenchmarks:
    """A class for testing the load and save times of multiple
    scoria-supported file formats."""

    molecule = None
    times = []
    sample_structures_dir = os.path.dirname(inspect.stack()[0][1]) + os.sep + "sample-files" + os.sep
    load_filename = ""
    save_filename = ""

    def __init__(self):
        self.load_filename = self.sample_structures_dir + "single_frame.pdb"
        self.save_filename = "tmptmp.pdb"

        # First, test PDB loading.
        print("Load PDB, no Bonds by Distance")
        self.timeit(self.load_pdb_no_bonds_by_distance, self.reset_test_vars_standard, self.new_molecule)
        print()

        print("Load PDB, with Bonds by Distance")
        self.timeit(self.load_pdb_with_bonds_by_distance, self.reset_test_vars_standard, self.new_molecule)
        print()

        print("Save PDB")
        self.timeit(self.save_pdb, self.reset_test_vars_before_saving, self.nothing)
        print()

        print("Save PYM")
        self.save_filename = "tmptmp.pym"
        self.timeit(self.save_pym, self.reset_test_vars_before_saving, self.nothing)
        print()

        print("Load PYM")
        self.load_filename = "tmptmp.pym"
        self.timeit(self.load_pym, self.reset_test_vars_standard, self.new_molecule)
        print()

        print("Load DCD (single frame)")
        self.timeit(self.load_dcd, self.reset_test_vars_standard, self.new_molecule)
        print()

        print("Load DCD (100 frames)")
        self.timeit(self.load_dcd_100_frames, self.reset_test_vars_standard, self.new_molecule)
        print()

        # Clean up
        if os.path.exists("tmptmp.pdb"):
            os.unlink("tmptmp.pdb")
        if os.path.exists("tmptmp.pym"):
            shutil.rmtree("tmptmp.pym")
    
    def mean(self, nums):
        asum = 0.0
        for num in nums:
            asum = asum + num
        return asum / len(nums)
    
    def std(self, nums):
        mn = self.mean(nums)
        sum_of_square_diffs = 0.0
        for num in nums:
            sum_of_square_diffs = sum_of_square_diffs + (num - mn)**2
        variance = sum_of_square_diffs / len(nums)
        std = math.sqrt(variance)
        return std

    def timeit(self, func, reset_test_vars, before_each_timing):
        reset_test_vars()
        for i in range(100):
            before_each_timing()
            t1 = time()
            func()
            t2 = time()
            self.times.append(t2 - t1)
        print((self.mean(self.times), "+/-", self.std(self.times)))

    def reset_test_vars_standard(self):
        self.times = []

    def reset_test_vars_before_saving(self):
        self.times = []
        self.molecule = Molecule()
        self.molecule.load_pdb_into(self.load_filename, False, False, False)
        #if os.path.exists(self.load_filename):
        #    os.unlink(self.load_filename)
        if os.path.exists(self.save_filename):
            try:
                shutil.rmtree(self.save_filename)
            except:
                os.unlink(self.save_filename)
    
    def new_molecule(self):
        self.molecule = Molecule()

    def nothing(self): pass

    def load_pdb_no_bonds_by_distance(self):
        self.molecule.load_pdb_into(self.load_filename, False, False, False)
        
    def load_pdb_with_bonds_by_distance(self):
        self.molecule.load_pdb_into(self.load_filename, True, False, False)
    
    def save_pdb(self):
        self.molecule.save_pdb("tmptmp.pdb", False, False, False)
    
    def save_pym(self):
        self.molecule.save_pym("tmptmp.pym", False, False, False, False, False)

    def load_pym(self):
        self.molecule.load_pym_into(self.load_filename)
    
    def load_dcd(self):
        if "MDANALYSIS" in numpy.dependencies_available:
            self.molecule.load_MDAnalysis_into(self.sample_structures_dir + "single_frame.psf", self.sample_structures_dir + "single_frame.dcd")
    
    def load_dcd_100_frames(self):
        if "MDANALYSIS" in numpy.dependencies_available:
            self.molecule.load_MDAnalysis_into(self.sample_structures_dir + "single_frame.psf", self.sample_structures_dir + "single_frame.100.dcd")
        
