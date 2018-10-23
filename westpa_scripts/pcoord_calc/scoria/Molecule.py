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
from scoria import dumbpy as numpy
from .FileIO import FileIO
from .AtomsAndBonds import AtomsAndBonds
from .Selections import Selections
from .Manipulation import Manipulation
from .Information import Information
from .OtherMolecules import OtherMolecules
from .Geometry import Geometry
import copy


class Molecule: # here's the actual Molecule class
    """
    Loads, saves, and manupulates molecuar models. The main scoria
    class. Contains Wrapper functions for subclasses.

    Examples assume::

        >>> import scoria
        >>> PDB = "./test_file.pdb"
        >>> mol = scoria.Molecule()
        >>> mol.load_pdb_into(PDB)
    """

    def __init__(self, *args):
        """
        Initializes the variables of the Molecule class.

        :param files args: Additional files added to the end of the paramter
            will be input with the method indicated by the fileType parameter.
        """

        self.fileio = FileIO(self)
        self.atoms_and_bonds = AtomsAndBonds(self)
        self.selections = Selections(self)
        self.manipulation = Manipulation(self)
        self.information = Information(self)
        self.other_molecules = OtherMolecules(self)
        self.geometry = Geometry(self)

        # Based on the file type, we will attempt to open the file.
        if len(args) > 0:
            if len(args) == 1:
                file = args[0]
                file_type = file.split('.')[-1].upper()

                if file_type == 'PDB':
                    self.load_pdb_trajectory_into(file)
                elif file_type == 'PDBQT':
                    self.load_pdbqt_trajectory_into(file)
                elif file_type == 'PYM':
                    self.load_pym_into(file)

    # Information methods
    ### Wrappers ###
    # Gets
    def get_coordinates(self, frame = None):
        """
        Returns the set of coordinates from the specified frame.

        Wrapper function for :meth:`~scoria.Information.Information.get_coordinates`

        :param int frame: The timestep from which the coordinates shoule be
                        returned. If ommitted, it defaults to the first
                        frame of the trajectory.
        
        :returns: The set of coordinates from the specified frame.
                    ::

                    [[x1, y1, z1], ... [xn, yn, zn]]

        :rtype: :any:`numpy.array`

        ::

            >>> print(mol.get_coordinates())
            [[ -30.85199928  -81.45800018  365.05499268]
             [ -31.99500084  -80.69300079  365.66900635]
             [ -32.0530014   -81.13200378  367.18200684]
             ...,
             [ -27.54199982  -96.25099945  402.83700562]
             [ -23.54199982  -94.7539978   400.41900635]
             [ -22.86100006  -93.72499847  400.55300903]]

            >>> print(mol.get_coordinates(2))
            [[ -28.88899994  -80.45700073  365.51699829]
             [ -30.20000076  -79.73699951  365.99700928]
             [ -30.90699959  -80.5510025   367.13000488]
             ...,
             [ -26.0189991   -97.28099823  403.52600098]
             [ -23.2140007   -94.73999786  400.94699097]
             [ -22.52899933  -93.73300171  400.81399536]]
        """
        
        return self.information.get_coordinates(frame)

    def get_trajectory_coordinates(self):
        """
        Returns the trajectory for the molecule.

        Wrapper function for :meth:`~scoria.Information.Information.get_trajectory_coordinates`
        
        :returns: The set of all coordinates.
                    ::
                
                        [[[x11, y11, z11], ... [x1n, y1n, z1n]],
                         ...,
                         [[xm1, ym1, zm1], ... [xmn, ymn, zmn]]] 

        :rtype: :any:`numpy.array`

        ::

            >>> for coord in mol.get_trajectory_coordinates():
            >>>     print(coord)
            >>>     print()
            [[ -30.85199928  -81.45800018  365.05499268]
             [ -31.99500084  -80.69300079  365.66900635]
             [ -32.0530014   -81.13200378  367.18200684]
             ..., 
             [ -27.54199982  -96.25099945  402.83700562]
             [ -23.54199982  -94.7539978   400.41900635]
             [ -22.86100006  -93.72499847  400.55300903]]

            [[ -30.6779995   -81.32499695  365.73199463]
             [ -31.88100052  -80.38600159  366.0289917 ]
             [ -32.40399933  -80.62799835  367.45700073]
             ..., 
             [ -27.44400024  -96.71099854  402.64700317]
             [ -23.79199982  -94.58899689  400.63598633]
             [ -23.10700035  -93.56300354  400.79598999]]
             <more>
        """

        return self.information.get_trajectory_coordinates()

    def get_filename(self):
        """
        Returns the filename that the molecule was originally loaded from.

        Wrapper function for :meth:`~scoria.Information.Information.get_filename`
        
        :returns: The name of the file.

        :rtype: :any:`str`

        ::

            >>> mol = scoria.Molecule()
            >>> mol.load_pdb_into("single_frame.pdb")
            >>> print(mol.get_filename())
            single_frame.pdb
        """
        
        return self.information.get_filename()

    def get_remarks(self):
        """
        Returns the remarks from the file the molecule was loaded from.

        Wrapper function for :meth:`~scoria.Information.Information.get_remarks`
        
        :returns: The remarks from the file an a list of strings.

        :rtype: :any:`list`

        ::

            >>> mol = scoria.Molecule()
            >>> mol.load_pdb_into("single_frame.pdb")
            >>> print(mol.get_remarks())
            [' This is a remark.']
        """
        
        return self.information.get_remarks()

    def get_atom_information(self):
        """
        Retreives the atomic information for the molecule.

        Wrapper function for :meth:`~scoria.Information.Information.get_atom_information`

        :returns: A masked array containing the atom information.

        :rtype: :any:`numpy.ma.MaskedArray` 

        The contents of the array are as follows:

        ================ ===== ================= ==============================
        member name      dtype Full Type         Description
        ================ ===== ================= ==============================
        record_name      S6    six char string   What the atom belongs to
        serial           <i8   64-bit integer    The index of the atom
        name             S5    five char string  The atom name
        resname          S5    five char string  The residue name
        chainid          S1    one char string   The chain identifier
        resseq           <i8   64-bit integer    The Residue sequence number   
        occupancy        <f8   64-bit float      Occupancy of atom
        tempfactor       <f8   64-bit float      Tempature Factor
        element          S2    two char string   The element symbol
        charge           S3    three char string Charge on the atom
        name_stripped    S5    five char string  Atom name without space
        resname_stripped S5    five char string  Residue name without space
        chainid_stripped S1    one char string   Chain identifier without space
        element_stripped S2    two char string   Element symbol without space
        ================ ===== ================= ==============================

        An example for printing the elemental symbols of the first five atoms::

            >>> atom_info = mol.get_atom_information()
            >>> print(atom_info['element'][0:5])
            ['N' 'C' 'C' 'O' 'C']
        """

        return self.information.get_atom_information()

    def get_coordinates_undo_point(self):
        """
        NEEDS CLARIFICATION.
        Retreives a previously save set of coordinates to revert to.

        Wrapper function for :meth:`~scoria.Information.Information.get_coordinates_undo_point`

        :returns: A set of coordinates from which to return to.

        :rtype: :any:`numpy.array` or :any:`None`
        """

        return self.information.get_coordinates_undo_point()

    def get_bonds(self):
        """
        Retreives the bonds beteween atoms as a n x n matrix.
        
        Wrapper function for :meth:`~scoria.Information.Information.get_bonds`

        :returns: A binary n x n matrix, where bonds are represented by 1.

        :rtype: :any:`numpy.array`

        An example for finding all atoms bonded with atom 153::

            >>> bonds = mol.get_bonds()
            >>> for i in range(0,len(bonds)):
            ...     if bonds[153][i] == 1:
            ...             print(153,"-",i)
            153 - 152
            153 - 154
            153 - 155
        """
        
        return self.information.get_bonds()

    def get_hierarchy(self):
        """
        NEEDS CLARIFICATION.

        Wrapper function for :meth:`~scoria.Information.Information.get_hierarchy`

        :returns: A dictionary?

        :rtype: :any:`dict`
        """
        
        return self.information.get_hierarchy()

    def get_constants(self):
        """
        Returns a dictionary containing the constants assumed for the molecular model.
        
        Wrapper function for :meth:`~scoria.Information.Information.get_constants`

        :returns: The constants assumed by the model.

        :rtype: :any:`dict`

        ============================== =============== ===============================
        Dictionary Keys                Value Type      Contains
        ============================== =============== ===============================
        mass_dict                      dict{str:float} The mass of elements
        rna_residues                   list(str)       RNA residue names
        f8_fields                      list(str)       Atom Information floats
        vdw_dict                       dict{str:float} Van der Waals force of elements
        i8_fields                      list(str)       Atom Information integers
        protein_residues               list(str)       Protein residue names
        bond_length_dict               dict{str:float} Element-pair bond length 
        element_names_with_two_letters list(str)       Element symbols with 2 letters
        max_number_of_bonds_permitted  dict{str:int}   Max bonds per element
        dna_residues                   list(str)       DNA reside names
        ============================== =============== ===============================

        """
        
        return self.information.get_constants()

    def get_center_of_mass(self, selection = None, frame = None):
        """
        Determines the center of mass.

        Wrapper function for :meth:`~scoria.Information.Information.get_center_of_mass`

        :param numpy.array selection: The indices of
                          the atoms to consider when calculating the center of mass. 
                          If ommitted, all atoms of the scoria.Molecule object 
                          will be considered.

        :param int frame: The timestep at which the center of mass
                      should be calculated. If ommitted, it defaults to the first 
                      frame of the trajectory.
        
        :returns: The x, y, and z coordinates of the center of mass.

        :rtype: :any:`numpy.ma.MaskedArray`

        ::

            >>> mol = scoria.Molecule()
            >>> mol.load_pdb_into("single_frame.pdb")
            >>> print(mol.get_center_of_mass())
            [33.0643089093134 19.135747088722564 16.05629867850796]
        """
        
        return self.information.get_center_of_mass(selection, frame)

    def get_geometric_center(self, selection = None, frame = None):
        """
        Determines the geometric center of the molecule.

        Wrapper function for :meth:`~scoria.Information.Information.get_geometric_center`

        :param numpy.array selection: The indices of
                          the atoms to consider when calculating the geometric.
                          If ommitted, all atoms of the scoria.Molecule object
                          will be considered.

        :param int frame: The timestep at which the geometric center
                      should be calculated. If ommitted, it defaults to the first
                      frame of the trajectory.
        
        :returns: The x, y, and z coordinates of the geometric center.

        :rtype: :any:`numpy.array`

        ::

            >>> mol = scoria.Molecule()
            >>> mol.load_pdb_into("single_frame.pdb")
            >>> print(mol.get_geometric_center())
            [ 33.09860848  19.1221197   16.0426808 ]
        """

        return self.information.get_geometric_center(selection, frame)

    def get_total_mass(self, selection = None):
        """
        Returns the total mass of all atoms within the molecule, or of a given 
        selection. 

        Wrapper function for :meth:`~scoria.Information.Information.get_total_mass`

        :param numpy.array selection: The indices of
                        the atoms to consider when calculating the geometric. 
                        If ommitted, all atoms of the scoria.Molecule object
                        will be considered.

        :returns: The total mass of the atom or selection

        :rtype: :any:`float`

        ::

            >>> print(mol.get_total_mass())
            5289.1729999999998

        """
        
        return self.information.get_total_mass(selection)

    def get_total_number_of_atoms(self, selection = None):
        """
        Counts the number of atoms.
        
        Wrapper function for 
        :meth:`~scoria.Information.Information.get_total_number_of_atoms`

        :param numpy.array selection: An optional numpy.array containing the indices of
                    the atoms to count. If ommitted, all atoms of the
                    scoria.Molecule object will be considered.
        :param int frame: An integer indicating at which timestep the center of
                    mass should be calculated. If ommitted, it defaults to the 
                    first frame of the trajectory.

        :returns:  The total number of atoms.
        :rtype: :any:`int`
        """
        
        return self.information.get_total_number_of_atoms(selection)

    def get_total_number_of_heavy_atoms(self, selection = None):
        """
        Counts the number of heavy atoms (i.e., atoms that are not
        hydrogens). 

        Wrapper function for 
        :meth:`~scoria.Information.Information.get_total_number_of_heavy_atoms`

        :param numpy.array selection: An optional numpy.array containing the indices of
                    the atoms to count. If ommitted, all atoms of the
                    scoria.Molecule object will be considered.

        :returns: The total number of heavy (non-hydrogen) atoms.
        :rtype: :any:`int`
        """
        
        return self.information.get_total_number_of_heavy_atoms(selection)

    def get_bounding_box(self, selection = None, padding = 0.0, frame = None):
        """
        Calculates a box that bounds (encompasses) a set of atoms.

        Wrapper function for :meth:`~scoria.Information.Information.get_bounding_box`

        :param numpy.array selection: An optional numpy.array containing the indices of
                    the atoms to consider. If ommitted, all atoms of the
                    scoria.Molecule object will be considered.
        :param float padding: An optional float. The bounding box will extend this
                    many angstroms beyond the atoms being considered.
        :param int frame: An integer indicating at which timestep the center of
                    mass should be calculated. If ommitted, it defaults to the 
                    first frame of the trajectory.

        :returns: A numpy array representing two 3D points, (min_x, min_y, min_z)
                    and (max_x, max_y, max_z), that bound the molecule.
        :rtype: :any:`numpy.array`
        """
        
        return self.information.get_bounding_box(selection, padding, frame)

    def get_bounding_sphere(self, selection = None, padding = 0.0, frame = None):
        """
        Calculates a sphere that bounds (encompasses) a set of atoms.

        Requires the :any:`numpy` and :any:`scipy<scipy.spatial>` libraries.

        Wrapper function for :meth:`~scoria.Information.Information.get_bounding_sphere`

        :param numpy.array selection: An optional numpy.array containing the indices of
                    the atoms to consider. If ommitted, all atoms of the
                    scoria.Molecule object will be considered.
        :param float padding: An optional float. The bounding sphere will extend
                    this many angstroms beyond the atoms being considered.
        :param int frame: An integer indicating at which timestep the center of
                    mass should be calculated. If ommitted, it defaults to the 
                    first frame of the trajectory.

        :returns: A tuple containing two elements. The first is a numpy.array
                    representing a 3D point, the (x, y, z) center of the
                    sphere. The second is a float, the radius of the sphere.
        :rtype: :any:`tuple` (:any:`numpy.array`, :any:`float`)
        """
        
        return self.information.get_bounding_sphere(selection, padding, frame)

    def get_default_trajectory_frame(self):
        """
        Retreives the default trajectory frame index.

        Wrapper function for :meth:`~scoria.Information.Information.get_default_trajectory_frame`

        :returns: An *int* representing the index of the default trajectory frame.
        """

        return self.information.get_default_trajectory_frame()


    # Set
    def set_filename(self, filename):
        """
        Sets the __filename variable. Note: this does not reload or modify the
        molecule in anyway.
        
        Wrapper function for :meth:`~scoria.Information.Information.set_filename`
         
        :param str filename: String representation of the filename.
        """
        
        self.information.set_filename(filename)

    def set_remarks(self, remarks):
        """
        Sets the __remarks variable.
        
        Wrapper function for :meth:`~scoria.Information.Information.set_remarks`

        :param list(str) remarks: List containing remarks.
        """
        
        self.information.set_remarks(remarks)

    def set_atom_information(self, atom_information):
        """
        Sets the __atom_information variable. See 
        :meth:`~scoria.Molecule.Molecule.get_atom_information` for
        information on the numpy.array structure.
        
        Wrapper function for :meth:`~scoria.Information.Information.set_atom_information`

        :param numpy.array atom_information: An array containing details
                            on the constituent atoms. 
        """
        
        self.information.set_atom_information(atom_information)

    def set_coordinates(self, coordinates, frame = None):
        """
        Sets a specified frame of the __trajectory variable.
        
        Wrapper function for :meth:`~scoria.Information.Information.set_coordinates`
        
        :param numpy.array coordinates: An array of atomic coordinates.
        :param int frame: An integer represeting the frame of the trajectory to be modified
        """
        
        self.information.set_coordinates(coordinates, frame)

    def set_trajectory_coordinates(self, trajectory):
        """
        Sets the __trajectory variable.
        
        Wrapper function for :meth:`~scoria.Information.Information.set_trajectory_coordinates`

        :param numpy.array trajectory: An array of atomic coordinates.
        """

        self.information.set_trajectory_coordinates(trajectory)

    def set_coordinates_undo_point(self, coordinates_undo_point):
        """
        Sets the __coordinates_undo_point variable.
        
        Wrapper function for :meth:`~scoria.Information.Information.set_coordinates_undo_point`
        
        :param numpy.array coordinates_undo_point: A coordinate set to revert 
            to after modification.
        """
        
        self.information.set_coordinates_undo_point(coordinates_undo_point)

    def set_bonds(self, bonds):
        """
        Sets the __bonds variable. See 
        :meth:`~scoria.Molecule.Molecule.get_bonds` for additional 
        information.
        
        Wrapper function for :meth:`~scoria.Information.Information.set_bonds`
        
        :param numpy.array bonds: A binary n x n matrix containing bonding 
            information.
        """
        
        self.information.set_bonds(bonds)

    def set_hierarchy(self, hierarchy):
        """
        DEPRECIATED?
        
        Wrapper function for :meth:`~scoria.Information.Information.set_hierarchy`
        """
        
        self.information.set_hierarchy(hierarchy)

    def set_default_trajectory_frame(self, frame):
        """
        Set's the default trajectory frame for various calculations. 

        Wrapper function for :meth:`~scoria.Information.Information.set_default_trajectory_frame`

        :param int frame: The default frame for coordinate selection.
        """

        self.information.set_default_trajectory_frame(frame)

    # Information functions
    def assign_masses(self):
        """
        Assigns masses to the atoms of the scoria.Molecule object. 

        Wrapper function for :meth:`~scoria.Information.Information.assign_masses`

        ``Note``:
        This will autopopulate the masses according to their element 
        identification and takes no input.
        """
        
        self.information.assign_masses()

    def assign_elements_from_atom_names(self, selection = None):
        """
        Determines the elements of all atoms from the atom names. Note that
        this will overwrite any existing element assignments, including those
        explicitly specified in loaded files. Note that this doesn't populate
        elements_stripped.

        Wrapper function for :meth:`~scoria.Information.Information.assign_elements_from_atom_names`

        :param numpy.array selection: An optional numpy.array containing the indices of
                    the atoms to consider when calculating the center of mass.
                    If ommitted, all atoms of the scoria.Molecule object
                    will be considered.
        """
        
        self.information.assign_elements_from_atom_names(selection)

    def define_molecule_chain_residue_spherical_boundaries(self):
        """
        Identifies spheres that bound (encompass) the entire molecule, the
        chains, and the residues. This information is stored in
        scoria.Information.Information.hierarchy.

        Requires the :any:`numpy` and :any:`scipy<scipy.spatial>` libraries.

        Wrapper function for 
        :meth:`~scoria.Information.Information.define_molecule_chain_residue_spherical_boundaries`
        """
        
        self.information.define_molecule_chain_residue_spherical_boundaries()

    def serial_reindex(self):
        """
        Reindexes the serial field of the atoms in the molecule, starting
        with 1.
        
        Wrapper function for :meth:`~scoria.Information.Information.serial_reindex`
        """
        
        self.information.serial_reindex()

    def resseq_reindex(self):
        """
        Reindexes the resseq field of the atoms in the molecule, starting
        with 1.

        Wrapper function for :meth:`~scoria.Information.Information.resseq_reindex`
        """
        
        self.information.resseq_reindex()
    
    def belongs_to_protein(self, atom_index):
        """
        Checks if the atom is part of a protein. Taken primarily from Amber
        residue names.

        Wrapper function for :meth:`~scoria.Information.Information.belongs_to_protein`
        
        :param int atom_index: An int, the index of the atom to consider.

        :returns: A boolean. True if part of protein, False if not.
        """
        
        return self.information.belongs_to_protein(atom_index)

    def belongs_to_rna(self, atom_index):
        """
        Checks if the atom is part of RNA.

        Wrapper function for :meth:`~scoria.Information.Information.belongs_to_rna`

        :param int atom_index: An int, the index of the atom to consider.

        :returns: A boolean. True if part of rna, False if not.

        """
        
        return self.information.belongs_to_rna(atom_index)

    def belongs_to_dna(self, atom_index):
        """
        Checks if the atom is part of DNA.

        Wrapper function for :meth:`~scoria.Information.Information.belongs_to_dna`

        :param int atom_index: An int, the index of the atom to consider.

        :returns: A boolean. True if part of dna, False if not.
        """
        
        return self.information.belongs_to_dna(atom_index)

    def insert_trajectory_frame(self, index, coordinates):
        """
        Inserts a new coordinate frame at the end of the trajectory.

        Wrapper function for :meth:`~scoria.Information.Information.insert_trajectory_frame`

        :param numpy.array coordinates: A single frame of coordinates to append.
        :param int index: The location where the frame should be added.
        """

        return self.information.insert_trajectory_frame(index, coordinates)
    
    def delete_trajectory_frame(self, index):
        """
        Removes a given frame from the trajectory.
  
        Wrapper function for :meth:`~scoria.Information.Information.delete_trajectory_frame`

        :param int index: Integer of the frame to remove.
        """

        return self.information.delete_trajectory_frame(index)

    def get_trajectory_frame_count(self):
        """
        Returns the number of frames in __trajectory.

        Wrapper function for :meth:`~scoria.Information.Information.get_trajectory_frame_count`

        :returns: The number of frames in the trajectory.
        :rtype: :any:`int`
        """

        return self.information.get_trajectory_frame_count()

    # File I/O class methods
    def load_pym_into(self, filename):
        """
        Loads the molecular data contained in a pym file into the current
        scoria.Molecule object.

        Requires the :any:`numpy` library.

        Wrapper function for :meth:`~scoria.FileIO.FileIO.load_pym_into`

        :param str filename: A string, the filename of the pym file.
        """
        
        self.fileio.load_pym_into(filename)

    def load_pdb_into(self, filename, bonds_by_distance = True,
                      serial_reindex = True, resseq_reindex = False,
                      is_trajectory = False):
        """
        Loads the molecular data contained in a pdb file into the current
        scoria.Molecule object.

        Wrapper function for :meth:`~scoria.FileIO.FileIO.load_pdb_into`

        :param str filename: A string, the filename of the pdb file.
        :param bool bonds_by_distance: An optional boolean, whether or not to
                    determine atomic bonds based on atom proximity. True by
                    default.
        :param bool serial_reindex: An optional boolean, whether or not to
                    reindex the pdb serial field. True by default.
        :param bool resseq_reindex: An optional boolean, whether or not to
                    reindex the pdb resseq field. False by default.
        :param bool is_trajectory: An optional boolean, whether or not the PDB
                    is multi-frame.
        """

        self.fileio.load_pdb_into(
            filename, bonds_by_distance, serial_reindex,
            resseq_reindex, is_trajectory
        )

    def load_pdb_into_using_file_object(self, file_obj,
                                        bonds_by_distance = True,
                                        serial_reindex = True,
                                        resseq_reindex = False,
                                        is_trajectory = False):
        """
        Loads molecular data from a python file object (pdb formatted) into
        the current scoria.Molecule object. Note that most users will want
        to use the load_pdb_into() function instead, which is identical except
        that it accepts a filename string instead of a python file object.

        Requires the :any:`numpy` library.

        Wrapper function for :meth:`~scoria.FileIO.FileIO.load_pdb_into_using_file_object`

        :param file file_obj: A python file object, containing pdb-formatted
                    data.
        :param bool bonds_by_distance: An optional boolean, whether or not to
                    determine atomic bonds based on atom proximity. True by
                    default.
        :param bool serial_reindex: An optional boolean, whether or not to
                    reindex the pdb serial field. True by default.
        :param bool resseq_reindex: An optional boolean, whether or not to
                    reindex the pdb resseq field. False by default.
        :param bool is_trajectory: An optional boolean, whether or not the PDB
                    is multi-frame.
        """

        self.fileio.load_pdb_into_using_file_object(
            file_obj, bonds_by_distance, serial_reindex,
            resseq_reindex, is_trajectory
        )

    def load_pdbqt_into(self, filename, bonds_by_distance = False,
                      serial_reindex = True, resseq_reindex = False,
                      is_trajectory = False):
        """
        Loads the molecular data contained in a pdbqt file into the current
        scoria.Molecule object. Note that this implementation is
        incomplete. It doesn't save atomic charges, for example. The atom
        types are stored in the "element_padded" and "element" columns.

        Wrapper function for :meth:`~scoria.FileIO.FileIO.load_pdbqt_into`

        :param str filename: A string, the filename of the pdbqt file.
        :param bool bonds_by_distance: An optional boolean, whether or not to
                    determine atomic bonds based on atom proximity. False by
                    default, unlike for PDB.
        :param bool serial_reindex: An optional boolean, whether or not to
                    reindex the pdb serial field. True by default.
        :param bool resseq_reindex: An optional boolean, whether or not to
                    reindex the pdbqt resseq field. False by default.
        :param bool is_trajectory: An optional boolean, whether or not the PDB
                    is multi-frame. Defaults of False.
        """

        self.fileio.load_pdbqt_into(
            filename, bonds_by_distance, serial_reindex, resseq_reindex,
            is_trajectory = is_trajectory
        )

    def load_pdbqt_into_using_file_object(self, file_obj,
                                        bonds_by_distance = False,
                                        serial_reindex = True,
                                        resseq_reindex = False,
                                        is_trajectory = False):
        """
        Loads molecular data from a python file object (pdbqt formatted)
        into the current scoria.Molecule object. Note that most users will
        want to use the load_pdb_into() function instead, which is identical
        except that it accepts a filename string instead of a python file
        object.

        Requires the :any:`numpy` library.

        Wrapper function for :meth:`~scoria.FileIO.FileIO.load_pdbqt_into_using_file_object`

        :param file file_obj: A python file object, containing pdb-formatted
                    data.
        :param bool bonds_by_distance: An optional boolean, whether or not to
                    determine atomic bonds based on atom proximity. False by
                    default, unlike for PDB.
        :param bool serial_reindex: An optional boolean, whether or not to
                    reindex the pdb serial field. True by default.
        :param bool resseq_reindex: An optional boolean, whether or not to
                    reindex the pdb resseq field. False by default.
        :param bool is_trajectory: An optional boolean, whether or not the PDB
                    is multi-frame. Defaults of False.
        """

        self.fileio.load_pdbqt_into_using_file_object(
            file_obj, bonds_by_distance, serial_reindex, resseq_reindex,
            is_trajectory = is_trajectory
        )

    def save_pym(self, filename, save_bonds = False, save_filename = False,
                 save_remarks = False, save_hierarchy = False,
                 save_coordinates_undo_point = False):
        """
        Saves the molecular data contained in a scoria.Molecule object
        to a pym file.

        Requires the :any:`numpy` library.

        Wrapper function for :meth:`~scoria.FileIO.FileIO.save_pym`

        :param str filename: An string, the filename to use for saving. (Note
                    that this is actually a directory, not a file.)
        :param bool save_bonds: An optional boolean, whether or not to save
                    information about atomic bonds. False by default.
        :param bool save_filename: An optional boolean, whether or not to save
                    the original (pdb) filename. False by default.
        :param bool save_remarks: An optional boolean, whether or not to save
                    remarks associated with the molecule. False by default.
        :param bool save_hierarchy: An optional boolean, whether or not to save
                    information about spheres the bound (encompass) the whole
                    molecule, the chains, and the residues. False by default.
        :param bool save_coordinates_undo_point: An optional boolean, whether or
                    not to save the last coordinate undo point. False by
                    default.
        """

        self.fileio.save_pym(
            filename, save_bonds, save_filename, save_remarks,
            save_hierarchy, save_coordinates_undo_point
        )

    def save_pdb(self, filename = "", serial_reindex = True,
                 resseq_reindex = False, return_text = False, frame = None):
        """
        Saves the molecular data contained in a scoria.Molecule object
        to a pdb file.

        Wrapper function for :meth:`~scoria.FileIO.FileIO.save_pdb`

        :param str filename: An string, the filename to use for saving.
        :param bool serial_reindex: An optional boolean, whether or not to
                    reindex the pdb serial field. True by default.
        :param bool resseq_reindex: An optional boolean, whether or not to
                    reindex the pdb resseq field. False by default.
        :param bool return_text: An optional boolean, whether or not to return
                    text instead of writing to a file. If True, the filename
                    variable is ignored.
        :param int frame: If specified, a single-frame PDB will be generated.
                    If not specified, a multi-frame PDB will be generated if 
                    the Molecule has multiple frames. Otherwise, the single
                    existing frame will be used.


        :returns: If return_text is True, a PDB-formatted string. Otherwise,
                returns nothing.
        :rtype: :any:`str` or :any:`None`
        """

        return self.fileio.save_pdb(
            filename, serial_reindex, resseq_reindex, return_text, frame
        )

    def load_pdbqt_trajectory_into(self, filename, bonds_by_distance = True,
                                   serial_reindex = True, 
                                   resseq_reindex = False):
        """
        Loads the molecular data contained in a pdbqt trajectoy file (e.g., an
        AutoDock Vina output file) into the current scoria.Molecule
        object.

        Should be called via the wrapper function :meth:`scoria.Molecule.Molecule.load_pdbqt_trajectory_into`

        :param str filename: A string, the filename of the pdbqt file.
        :param bool bonds_by_distance: An optional boolean, whether or not to
                    determine atomic bonds based on atom proximity. True by
                    default.
        :param bool serial_reindex: An optional boolean, whether or not to
                    reindex the pdb serial field. True by default.
        :param bool resseq_reindex: An optional boolean, whether or not to
                    reindex the pdb resseq field. False by default.
        """

        return self.fileio.load_pdbqt_trajectory_into(filename, 
                                             bonds_by_distance, 
                                             serial_reindex, 
                                             resseq_reindex)
        
    def load_pdbqt_trajectory_into_using_file_object(self, file_obj,
                                                     bonds_by_distance = True,
                                                     serial_reindex = True,
                                                     resseq_reindex = False):
        """
        Loads molecular data from a python file object (pdbqt trajectory
        formatted) into the current scoria.Molecule object. Note that most
        users will want to use the load_pdbqt_trajectory_into() function
        instead, which is identical except that it accepts a filename string
        instead of a python file object.

        Wrapper function for
        :meth:`~scoria.FileIO.FileIO.load_pdbqt_trajectory_into_using_file_object`

        :param file file_obj: A python file object, containing pdbqt-formatted
                    trajectory data.
        :param bool bonds_by_distance: An optional boolean, whether or not to
                    determine atomic bonds based on atom proximity. True by
                    default.
        :param bool serial_reindex: An optional boolean, whether or not to
                    reindex the pdb serial field. True by default.
        :param bool resseq_reindex: An optional boolean, whether or not to
                    reindex the pdb resseq field. False by default.
        """

        return self.fileio.load_pdbqt_trajectory_into_using_file_object(file_obj,
                                                     bonds_by_distance,
                                                     serial_reindex,
                                                     resseq_reindex)

    def load_pdb_trajectory_into(self, filename, bonds_by_distance = True,
                                 serial_reindex = True, 
                                 resseq_reindex = False):
        """
        Loads the molecular data contained in a pdb trajectory file into the
        current scoria.Molecule object.

        Should be called via the wrapper function :meth:`scoria.Molecule.Molecule.load_pdb_trajectory_into`

        :param str filename: A string, the filename of the pdb trajectory
                   file.
        :param bool bonds_by_distance: An optional boolean, whether or not to
                    determine atomic bonds based on atom proximity. True by
                    default.
        :param bool serial_reindex: An optional boolean, whether or not to
                    reindex the pdb serial field. True by default.
        :param bool resseq_reindex: An optional boolean, whether or not to
                    reindex the pdb resseq field. False by default.
        """

        return self.fileio.load_pdb_trajectory_into(filename, bonds_by_distance,
                                        serial_reindex, resseq_reindex)
                                
    def load_pdb_trajectory_into_using_file_object(self, file_obj,
                                                   bonds_by_distance = True,
                                                   serial_reindex = True,
                                                   resseq_reindex = False):
        """
        Loads molecular data from a python file object (pdb trajectory
        formatted) into the current scoria.Molecule object. Note that most
        users will want to use the load_pdb_trajectory_into() function
        instead, which is identical except that it accepts a filename string
        instead of a python file object.

        Should be called via the wrapper function :meth:`scoria.Molecule.Molecule.load_pdb_trajectory_into_using_file_object`

        :param file file_obj: A python file object, containing pdb-formatted
                    trajectory data.
        :param bool bonds_by_distance: An optional boolean, whether or not to
                    determine atomic bonds based on atom proximity. True by
                    default.
        :param bool serial_reindex: An optional boolean, whether or not to
                    reindex the pdb serial field. True by default.
        :param bool resseq_reindex: An optional boolean, whether or not to
                    reindex the pdb resseq field. False by default.
        """

        return self.fileio.load_pdb_trajectory_into_using_file_object(file_obj,
                                                                    bonds_by_distance,
                                                                    serial_reindex,
                                                                    resseq_reindex)

    def load_MDAnalysis_into(self, *args):
        """
        Allows import of molecular structure with MDAnalysis
    
        Requires the :any:`MDAnalysis <MDAnalysis.core.AtomGroup>` library.
    
        Wrapper function for 
        :meth:`~scoria.FileIO.FileIO.load_MDAnalysis_into`
         
        :param \*args: Filename, filenames, or list of file names. Used to
            inizalize a MDAnalysis.Universe object.
        """
    
        self.fileio.load_MDAnalysis_into(*args)
    
    
    def load_MDAnalysis_into_using_universe_object(self, universe):
        """
        Allows import of molecular structure with MDAnalysis
    
        Requires the :any:`MDAnalysis <MDAnalysis.core.AtomGroup>` library.
    
        Wrapper function for 
        :meth:`~scoria.FileIO.FileIO.load_MDAnalysis_into_using_universe_object`
         
        :param MDAnalysis.core.Universe universe: MDAnalysis Universe object.
        """
    
        self.fileio.load_MDAnalysis_into_using_universe_object(universe)


    # Atoms and Bonds class methods
    def get_number_of_bond_partners_of_element(self, atom_index, the_element):
        """
        Counts the number of atoms of a given element bonded to a specified
        atom of interest.

        Requires the :any:`numpy` library.
        
        Wrapper function for 
        :meth:`~scoria.AtomsAndBonds.AtomsAndBonds.get_number_of_bond_partners_of_element`

        :param int atom_index: An int, the index of the atom of interest.
        :param str the_element: A string describing the element of the neighbors
                    to be counted.

        :returns: An int, the number of neighboring atoms of the specified
                element.
        :rtype: :any:`int`
        """

        return self.atoms_and_bonds.get_number_of_bond_partners_of_element(
            atom_index, the_element
        )

    def get_index_of_first_bond_partner_of_element(self, atom_index,
                                                   the_element):
        """
        For a given atom of interest, returns the index of the first
        neighbor of a specified element.

        Wrapper function for :meth:`~scoria.AtomsAndBonds.AtomsAndBonds.get_index_of_first_bond_partner_of_element`

        :param int atom_index: An int, the index of the atom of interest.
        :param str the_element: A string specifying the desired element of the
                    neighbor.

        :returns: An int, the index of the first neighbor atom of the specified
                element. If no such neighbor exists, returns -1.
        :rtype: :any:`int`
        """

        return self.atoms_and_bonds.get_index_of_first_bond_partner_of_element(
            atom_index, the_element
        )

    def create_bonds_by_distance(self, remove_old_bond_data = True,
                                 delete_excessive_bonds = True):
        """
        Determines which atoms are bound to each other based on their
        proximity.

        Requires the :any:`numpy` and :any:`scipy<scipy.spatial>` libraries.

        Wrapper function for 
        :meth:`~scoria.AtomsAndBonds.AtomsAndBonds.create_bonds_by_distance`

        :param bool remove_old_bond_data: An optional boolean, whether or not to
                    discard old bond data before adding in bonds determined by
                    distance. True by default.
        :param bool delete_excessive_bonds: An optional boolean, whether or not
                    to check for and delete excessive bonds. True by default.
        """

        self.atoms_and_bonds.create_bonds_by_distance(
            remove_old_bond_data, delete_excessive_bonds
        )

    def delete_bond(self, index1, index2):
        """
        Deletes a bond.

        Wrapper function for :meth:`~scoria.AtomsAndBonds.AtomsAndBonds.delete_bond`
        
        :param int index1: An int, the index of the first atom of the bonded
                    pair.
        :param int index2: An int, the index of the second atom of the bonded
                    pair.
        """
        
        self.atoms_and_bonds.delete_bond(index1, index2)

    def add_bond(self, index1, index2, order = 1):
        """
        Adds a bond.

        Wrapper function for :meth:`~scoria.AtomsAndBonds.AtomsAndBonds.add_bond`

        :param int index1: An int, the index of the first atom of the bonded
                    pair.
        :param int index2: An int, the index of the second atom of the bonded
                    pair.
        :param int order: An optional int, the order of the bond. 1 by default.
        """
        
        self.atoms_and_bonds.add_bond(index1, index2, order)

    def delete_atom(self, index):
        """
        Deletes an atom.

        Wrapper function for :meth:`~scoria.AtomsAndBonds.AtomsAndBonds.delete_atom`

        :param int index: An int, the index of the atom to delete.
        """
        
        self.atoms_and_bonds.delete_atom(index)

    def add_atom(self, record_name = "ATOM", serial = 1, name = "X",
                 resname = "XXX", chainid = "X", resseq = 1, occupancy = 0.0,
                 tempfactor = 0.0, charge = '', element = "X",
                 coordinates = numpy.array([0.0, 0.0, 0.0]), autoindex = True):
        """
        Adds an atom.

        Wrapper function for :meth:`~scoria.AtomsAndBonds.AtomsAndBonds.add_atom`

        :param str record_name: An optional string, the record name of the atom.
                    "ATOM" is the default.
        :param int serial: An optional int, the serial field of the atom. 1 is
                    the default.
        :param str name: An optional string, the name of the atom. 'X' is the
                    default.
        :param str resname: An optional string, the resname of the atom. 'XXX'
                    is the default.
        :param str chainid: An optional string, chainid of the atom. 'X' is the
                    default.
        :param int resseq: An optional int, the resseq field of the atom. 1 is
                    the default.
        :param float occupancy: An optional float, the occupancy of the atom. 0.0
                    is the default.
        :param float tempfactor: An optional float, the tempfactor of the atom.
                    0.0 is the default.
        :param str charge: An optional string, the charge of the atom. '' is the
                    default.
        :param str element: An optional string, the element of the atom. 'X' is
                    the default.
        :param numpy.array coordinates: An optional numpy.array, the (x, y, z)
                    coordinates of the atom. numpy.array([0.0, 0.0, 0.0]) is
                    the default.
        """

        self.atoms_and_bonds.add_atom(
            record_name, serial, name, resname, chainid, resseq, occupancy,
            tempfactor, charge, element, coordinates, autoindex
        )

    # Selections class
    def get_molecule_from_selection(self, selection, serial_reindex = True,
                                    resseq_reindex = False):
        """
        Creates a scoria.Molecule from a user-defined atom selection.

        Wrapper function for :meth:`~scoria.Selections.Selections.get_molecule_from_selection`

        :param numpy.array selection: A numpy.array containing the indices of the atoms
                    in the user-defined selection.
        :param bool serial_reindex: An optional boolean, whether or not to
                    reindex the atom serial fields. Default is True.
        :param bool resseq_reindex: An optional boolean, whether or not to
                    reindex the atom resseq fields. Default is False.

        :returns: A scoria.Molecule object containing the atoms of the
                    user-defined selection.
        """

        return self.selections.get_molecule_from_selection(
            selection, serial_reindex, resseq_reindex
        )

    def select_atoms(self, selection_criteria):
        """
        Select a set of atoms based on user-specified criteria.

        Wrapper function for :meth:`~scoria.Selections.Selections.select_atoms`

        :param dict selection_criteria: A dictionary, where the keys correspond
                    to keys in the
                    self.__parent_Information.Information.get_atom_information()
                    structured numpy array, and the values are lists of
                    acceptable matches. The selection is a logical "AND"
                    between dictionary entries, but "OR" within the value lists
                    themselves. For example: {'atom':['CA', 'O'], 'chain':'A',
                    'resname':'PRO'} would select all atoms with the names CA
                    or O that are located in the PRO residues of chain A.

        :returns: A numpy.array containing the indices of the atoms of the
                    selection.
        """
        
        return self.selections.select_atoms(selection_criteria)

    def select_atoms_in_bounding_box(self, bounding_box):
        """
        Selects all the atoms that are within a bounding box.

        Requires the :any:`numpy` library.

        Wrapper function for :meth:`~scoria.Selections.Selections.select_atoms_in_bounding_box`
    
        :param numpy.array bounding_box: A 2x3 numpy.array containing the minimum and
                    maximum points of the bounding box. Example:
                    numpy.array( [[min_x, min_y, min_z], [max_x, max_y, max_z]] ).

        :returns: A numpy.array containing the indices of the atoms that are
                    within the bounding box.
        """
        
        return self.selections.select_atoms_in_bounding_box(bounding_box)

    def select_branch(self, root_atom_index, directionality_atom_index):
        """
        Identify an isolated "branch" of a molecular model. Assumes the
        atoms with indices root_atom_index and directionality_atom_index are
        bound to one another and that the branch starts at root_atom_index one
        and "points" in the direction of directionality_atom_index.

        Requires the :any:`numpy` library.

        Wrapper function for :meth:`~scoria.Selections.Selections.select_branch`

        :param int root_atom_index: An int, the index of the first atom in the
                branch (the "root").
        :param int directionality_atom_index: An int, the index of the second atom
                in the branch, used to establish directionality

        :returns: A numpy array containing the indices of the atoms of the branch.
        """
        
        return self.selections.select_branch(
            root_atom_index, directionality_atom_index
        )

    def select_all_atoms_bound_to_selection(self, selections):
        """
        Selects all the atoms that are bound to a user-specified selection.

        Requires the :any:`numpy` library.

        Wrapper function for :meth:`~scoria.Selections.Selections.select_all_atoms_bound_to_selection`
        
        :param numpy.array selection: A numpy.array containing the indices of the
                    user-specified selection.

        :returns: A numpy.array containing the indices of the atoms that are
                bound to the user-specified selection. Note that this new
                selection does not necessarily include the indices of the
                original user-specified selection.
        """
        
        return self.selections.select_all_atoms_bound_to_selection(selections)

    def select_atoms_from_same_molecule(self, selection):
        """
        Selects all the atoms that belong to the same molecule as a
        user-defined selection, assuming that the scoria.Molecule object
        actually contains multiple physically distinct molecules that are not
        bound to each other via covalent bonds.

        Requires the :any:`numpy` library.

        Wrapper function for :meth:`~scoria.Selections.Selections.select_atoms_from_same_molecule`
        
        :param numpy.array selection: A numpy.array containing the indices of the
                    user-defined selection.

        :returns: A numpy.array containing the indices of the atoms belonging to
                    the same molecules as the atoms of the user-defined
                    selection.
        """
        
        return self.selections.select_atoms_from_same_molecule(selection)

    def selections_of_constituent_molecules(self):
        """
        Identifies the indices of atoms belonging to separate molecules,
        assuming that the scoria.Molecule object actually contains multiple
        physically distinct molecules that are not bound to each other via
        covalent bonds.

        Requires the :any:`numpy` library.

        Wrapper function for :meth:`~scoria.Selections.Selections.selections_of_constituent_molecules`
        
        :Returns: A python list of numpy.array objects containing the indices of
                    the atoms belonging to each molecule of the composite
                    scoria.Molecule object.
        """
        
        return self.selections.selections_of_constituent_molecules()

    def select_atoms_near_other_selection(self, selection, cutoff):
        """
        Selects all atoms that are near the atoms of a user-defined
        selection.

        Requires the :any:`numpy` and :any:`scipy<scipy.spatial>` libraries.

        Wrapper function for :meth:`~scoria.Selections.Selections.select_atoms_near_other_selection`
        
        :param numpy.array selection: A numpy.array containing the indices of the
                    user-defined selection.
        :param float cutoff: A float, the distance cutoff (in Angstroms).

        :returns: A numpy.array containing the indices of all atoms near the
                    user-defined selection, not including the atoms of the
                    user-defined selection themselves.
        """
        
        return self.selections.select_atoms_near_other_selection(
            selection, cutoff
        )

    def select_atoms_in_same_residue(self, selection):
        """
        Selects all atoms that are in the same residue as any of the atoms
        of a user-defined seleciton. Residues are considered unique if they
        have a unique combination of resname, resseq, and chainid fields.

        Wrapper function for :meth:`~scoria.Selections.Selections.select_atoms_in_same_residue`
        
        :param numpy.array selection: A numpy.array containing the indices of the
                    user-defined selection.

        :returns: A numpy.array containing the indices of all atoms in the same
                    residue as any of the atoms of the user-defined selection.
        """
        
        return self.selections.select_atoms_in_same_residue(selection)

    def invert_selection(self, selection):
        """
        Inverts a user-defined selection (i.e., identifies all atoms that
        are not in the seleciton).

        Wrapper function for :meth:`~scoria.Selections.Selections.invert_selection`
        
        :param numpy.array selection: A numpy.array containing the indices of the
                    user-defined selection.

        :returns: A numpy.array containing the indices of all atoms that are not
                    in the user-defined seleciton.
        """
        
        return self.selections.invert_selection(selection)

    def select_all(self):
        """
        Selects all the atoms in a scoria.Molecule object.

        Wrapper function for :meth:`~scoria.Selections.Selections.select_all`
        
        :returns: A numpy.array containing the indices of all atoms in the
                    scoria.Molecule object.
        """
        
        return self.selections.select_all()

    def select_close_atoms_from_different_molecules(self, other_mol, cutoff,
                                                    pairwise_comparison = True,
                                                    terminate_early = False):
        """
        Effectively detects steric clashes between self and another
        scoria.Molecule.

        Requires the :any:`numpy` and :any:`scipy<scipy.spatial>` libraries.

        Wrapper function for :meth:`~scoria.Selections.Selections.select_close_atoms_from_different_molecules`
        
        :param scoria.Molecule other_mol: A scoria.Molecule object of the other
                    molecule.
        :param float cutoff: A float, the user-defined distance cutoff in
                    Angstroms.
        :param bool pairwise_comparison: An optional boolean, whether or not to
                    perform a simple pairwise distance comparison (if True) or
                    to use a more sophisitcated method (if False). True by
                    default.
        :param bool terminate_early: An optional boolean, whether or not to stop
                    looking for steric clashes once one is found. False by
                    default.

        :returns: A tuple containing two elements. The first is a numpy.array
                    containing the indices of all nearby atoms from this
                    scoria.Molecule object (self). The second is a
                    numpy.array containing the indices of all nearby atoms from
                    the other molecule.
        """

        return self.selections.select_close_atoms_from_different_molecules(
            other_mol, cutoff, pairwise_comparison, terminate_early
        )

    def selections_of_chains(self):
        """
        Identifies the atom selections of each chain.

        Requires the :any:`numpy` library.

        Wrapper function for :meth:`~scoria.Selections.Selections.selections_of_chains`
        
        :returns: A dictionary. The keys of the dictionary correspond to the
                    chainids, and the values are numpy.array objects containing
                    the indices of the associated chain atoms.
        """
        
        return self.selections.selections_of_chains()

    def selections_of_residues(self):
        """
        Identifies the atom selections of each residue.

        Requires the :any:`numpy` library.

        Wrapper function for :meth:`~scoria.Selections.Selections.selections_of_residues`
        
        :returns: A dictionary. The keys of this dictionary correspond to the
                    unique resname-resseq-chainid residue identifiers, and the
                    values are numpy.array objects containing the indices of
                    the associated residue atoms.
        """
        
        return self.selections.selections_of_residues()

    # Manipulation class
    def set_atom_location(self, atom_index, new_location):
        """
        Translates the entire molecular model (without rotating) so that the
        atom with the specified index is located at the specified coordinate.

        Wrapper function for :meth:`~scoria.Manipulation.Manipulation.set_atom_location`
        
        :param int atom_index: An int, the index of the target atom.
        :param numpy.array new_location: A numpy.array specifying the new (x, y, z)
                    coordinate of the specified atom.

        :returns: A numpy.array specifying the (delta_x, delta_y, delta_z) vector
                by which the pmolecule.Molecule was translated.
        """
        
        return self.manipulation.set_atom_location(atom_index, new_location)

    def set_coordinate_undo_point(self):
        """
        Sets ("saves") the undo point of the atom coordinates. Any
        subsequent manipulations of atomic coordinates can be "undone" by
        reseting to this configuration via the coordinate_undo function.

        Wrapper function for :meth:`~scoria.Manipulation.Manipulation.set_coordinate_undo_point`
        """

        self.manipulation.set_coordinate_undo_point()

    def coordinate_undo(self):
        """
        Resets the coordinates of all atoms to those saved using the
        set_coordinate_undo_point function.
        
        Wrapper function for :meth:`~scoria.Manipulation.Manipulation.coordinate_undo`
        """
        
        self.manipulation.coordinate_undo()

    def translate_molecule(self, delta):
        """
        Translate all the atoms of the molecular model by a specified
        vector.

        Wrapper function for :meth:`~scoria.Manipulation.Manipulation.translate_molecule`

        :param numpy.array delta: A numpy.array (delta_x, delta_y, delta_z) specifying the
            amount to move each atom along the x, y, and z coordinates.
        """
        
        self.manipulation.translate_molecule(delta)

    def rotate_molecule_around_a_line_between_points(self, line_point1,
                                                     line_point2, rotate):
        """
        Rotate the molecular model about a line segment. The end points of
        the line segment are explicitly specified coordinates.

        Wrapper function for :meth:`~scoria.Manipulation.Manipulation.rotate_molecule_around_a_line_between_points`
        
        :param numpy.array line_point1: A numpy.array (x, y, z) corresponding to one end
                    of the line segment.
        :param numpy.array line_point2: A numpy.array (x, y, z) corresponding to the
                    other end of the line segment.
        :param float rotate: A float, the angle of rotation, in radians.
        """

        self.manipulation.rotate_molecule_around_a_line_between_points(
            line_point1, line_point2, rotate
        )

    def rotate_molecule_around_a_line_between_atoms(self, line_point1_index,
                                                    line_point2_index, rotate):
        """
        Rotate the molecular model about a line segment. The end points of
        the line segment are atoms of specified indices.

        Wrapper function for :meth:`~scoria.Manipulation.Manipulation.rotate_molecule_around_a_line_between_atoms`

        :param int line_point1_index: An int, the index of the first atom at one
                    end of the line segment.
        :param int line_point2_index: An int, the index of the second atom at
                    the other end of the line segment.
        :param float rotate: A float, the angle of rotation, in radians.
        """

        self.manipulation.rotate_molecule_around_a_line_between_atoms(
            line_point1_index, line_point2_index, rotate
        )

    def rotate_molecule_around_pivot_point(self, pivot, thetax,
                                           thetay, thetaz):
        """
        Rotate the molecular model around a specified atom.

        Requires the :any:`numpy` library.

        Wrapper function for :meth:`~scoria.Manipulation.Manipulation.rotate_molecule_around_pivot_point`
        
        :param numpy.array pivot: A numpy.array, the (x, y, z) coordinate about which
                    the molecular model will be rotated.
        :param float thetax: A float, the angle to rotate relative to the x axis,
                    in radians.
        :param float thetay: A float, the angle to rotate relative to the y axis,
                    in radians.
        :param float thetaz: A float, the angle to rotate relative to the z axis,
                    in radians.
        """

        self.manipulation.rotate_molecule_around_pivot_point(
            pivot, thetax, thetay, thetaz
        )

    def rotate_molecule_around_pivot_atom(self, pivot_index, thetax,
                                          thetay, thetaz):
        """
        Rotate the molecular model around a specified atom.

        Requires the :any:`numpy` library.

        Wrapper function for :meth:`~scoria.Manipulation.Manipulation.rotate_molecule_around_pivot_atom`
        
        :param int pivot_index: An int, the index of the atom about which the
                    molecular model will be rotated.
        :param float thetax: A float, the angle to rotate relative to the x axis,
                    in radians.
        :param float thetay: A float, the angle to rotate relative to the y axis,
                    in radians.
        :param float thetaz: A float, the angle to rotate relative to the z axis,
                    in radians.
        """

        self.manipulation.rotate_molecule_around_pivot_atom(
            pivot_index, thetax, thetay, thetaz
        )

    # Geometry class
    def get_angle_between_three_points(self, pt1, pt2, pt3):
        """
        Computes the angle (in radians) formed by three points (numpy.array
        objects).
            
        Wrapper function for :meth:`~scoria.Geometry.Geometry.get_angle_between_three_points`
        
        :param numpy.array pt1: A numpy.array (x, y, z) representing the first of the
                    three 3D points.
        :param numpy.array pt2: A numpy.array (x, y, z) representing the second of the
                    three 3D points.
        :param numpy.array pt3: A numpy.array (x, y, z) representing the third of the
                    three 3D points.

        :returns: A float containing the angle between the three points, in
                    radians.
        """
        
        return self.geometry.get_angle_between_three_points(pt1, pt2, pt3)

    def get_dihedral_angle(self, pt1, pt2, pt3, pt4):
        """
        Calculates the dihedral angle formed by four points (numpy.array
        objects).

        Wrapper function for :meth:`~scoria.Geometry.Geometry.get_dihedral_angle`
        
        :param numpy.array pt1: A numpy.array (x, y, z) representing the first 3D
                    point.
        :param numpy.array pt2: A numpy.array (x, y, z) representing the second 3D
                    point.
        :param numpy.array pt3: A numpy.array (x, y, z) representing the third 3D
                    point.
        :param numpy.array pt4: A numpy.array (x, y, z) representing the fourth 3D
                    point.

        :returns: A float containing the dihedral angle between the four points,
                    in radians.
        """
        
        return self.geometry.get_dihedral_angle(pt1, pt2, pt3, pt4)

    def get_planarity_deviation(self, pt1, pt2, pt3, pt4):
        """
        Determines how close four points (numpy.array objects) come to lying
        in a common plane.

        Wrapper function for :meth:`~scoria.Geometry.Geometry.get_planarity_deviation`
        
        :param numpy.array pt1: A numpy.array (x, y, z) representing a 3D point.
        :param numpy.array pt2: A numpy.array (x, y, z) representing a 3D point.
        :param numpy.array pt3: A numpy.array (x, y, z) representing a 3D point.
        :param numpy.array pt4: A numpy.array (x, y, z) representing a 3D point.

        :returns: A float, the minimum distance between one point and the plane
                    formed by the other three.
        """
        
        return self.geometry.get_planarity_deviation(pt1, pt2, pt3, pt4)

    def is_planar(self, pt1, pt2, pt3, pt4, planarity_cutoff = 0.2):
        """
        Checks whether four points (numpy.array) lie in a common plane.

        Wrapper function for :meth:`~scoria.Geometry.Geometry.is_planar`
        
        :param numpy.array pt1: A numpy.array (x, y, z) representing a 3D point.
        :param numpy.array pt2: A numpy.array (x, y, z) representing a 3D point.
        :param numpy.array pt3: A numpy.array (x, y, z) representing a 3D point.
        :param numpy.array pt4: A numpy.array (x, y, z) representing a 3D point.
        :param float planarity_cutoff: An optional float. How much the points can
                    deviate (in Angstroms) and still be considered planar. The
                    default is 0.2.

        :returns: A boolean, whether the 4 points can be considered planar.
        """
        
        return self.geometry.is_planar(pt1, pt2, pt3, pt4, planarity_cutoff)

    # Other molecule class
    def get_other_molecules_aligned_to_this(self, other_mol, tethers):
        """
        Aligns a molecule to self (this scoria.Molecule object) using a
        quaternion RMSD alignment.

        Requires the :any:`numpy` library.

        Wrapper function for :meth:`~scoria.OtherMolecules.OtherMolecules.get_other_molecules_aligned_to_this`
                
        :param scoria.Molecule other_mol: A scoria.Molecule that is to be aligned to
                    this one.
        :param tuple tethers: A tuple of two numpy.array objects, where each array
                    contains the indices of self and other_mol, respectively,
                    such that equivalent atoms are listed in the same order.
                    So, for example, if (atom 1, self = atom 3, other) and
                    (atom2, self = atom6, other) than the tethers would be
                    (numpy.array([1, 2]), numpy.array([3, 6])).

        :returns: The new molecule.
        """

        # Add Weight Matrix
        return self.other_molecules.get_other_molecules_aligned_to_this(
            other_mol, tethers
        )

    def get_distance_to_another_molecule(self, other_molecules,
                                         pairwise_comparison = True):
        """
        Computes the minimum distance between any of the atoms of this
        molecular model and any of the atoms of a second specified model.

        Requires the :any:`numpy` and :any:`scipy<scipy.spatial>` libraries.

        Wrapper function for :meth:`~scoria.OtherMolecules.OtherMolecules.get_distance_to_another_molecule`
        
        :param scoria.Molecule other_molecules: a scoria.Molecule, the other molecular
                    model.
        :param bool pairwise_comparison: An optional boolean, whether or not to
                    perform a simple pairwise distance comparison (if True) or
                    to use a more sophisitcated method (if False). True by
                    default.

        :returns: A float, the minimum distance between any two atoms of the two
                specified molecular models (self and other_molecules).
        """

        return self.other_molecules.get_distance_to_another_molecule(
            other_molecules, pairwise_comparison
        )

    def get_distance_to_another_molecules(self, other_molecules,
                                         pairwise_comparison = True):
        """
        DEPRECATION WARNING: Please use :meth:`~scoria.Molecule.Molecule.get_distance_to_another_molecule`

        Computes the minimum distance between any of the atoms of this
        molecular model and any of the atoms of a second specified model.

        Requires the :any:`numpy` and :any:`scipy<scipy.spatial>` libraries.

        Wrapper function for :meth:`~scoria.OtherMolecules.OtherMolecules.get_distance_to_another_molecule`
        
        :param scoria.Molecule other_molecules: a scoria.Molecule, the other molecular
                    model.
        :param bool pairwise_comparison: An optional boolean, whether or not to
                    perform a simple pairwise distance comparison (if True) or
                    to use a more sophisitcated method (if False). True by
                    default.

        :returns: A float, the minimum distance between any two atoms of the two
                specified molecular models (self and other_molecules).
        """

        return self.other_molecules.get_distance_to_another_molecule(
            other_molecules, pairwise_comparison
        )

    def get_rmsd_equivalent_atoms_specified(self, other_mol, tethers):
        """
        Calculates the RMSD between this scoria.Molecle object and
        another, where equivalent atoms are explicitly specified.

        Wrapper function for :meth:`~scoria.OtherMolecules.OtherMolecules.get_rmsd_equivalent_atoms_specified`
        
        :param scoria.Molecule other_mol: The other scoria.Molecule object.
        :param tuple tethers: A tuple of two numpy.array objects, where each array
                    contains the indices of self and other_mol, respectively,
                    such that equivalent atoms are listed in the same order.
                    So, for example, if (atom 1, self = atom 3, other) and
                    (atom2, self = atom6, other) than the tethers would be
                    (numpy.array([1, 2]), numpy.array([3, 6])).

        :returns: A float, the RMSD between self and other_mol.
        """

        return self.other_molecules.get_rmsd_equivalent_atoms_specified(
            other_mol, tethers
        )

    def get_rmsd_order_dependent(self, other_mol):
        """
        Calculates the RMSD between two structures, where equivalent atoms
        are listed in the same order.

        Wrapper function for :meth:`~scoria.OtherMolecules.OtherMolecules.get_rmsd_order_dependent`
        
        :param scoria.Molecule other_mol: The other scoria.Molecule object.

        :returns: A float, the RMSD between self and other_mol.
        """
        
        return self.other_molecules.get_rmsd_order_dependent(other_mol)

    def get_rmsd_heuristic(self, other_mol):
        """
        Caluclates the RMSD between two identical molecules with different
        conformations, per the definition given in "AutoDock Vina: Improving
        the speed and accuracy of docking with a new scoring function,
        efficient optimization, and multithreading,"" by Oleg Trott and Arthur
        J. Olson. Note: Identical means the order of the atoms is the same as
        well.

        Requires the :any:`numpy` library.

        Wrapper function for :meth:`~scoria.OtherMolecules.OtherMolecules.get_rmsd_heuristic`
        
        :param scoria.Molecule other_mol: The other scoria.Molecule object.
            
        :returns: A float, the RMSD between self and other_mol.
        """
        
        return self.other_molecules.get_rmsd_heuristic(other_mol)

    def steric_clash_with_another_molecules(self, other_mol, cutoff,
                                           pairwise_comparison = True):
        """
        DEPRECATION WARNING: Please use :meth:`~scoria.Molecule.Molecule.steric_clash_with_another_molecule`
        Detects steric clashes between the scoria.Molecule (self) and
        another scoria.Molecule.

        Requires the :any:`numpy` and :any:`scipy<scipy.spatial>` libraries.

        Wrapper function for :meth:`~scoria.OtherMolecules.OtherMolecules.steric_clash_with_another_molecule`
        
        :param scoria.Molecule other_mol: The scoria.Molecule object that will be
                    evaluated for steric clashes.
        :param float cutoff: A float, the user-defined distance cutoff in
                    Angstroms.
        :param bool pairwise_comparison: An optional boolean, whether or not to
                    perform a simple pairwise distance comparison (if True) or
                    to use a more sophisitcated method (if False). True by
                    default.

        :returns: A boolean. True if steric clashes are present, False if they
                    are not.
        """

        return self.other_molecules.steric_clash_with_another_molecule(
            other_mol, cutoff, pairwise_comparison
        )

    def steric_clash_with_another_molecule(self, other_mol, cutoff,
                                           pairwise_comparison = True):
        """
        Detects steric clashes between the scoria.Molecule (self) and
        another scoria.Molecule.

        Requires the :any:`numpy` and :any:`scipy<scipy.spatial>` libraries.

        Wrapper function for :meth:`~scoria.OtherMolecules.OtherMolecules.steric_clash_with_another_molecule`
        
        :param scoria.Molecule other_mol: The scoria.Molecule object that will be
                    evaluated for steric clashes.
        :param float cutoff: A float, the user-defined distance cutoff in
                    Angstroms.
        :param bool pairwise_comparison: An optional boolean, whether or not to
                    perform a simple pairwise distance comparison (if True) or
                    to use a more sophisitcated method (if False). True by
                    default.

        :returns: A boolean. True if steric clashes are present, False if they
                    are not.
        """

        return self.other_molecules.steric_clash_with_another_molecule(
            other_mol, cutoff, pairwise_comparison
        )

    def merge_with_another_molecules(self, other_molecules):
        """
        DEPRECATION WARNING: Please use :meth:`~scoria.Molecule.Molecule.merge_with_another_molecule`
        Merges two molecular models into a single model.

        Wrapper function for :meth:`~scoria.OtherMolecules.OtherMolecules.merge_with_another_molecule`
        
        :param scoria.Molecule other_molecules: A molecular model (scoria.Molecule
                    object).

        :returns: A single scoria.Molecule object containing the atoms of
                    this model combined with the atoms of other_molecules.
        """     
        
        return self.other_molecules.merge_with_another_molecule(other_molecules)

    def merge_with_another_molecule(self, other_molecules):
        """
        Merges two molecular models into a single model.

        Wrapper function for :meth:`~scoria.OtherMolecules.OtherMolecules.merge_with_another_molecule`
        
        :param scoria.Molecule other_molecules: A molecular model (scoria.Molecule
                    object).

        :returns: A single scoria.Molecule object containing the atoms of
                    this model combined with the atoms of other_molecules.
        """     
        
        return self.other_molecules.merge_with_another_molecule(other_molecules)

    ######## Supporting functions ########

    def numpy_structured_array_remove_field(self, narray, field_names):
        """
        Removes a specific field name from a structured numpy array.

        :param numpy.array narray: A structured numpy array.
        :param list(str) field_names: A list of strings, where each string is one of
                the field names of narray.

        :returns: A structured numpy array identical to narray, but with the
                field names in field_names removed.
        """

        # surprised this doesn't come with numpy

        # now remove the coordinates from the atom_information object to save
        # memory
        names = list(narray.dtype.names)
        for f in field_names:
            names.remove(f)
        
        return narray[names]

    def __is_number(self, s):
        """
        Determines whether or not a string represents a number.

        :param str s: A string (e.g., "5.4").

        :returns: A boolean, whether or not the string can be represented by a
                    float.
        """

        try:
            float(s)
            return True
        except ValueError:
            return False

    def copy(self):
        """
        Returns an exact copy (scoria.Molecule) of this Molecule object.
        Undo points are NOT copied.

        :returns: A scoria.Molecule, containing to the same atomic
                    information as this scoria.Molecule object.
        """

#        new_molecule = Molecule()
#        new_molecule.set_filename(self.get_filename()[:])
#        new_molecule.set_remarks(self.get_remarks()[:])
#        new_molecule.set_atom_information(self.get_atom_information().copy())
#        new_molecule.set_trajectory_coordinates(self.get_trajectory_coordinates().copy())
#
#        if not self.get_bonds() is None:
#            new_molecule.set_bonds(self.get_bonds().copy())
#        else:
#            new_molecule.set_bonds(None)
#
#        new_molecule.set_hierarchy(copy.deepcopy(self.get_hierarchy()))

        new_molecule = copy.deepcopy(self)

        return new_molecule
