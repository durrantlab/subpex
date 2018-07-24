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
from .six.moves import range
import copy


class Information():
    """
    A class for storing and accessing information about the elements of a
    scoria.Molecule object.
    """

    def __init__(self, parent_molecule_object):
        """
        Initializes the scoria.Information class.

        :param scoria.Molecule parent_molecule_object: The scoria.Molecule object
                    associated with this class.

        """

        self.__parent_molecule = parent_molecule_object

        self.__constants = {}
        # Removed HG from this list to avoid capturing gamma hydrogens
        self.__constants['element_names_with_two_letters'] = ['BR', 'CL', 'BI',
                                                              'AS', 'AG', 'LI',
                                                              'MG', 'RH', 'ZN',
                                                              'MN']

        # SHORTEN LENGTH OF BOND_LENGTH_DICT
        self.__constants['bond_length_dict'] = {
            'C-C': 1.53, 'N-N': 1.425, 'O-O': 1.469, 'S-S': 2.048,
            'C-H': 1.059, 'H-C': 1.059, 'C-N': 1.469, 'N-C': 1.469,
            'C-O': 1.413, 'O-C': 1.413, 'C-S': 1.819, 'S-C': 1.819,
            'N-H': 1.009, 'H-N': 1.009, 'N-O': 1.463, 'O-N': 1.463,
            'O-S': 1.577, 'S-O': 1.577, 'O-H': 0.967, 'H-O': 0.967,
            'S-H': 1.35, 'H-S': 1.35, 'S-N': 1.633, 'N-S': 1.633, 'C-F': 1.399,
            'F-C': 1.399, 'C-CL': 1.790, 'CL-C': 1.790, 'C-BR': 1.910,
            'BR-C': 1.910, 'C-I':2.162, 'I-C':2.162, 'S-BR': 2.321,
            'BR-S': 2.321, 'S-CL': 2.283, 'CL-S': 2.283, 'S-F': 1.640,
            'F-S': 1.640, 'S-I': 2.687, 'I-S': 2.687, 'P-BR': 2.366,
            'BR-P': 2.366, 'P-CL': 2.008, 'CL-P': 2.008, 'P-F': 1.495,
            'F-P': 1.495, 'P-I': 2.490, 'I-P': 2.490, 'P-C': 1.841,
            'C-P': 1.841, 'P-N': 1.730, 'N-P': 1.730, 'P-O': 1.662,
            'O-P': 1.662, 'P-S': 1.954, 'S-P': 1.954, 'N-BR': 1.843,
            'BR-N': 1.843, 'N-CL': 1.743, 'CL-N': 1.743, 'N-F': 1.406,
            'F-N': 1.406, 'N-I': 2.2, 'I-N': 2.2, 'SI-BR': 2.284,
            'BR-SI': 2.284, 'SI-CL': 2.072, 'CL-SI': 2.072, 'SI-F': 1.636,
            'F-SI': 1.636, 'SI-P': 2.264, 'P-SI': 2.264, 'SI-S': 2.145,
            'S-SI': 2.145, 'SI-SI': 2.359, 'SI-SI': 2.359, 'SI-C': 1.888,
            'C-SI': 1.888, 'SI-N': 1.743, 'N-SI': 1.743, 'SI-O': 1.631,
            'O-SI': 1.631, 'X-X': 1.53, 'X-C': 1.53, 'C-X': 1.53, 'X-H': 1.059,
            'H-X': 1.059, 'X-N': 1.469, 'N-X': 1.469, 'X-O': 1.413,
            'O-X': 1.413, 'X-S': 1.819, 'S-X': 1.819, 'X-F': 1.399,
            'F-X': 1.399, 'X-CL': 1.790, 'CL-X': 1.790, 'X-BR': 1.910,
            'BR-X': 1.910, 'X-I':2.162, 'I-X':2.162, 'SI-X': 1.888,
            'X-SI': 1.888, 'H-H':0.74}

        #INCLUDE ALL N AND C TERMINAL

        self.__constants['protein_residues'] = ["ALA", "ARG", "ASH", "ASN",
                                                "ASP", "CYM", "CYS", "CYX",
                                                "GLN", "GLH", "GLU", "GLY",
                                                "HID", "HIE", "HIP", "HIS",
                                                "ILE", "LEU", "LYN", "LYS",
                                                "MET", "PHE", "PRO", "SER",
                                                "THR", "TRP", "TYR", "VAL",
                                                "MSE", "TPO", "PTR", "SEP",
                                                'CALA', 'CARG', 'CASN', 'CASP',
                                                'CCYS', 'CCYX', 'CGLN', 'CGLU',
                                                'CGLY', 'CHID', 'CHIE', 'CHIP',
                                                'CHIS', 'CILE', 'CLEU', 'CLYS',
                                                'CMET', 'CPHE', 'CPRO', 'CSER',
                                                'CTHR', 'CTRP', 'CTYR', 'CVAL',
                                                'NALA', 'NARG', 'NASN', 'NASP',
                                                'NCYS', 'NCYX', 'NGLN', 'NGLU',
                                                'NGLY', 'NHID', 'NHIE', 'NHIP',
                                                'NHIS', 'NILE', 'NLEU', 'NLYS',
                                                'NMET', 'NPHE', 'NPRO', 'NSER',
                                                'NTHR', 'NTRP', 'NTYR', 'NVAL']

        self.__constants['dna_residues'] = ["A", "C", "G", "T", "DA", "DA3",
                                            "DA5", "DAN", "DC", "DC3", "DC4",
                                            "DC5", "DCN", "DG", "DG3", "DG5",
                                            "DGN", "DT", "DT3", "DT5", "DTN"]

        self.__constants['rna_residues'] = ["A", "C", "G", "U", "RA", "RA3",
                                            "RA5", "RAN", "RC", "RC3", "RC4",
                                            "RC5", "RCN", "RG", "RG3", "RG5",
                                            "RGN", "RU", "RU3", "RU5", "RUN"]

        self.__constants['mass_dict'] = {
            'H': 1.00794, 'C': 12.0107, 'N': 14.0067, 'O': 15.9994,
            'S': 32.065, 'P': 30.973762, 'NA': 22.9897, 'MG': 24.3050,
            'F': 18.9984032, 'CL': 35.453, 'K': 39.0983, 'CA': 40.078,
            'I': 126.90447, 'LI': 6.941, 'BE': 9.0122, 'B': 10.811,
            'AL': 26.9815, 'MN': 54.938, 'FE': 55.845, 'CO': 58.9332,
            'CU': 63.9332, 'ZN': 65.38, 'AS': 74.9216, 'BR': 79.904,
            'MO': 95.94, 'RH': 102.9055, 'AG': 107.8682, 'AU': 196.9655,
            'HG': 200.59, 'PB': 207.2, 'BI': 208.98040}

        self.__constants['vdw_dict'] = {
            'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'P': 1.8,
            'S': 1.8, 'B': 2.0, 'LI': 1.82, 'NA': 2.27, 'MG': 1.73, 'AL': 2.00,
            'CL': 1.75, 'CA': 2.00, 'MN': 2.00, 'FE': 2.00, 'CO': 2.00,
            'CU': 1.40, 'ZN': 1.39, 'AS': 1.85, 'BR': 1.85, 'MO': 2.00,
            'RH': 2.00, 'AG': 1.72, 'AU': 1.66, 'PB': 2.02, 'BI': 2.00,
            'K': 2.75, 'I': 1.98
        }

        self.__constants['i8_fields'] = ['serial', 'resseq']
        self.__constants['f8_fields'] = ['x', 'y', 'z', 'occupancy',
                                        'tempfactor']
        self.__constants['max_number_of_bonds_permitted'] = {
            "C": 4, "N": 4, "O": 2, "H": 1, "F": 1, "Cl": 1, "BR": 1, "CL": 1,
            "I": 1, "P": 5, "S": 6
        }

        self.__filename = ""
        self.__remarks = []
        self.__atom_information = None
        self.__trajectory = None
        self.__coordinates_undo_point = None
        self.__bonds = None
        self.__hierarchy = {}
        self.__max_ring_size = 50
        self.__default_frame = 0

    #### Aliases ####
    # Gets
    def get_filename(self):
        """
        Returns the filename that the molecule was originally loaded from.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_filename`
        
        :returns: The names of the files as a list.

        :rtype: :any:`list`

        ::

            >>> mol = scoria.Molecule()
            >>> mol.load_pdb_into("single_frame.pdb")
            >>> print(mol.get_filename())
            single_frame.pdb
        """
        
        #if len(self.__filename) == 1:
        #    return self.__filename[0]
        #else:
        #    return self.__filename

        return self.__filename

    def get_remarks(self):
        """
        Returns the remarks from the file the molecule was loaded from.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_remarks`
        
        :returns: The remarks from the file an a list of strings.

        :rtype: *list*

        ::

            >>> mol = scoria.Molecule()
            >>> mol.load_pdb_into("single_frame.pdb")
            >>> print(mol.get_remarks())
            [' This is a remark.']
        """

        return self.__remarks

    def get_atom_information(self):
        """
        Retreives the atomic information for the molecule.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_atom_information`

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

        return self.__atom_information

    def get_trajectory_coordinates(self):
        """
        Returns the trajectory for the molecule.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_trajectory_coordinates`
        
        :returns: The set of all coordinates.
                    ::
                
                        [[[x11, y11, z11], ... [x1n, y1n, z1n]],
                         ...,
                         [[xm1, ym1, zm1], ... [xmn, ymn, zmn]]] 

        :rtype: *numpy.array*

        ::

            >>> for coord in mol.get_trajectory_coordinates():
            >>>     print(coord)
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

        return copy.deepcopy(self.__trajectory)

    def get_coordinates(self, frame = None):
        """
        Returns the set of coordinates from the specified frame.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_coordinates`

        :param int frame: The timestep from which the coordinates shoule be 
                        returned. If ommitted, it defaults to the first 
                        frame of the trajectory.
        
        :returns: The set of coordinates from the specified frame.
                    ::

                    [[x1, y1, z1], ... [xn, yn, zn]]

        :rtype: *numpy.array*

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

        if frame == None:
            frame = self.__default_frame

        if (self.__trajectory is None) or (len(self.__trajectory) == 0):
            return None
        else:
            return self.__trajectory[frame]

    def get_coordinates_undo_point(self):
        """
        NEEDS CLARIFICATION.
        Retreives a previously save set of coordinates to revert to.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_coordinates_undo_point`

        :returns: A set of coordinates from which to return to.

        :rtype: *numpy.array* or *None*
        """

        return self.__coordinates_undo_point

    def get_bonds(self):
        """
        Retreives the bonds beteween atoms as a n x n matrix.
        
        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_bonds`

        :returns: A binary n x n matrix, where bonds are represented by 1.

        :rtype: *numpy.array*

        An example for finding all atoms bonded with atom 153::

            >>> bonds = mol.get_bonds()
            >>> for i in range(0,len(bonds)):
            ...     if bonds[153][i] == 1:
            ...             print(153,"-",i)
            153 - 152
            153 - 154
            153 - 155
        """

        #if self.__bonds is None:
        #    dim = len(self.get_coordinates())
        #     self.set_bonds(numpy.zeros((dim, dim)))

        return self.__bonds

    def get_hierarchy(self):
        """
        NEEDS CLARIFICATION.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_hierarchy`

        :returns: A dictionary?

        :rtype: *dict*
        """

        if not numpy.class_dependency("create a hierarchical organization of the molecule", "NUMPY"):
            return

        return self.__hierarchy

    def get_constants(self):
        """
        Returns a dictionary containing the constants assumed for the molecular model.
        
        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_constants`

        :returns: The constants assumed by the model.

        :rtype: *dict*

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

        return self.__constants

    # Sets
    def set_filename(self, filename):
        """
        Sets the __filename variable. Note: this does not reload or modify the
        molecule in anyway.
        
        Wrapper function for :meth:`~scoria.Molecule.Molecule.set_filename`
         
        :param str filename: String representation of the filename.
        """

        if type(filename) == str:
            self.__filename = [filename]
        elif type(filename) == list:
            self.__filename = filename

    def set_remarks(self, remarks):
        """
        Sets the __remarks variable.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.set_remarks`

        :param list(str) remarks: List containing remarks.
        """

        self.__remarks = remarks

    def set_atom_information(self, atom_information):
        """
        Sets the __atom_information variable. See
        :meth:`~scoria.Molecule.Molecule.get_atom_information` for
        information on the numpy.array structure.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.set_atom_information`

        :param numpy.array atom_information: An array containing details
                            on the constituent atoms.
        """

        self.__atom_information = atom_information

    def set_trajectory_coordinates(self, trajectory):
        """
        Sets the __trajectory variable.
        
        Wrapper function for :meth:`~scoria.Molecule.Molecule.set_trajectory_coordinates`

        :param numpy.array trajectory: An array of atomic coordinates.
        """

        self.__trajectory = trajectory

    def set_coordinates(self, coordinates, frame = None):
        """
        Sets a specified frame of the __trajectory variable.
        
        Wrapper function for :meth:`~scoria.Molecule.Molecule.set_coordinates`
        
        :param numpy.array coordinates: An array of atomic coordinates.
        :param int frame: An integer represeting the frame of the trajectory to be modified
        """

        if frame == None:
            frame = self.__default_frame

        if (self.__trajectory is None) or (len(self.__trajectory) <= frame):
            self.insert_trajectory_frame(frame, coordinates)
        else:
            self.__trajectory[frame] = coordinates

    def set_coordinates_undo_point(self, coordinates_undo_point):
        """
        Sets the __coordinates_undo_point variable.
        
        Wrapper function for :meth:`~scoria.Molecule.Molecule.set_coordinates_undo_point`
        
        :param numpy.array coordinates_undo_point: A coordinate set to revert 
            to after modification.
        """

        self.__coordinates_undo_point = coordinates_undo_point

    def set_bonds(self, bonds):
        """
        Sets the __bonds variable. See 
        :meth:`~scoria.Molecule.Molecule.get_bonds` for additional 
        information.
        
        Wrapper function for :meth:`~scoria.Molecule.Molecule.set_bonds`
        
        :param numpy.array bonds: A binary n x n matrix containing bonding 
            information.
        """

        self.__bonds = bonds

    def set_hierarchy(self, hierarchy):
        """Sets the __hierarchy variable.
        DEPRECIATED?
        
        Wrapper function for :meth:`~scoria.Molecule.Molecule.set_hierarchy`
        """

        self.__hierarchy = hierarchy

    def belongs_to_protein(self, atom_index):
        """
        Checks if the atom is part of a protein. Taken primarily from Amber
        residue names.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.belongs_to_protein`
        
        :param int atom_index: An int, the index of the atom to consider.

        :returns: A boolean. True if part of protein, False if not.
        """

        # this function is retained for legacy reasons. past versions of
        # scoria had this functionality.

        if (self.__atom_information['resname'][atom_index]
            in self.__constants['protein_residues']):
            return True
        return False

    def belongs_to_dna(self, atom_index):
        """
        Checks if the atom is part of DNA.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.belongs_to_dna`

        :param int atom_index: An int, the index of the atom to consider.

        :returns: A boolean. True if part of dna, False if not.
        """

        # this function is retained for legacy reasons. past versions of
        # scoria had this functionality.

        if (self.__atom_information['resname'][atom_index]
            in self.__constants['dna_residues']):
            return True

        return False

    def belongs_to_rna(self, atom_index):
        """
        Checks if the atom is part of RNA.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.belongs_to_rna`

        :param int atom_index: An int, the index of the atom to consider.

        :returns: A boolean. True if part of rna, False if not.

        """

        # this function is retained for legacy reasons. past versions of
        # scoria had this functionality.

        if (self.__atom_information['resname'][atom_index]
            in self.__constants['rna_residues']):
            return True

        return False

    def assign_masses(self):
        """
        Assigns masses to the atoms of the scoria.Molecule object. 

        Wrapper function for :meth:`~scoria.Molecule.Molecule.assign_masses`

        **Note**:
        This will autopopulate the masses according to their element 
        identification and takes no input.
        """

        if not "mass" in self.__atom_information.dtype.names:
            # only assign if not been assigned previously
            masses = numpy.empty((
                len(self.__atom_information['element'])
            ))

            for i in range(len(self.__atom_information['element'])):
                element = self.__atom_information['element'][i]
                mass = self.__constants['mass_dict'][element]
                masses[i] = mass

            self.__atom_information = numpy.append_fields(self.__atom_information,
                                                    'mass', data = masses)

    def assign_elements_from_atom_names(self, selection = None):
        """
        Determines the elements of all atoms from the atom names. Note that
        this will overwrite any existing element assignments, including those
        explicitly specified in loaded files. Note that this doesn't populate
        elements_stripped.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.assign_elements_from_atom_names`

        :param numpy.array selection: An optional numpy.array containing the indices of
                    the atoms to consider when calculating the center of mass.
                    If ommitted, all atoms of the scoria.Molecule object
                    will be considered.
        """

        if selection is None:
            selection = self.__parent_molecule.select_all()

        if len(selection) == 0:
            return

        # get the atom names
        fix_element_names = copy.deepcopy(self.__atom_information['name_padded']
                                          [selection].astype('<U5'))

        fix_element_names = numpy.defchararray_upper(fix_element_names)

        fix_element_names = numpy.defchararray_strip(fix_element_names)

        # first remove any numbers at the begining of these names
        fix_element_names = numpy.defchararray_lstrip(fix_element_names,
                                                           '0123456789')

        # remove any thing, letters or numbers, that follows a number,
        # including the number itself. so C2L becomes C, not CL.
        for num in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
            # I wish there was a more numpified way of doing this. :(
            tmp = numpy.defchararray_split(fix_element_names.astype('<U5'), num)
            fix_element_names = numpy.empty(len(fix_element_names),
                                            dtype = "S5")
            for i, item in enumerate(tmp):
                fix_element_names[i] = tmp[i][0]

        # take just first two letters of each item
        fix_element_names = numpy.array(fix_element_names, dtype = "|S2")

        # identify ones that are two-letter elements and one-letter elements
        cnsts = self.__constants

        one_tht_shd_b_2_lttrs = [False] * len(fix_element_names)

        for other_two_letter in cnsts['element_names_with_two_letters']:
            one_tht_shd_b_2_lttrs = numpy.logical_or(
                one_tht_shd_b_2_lttrs,
                (fix_element_names == other_two_letter)
            )

        indices_of_two_letter_elements = numpy.nonzero(
            one_tht_shd_b_2_lttrs
        )[0]


        indices_of_one_letter_elements = numpy.nonzero(
            numpy.logical_not(one_tht_shd_b_2_lttrs)
        )[0]


        # get ones that are one-letter elements
        fix_element_names[indices_of_one_letter_elements] = (
            numpy.defchararray_rjust(numpy.array(
                fix_element_names[indices_of_one_letter_elements],
                dtype = "|S1"
            ), 2)
        )

        # they should be capitalized for consistency
        fix_element_names = numpy.defchararray_upper(fix_element_names)
        stripped_element_names = numpy.defchararray_strip(fix_element_names)

        # now map missing element names back
        self.__atom_information['element_padded'][selection] = fix_element_names
        self.__atom_information['element'][selection] = stripped_element_names

        # element_stripped also needs to be updated try:
        # self.__parent_molecule.information.get_atom_information()
        # ['element'][selection] =
        # numpy.defchararray_strip(fix_element_names) except: # so
        # element_stripped hasn't been defined yet
        #    self.__parent_molecule.information.get_atom_information() =
        #    append_fields(self.__parent_molecule.
        #    information.get_atom_information(), 'element',
        #    data = numpy.defchararray_strip(
        #    self.__parent_molecule.information.
        #    get_atom_information()['element_padded']))

    def get_center_of_mass(self, selection = None, frame = None):
        """
        Determines the center of mass.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_center_of_mass`

        :param numpy.array selection: The indices of
                          the atoms to consider when calculating the center of mass.
                          If ommitted, all atoms of the scoria.Molecule object
                          will be considered.

        :param int frame: The timestep at which the center of mass
                      should be calculated. If ommitted, it defaults to the first
                      frame of the trajectory.

        :returns: The x, y, and z coordinates of the center of mass.

        :rtype: *numpy.ma*

        ::

            >>> mol = scoria.Molecule()
            >>> mol.load_pdb_into("single_frame.pdb")
            >>> print(mol.get_center_of_mass())
            [33.0643089093134 19.135747088722564 16.05629867850796]
        """

        if frame == None:
            frame = self.__default_frame

        if selection is None: 
            selection = self.__parent_molecule.select_all()

        # make sure the masses have been asigned
        self.assign_masses()

        # calculate the center of mass

        # multiply each coordinate by its mass
        coors = self.__trajectory[frame][selection]
        masses = numpy.vstack((
            self.__atom_information['mass'][selection],
            self.__atom_information['mass'][selection],
            self.__atom_information['mass'][selection]
        )).T

        center_of_mass = coors * masses

        # now sum all that
        center_of_mass = numpy.sum(center_of_mass, 0)

        # now divide by the total mass
        center_of_mass = center_of_mass / self.get_total_mass(selection)

        return center_of_mass

    def get_geometric_center(self, selection = None, frame = None):
        """
        Determines the geometric center of the molecule.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_geometric_center`

        :param numpy.array selection: The indices of
                          the atoms to consider when calculating the geometric.
                          If ommitted, all atoms of the scoria.Molecule object
                          will be considered.

        :param int frame: The timestep at which the geometric center 
                      should be calculated. If ommitted, it defaults to the first
                      frame of the trajectory.
        
        :returns: The x, y, and z coordinates of the geometric center.

        :rtype: *numpy.array*

        ::

            >>> mol = scoria.Molecule()
            >>> mol.load_pdb_into("single_frame.pdb")
            >>> print(mol.get_geometric_center())
            [ 33.09860848  19.1221197   16.0426808 ]
        """

        if frame == None:
            frame = self.__default_frame

        if selection is None:
            selection = self.__parent_molecule.select_all()

        return (numpy.sum(self.__trajectory[frame][selection], 0) /
                self.get_total_number_of_atoms(selection))

    def get_total_mass(self, selection = None):
        """
        Returns the total mass of all atoms within the molecule, or of a given
        selection.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_total_mass`

        :param numpy.array selection: The indices of
                        the atoms to consider when calculating the geometric.
                        If ommitted, all atoms of the scoria.Molecule object
                        will be considered.

        :returns: The total mass of the atom or selection

        :rtype: *float*

        ::

            >>> print(mol.get_total_mass())
            5289.1729999999998

        """

        if selection is None:
            selection = self.__parent_molecule.select_all()

        # assign masses if necessary
        self.assign_masses()

        # return total mass
        return numpy.sum(self.__atom_information['mass'][selection])

    def get_total_number_of_atoms(self, selection = None, frame = None):
        """
        Counts the number of atoms.
        
        Wrapper function for 
        :meth:`~scoria.Molecule.Molecule.get_total_number_of_atoms`

        :param numpy.array selection: An optional numpy.array containing the indices of
                    the atoms to count. If ommitted, all atoms of the
                    scoria.Molecule object will be considered.
        :param int frame: An integer indicating at which timestep the center of
                    mass should be calculated. If ommitted, it defaults to the 
                    first frame of the trajectory.

        :returns:  The total number of atoms.
        :rtype: *int*
        """

        if frame == None:
            frame = self.__default_frame

        if selection is None:
            selection = self.__parent_molecule.select_all()


        if (self.__trajectory is None) or (len(self.__trajectory) == 0):
            return 0
        else:
            return len(self.__trajectory[frame][selection])

    def get_total_number_of_heavy_atoms(self, selection = None):
        """
        Counts the number of heavy atoms (i.e., atoms that are not
        hydrogens). 

        Wrapper function for 
        :meth:`~scoria.Molecule.Molecule.get_total_number_of_heavy_atoms`

        :param numpy.array selection: An optional numpy.array containing the indices of
                    the atoms to count. If ommitted, all atoms of the
                    scoria.Molecule object will be considered.

        :returns: The total number of heavy (non-hydrogen) atoms.
        :rtype: *int*
        """

        if selection is None:
            selection = self.__parent_molecule.select_all()

        if (self.__trajectory is None) or (len(self.__trajectory) == 0):
            return 0

        all_hydrogens = self.__parent_molecule.select_atoms({
            'element': 'H'
        })

        return self.get_total_number_of_atoms() - len(all_hydrogens)

    def get_bounding_box(self, selection = None, padding = 0.0, frame = None):
        """
        Calculates a box that bounds (encompasses) a set of atoms.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_bounding_box`

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
        :rtype: *numpy.array*
        """

        if frame == None:
            frame = self.__default_frame

        if selection is None: 
            selection = self.__parent_molecule.select_all()

        return numpy.vstack(
            (numpy.min(self.__trajectory[frame][selection], 0) - padding,
             numpy.max(self.__trajectory[frame][selection], 0) + padding))

    def get_bounding_sphere(self, selection = None, padding = 0.0, frame = None):
        """
        Calculates a sphere that bounds (encompasses) a set of atoms.

        Requires the :any:`numpy` and :any:`scipy<scipy.spatial>` libraries.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_bounding_sphere`

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
        :rtype: *tuple* (*numpy.array*, *float*)
        """

        if not numpy.class_dependency("calculate a sphere that bounds a set of atoms", "NUMPY"):
            return

        if not numpy.class_dependency("calculate a sphere that bounds a set of atoms", "SCIPY"):
            return

        if frame == None:
            frame = self.__default_frame

        if selection is None:
            selection = self.__parent_molecule.select_all()

        # get center
        center_of_selection = numpy.array([
            self.get_geometric_center(selection)
        ])

        # get distance to farthest point in selection
        return (center_of_selection[0],
                numpy.max(
                    numpy.cdist(center_of_selection,
                          self.__trajectory[frame][selection])[0])
                )

    def define_molecule_chain_residue_spherical_boundaries(self):
        """
        Identifies spheres that bound (encompass) the entire molecule, the
        chains, and the residues. This information is stored in
        scoria.Molecule.Molecule.hierarchy.

        Requires the :any:`numpy` and :any:`scipy<scipy.spatial>` libraries.

        Wrapper function for 
        :meth:`~scoria.Molecule.Molecule.define_molecule_chain_residue_spherical_boundaries`
        """

        if not numpy.class_dependency("calculate the spherical boundaries around molecules, chains, and residues", "NUMPY"):
            return

        if not numpy.class_dependency("calculate the spherical boundaries around molecules, chains, and residues", "SCIPY"):
            return

        # first, check to see if it's already been defined
        if 'spheres' in list(self.__hierarchy.keys()):
            return

        # set up the new structure
        hrchy = self.__hierarchy
        hrchy['spheres'] = {}
        hrchy['spheres']['molecule'] = {}
        hrchy['spheres']['chains'] = {}
        hrchy['spheres']['residues'] = {}

        # get all the chains and residues
        chains = self.__parent_molecule.selections_of_chains()
        residues = self.__parent_molecule.selections_of_residues()

        # do calcs for the whole molcules
        whole_mol_calc = self.get_bounding_sphere()
        hrchy['spheres']['molecule']['center'] = numpy.array(
            [whole_mol_calc[0]]
        )

        hrchy['spheres']['molecule']['radius'] = whole_mol_calc[1]

        # do calcs for the chains
        hrchy['spheres']['chains']['keys'] = numpy.array(
            # numpy string array e.g. ['a', 'b', 'c']
            list(hrchy['chains']['indices'].keys())
        )

        hrchy['spheres']['chains']['centers'] = numpy.empty(
            (len(hrchy['spheres']['chains']['keys']), 3)
        )

        hrchy['spheres']['chains']['radii'] = numpy.empty(
            len(hrchy['spheres']['chains']['keys'])
        )

        for index, chainid in enumerate(hrchy['spheres']['chains']['keys']):
            asphere = self.get_bounding_sphere(
                selection = hrchy['chains']['indices'][chainid]
            )

            hrchy['spheres']['chains']['centers'][index][0] = asphere[0][0]
            hrchy['spheres']['chains']['centers'][index][1] = asphere[0][1]
            hrchy['spheres']['chains']['centers'][index][2] = asphere[0][2]
            hrchy['spheres']['chains']['radii'][index] = asphere[1]

        # do calcs for the residues
        hrchy['spheres']['residues']['keys'] = numpy.array(
            list(hrchy['residues']['indices'].keys())
        )

        hrchy['spheres']['residues']['centers'] = numpy.empty(
            (len(hrchy['spheres']['residues']['keys']), 3)
        )

        hrchy['spheres']['residues']['radii'] = numpy.empty(
            len(hrchy['spheres']['residues']['keys'])
        )

        for index, resid in enumerate(hrchy['spheres']['residues']['keys']):
            asphere = self.get_bounding_sphere(
                selection = hrchy['residues']['indices'][resid]
            )
            hrchy['spheres']['residues']['centers'][index][0] = asphere[0][0]
            hrchy['spheres']['residues']['centers'][index][1] = asphere[0][1]
            hrchy['spheres']['residues']['centers'][index][2] = asphere[0][2]
            hrchy['spheres']['residues']['radii'][index] = asphere[1]

    def serial_reindex(self):
        """
        Reindexes the serial field of the atoms in the molecule, starting
        with 1.
        
        Wrapper function for :meth:`~scoria.Molecule.Molecule.serial_reindex`
        """

        for i in range(len(self.__atom_information['serial'])):
            self.__atom_information['serial'][i] = i + 1

    def resseq_reindex(self):
        """
        Reindexes the resseq field of the atoms in the molecule, starting
        with 1.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.resseq_reindex`
        """

        keys = numpy.defchararray_add(
            self.__atom_information['resname'], '-'
        )

        keys = numpy.defchararray_add(
            keys,
            numpy.array([str(t) for t in self.__atom_information['resseq']])
        )

        keys = numpy.defchararray_add(keys, '-')

        keys = numpy.defchararray_add(
            keys, self.__atom_information['chainid']
        )

        keys2 = numpy.insert(keys, 0, '')[:-1]
        index_of_change = numpy.nonzero(numpy.logical_not(keys == keys2))[0]
        index_of_change = numpy.append(index_of_change,
                                       len(self.__atom_information))
        

        count = 1
        for t in range(len(index_of_change[:-1])):
            start = index_of_change[t]
            end = index_of_change[t + 1]

            self.__atom_information['resseq'][numpy.arange(
                start, end, 1, dtype = 'int'
            )] = count

            count = count + 1

    def insert_trajectory_frame(self, index, coordinates):
        """
        Inserts a new coordinate frame at the end of the trajectory.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.insert_trajectory_frame`

        :param numpy.array coordinates: A single frame of coordinates to append.
        :param int index: The location where the frame should be added.
        """
        if self.__trajectory is None:
            self.__trajectory = [coordinates]
        else:
            self.__trajectory.insert(index, coordinates)

    def delete_trajectory_frame(self, index):
        """
        Removes a given frame from the trajectory.
  
        Wrapper function for :meth:`~scoria.Molecule.Molecule.delete_trajectory_frame`

        :param int index: Integer of the frame to remove.
        """

        del self.__trajectory[index]

    def get_trajectory_frame_count(self):
        """
        Returns the number of frames in __trajectory.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_trajectory_frame_count`

        :returns: The number of frames in the trajectory.
        :rtype: *int*
        """

        if self.__trajectory is None:
            return 0
        else:
            return len(self.__trajectory)

    def get_default_trajectory_frame(self):
        """
        Retreives the default trajectory frame index.

        :returns: An *int* representing the index of the default trajectory frame.
        """

        return self.__default_frame

    def set_default_trajectory_frame(self, frame):
        """
        Se's the default trajectory frame index for various calculations. 

        :param int frame: The default frame for coordinate selection.
        """

        self.__default_frame = frame