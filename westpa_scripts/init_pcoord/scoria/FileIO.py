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
from scoria import dumbpy as numpy
import os
import sys
from .six.moves import range
from .six.moves import zip

try: from .six.moves import cPickle as pickle  # python2
except: import pickle  # python3

import shutil
import tempfile

try: import cStringIO as StringIO  # python2
except: from io import StringIO  # python3

import scoria

class FileIO():
    """A class for saving and loading molecular data into a scoria.Molecule
    object."""

    def __init__(self, parent_molecule_object):
        """
        Initializes the scoria.FileIO class.

        :param scoria.Molecule parent_molecule_object: The scoria.Molecule object
                    associated with this class.
        """

        self.__parent_molecule = parent_molecule_object
        self.__u = None

    def load_pym_into(self, filename):
        """
        Loads the molecular data contained in a pym file into the current
        scoria.Molecule object.

        Requires the :any:`numpy` library.

        Should be called via the wrapper function :meth:`~scoria.Molecule.Molecule.load_pym_into`

        :param str filename: A string, the filename of the pym file.
        """

        if not numpy.class_dependency("load pym files", "NUMPY"):
            return

        if filename[-1:] != os.sep: filename = filename + os.sep
    
        if numpy.python_version == 2:
            # first, get the files that must exist
            self.__parent_molecule.set_atom_information(
                pickle.load(open(filename + 'atom_information', "rb"))
            )

            self.__parent_molecule.set_coordinates(
                numpy.load(filename + "coordinates.npz")['arr_0']
            )

            # now look for other possible files (optional output)
            prnt = self.__parent_molecule
            if os.path.exists(filename + 'remarks'):
                prnt.set_remarks(pickle.load(open(filename + 'remarks', "rb")))

            if os.path.exists(filename + 'hierarchy'):
                prnt.set_hierarchy(pickle.load(open(filename + 'hierarchy', "rb")))

            if os.path.exists(filename + 'filename'):
                prnt.set_filename(pickle.load(open(filename + 'filename', "rb")))
            
            if prnt.get_filename() == []:  # If still no filename, set it to the one used as a parameter.
                prnt.set_filename(filename)

            if os.path.exists(filename + "bonds.npz"):
                prnt.set_bonds(numpy.load(filename + "bonds.npz")['arr_0'])

            if os.path.exists(filename + "coordinates_undo_point.npz"):
                prnt.set_coordinates_undo_point(
                    numpy.load(filename + "coordinates_undo_point.npz")['arr_0']
                )
        else:
            # first, get the files that must exist
            self.__parent_molecule.set_atom_information(
                pickle.load(open(filename + 'atom_information', "rb"),
                            encoding='latin1')
            )

            self.__parent_molecule.set_coordinates(
                numpy.load(filename + "coordinates.npz")['arr_0']
            )

            # now look for other possible files (optional output)
            prnt = self.__parent_molecule
            if os.path.exists(filename + 'remarks'):
                prnt.set_remarks(pickle.load(open(filename + 'remarks', "rb")))

            if os.path.exists(filename + 'hierarchy'):
                prnt.set_hierarchy(pickle.load(open(filename + 'hierarchy', "rb")))

            if os.path.exists(filename + 'filename'):
                prnt.set_filename(pickle.load(open(filename + 'filename', "rb")))

            if prnt.get_filename() == []:  # If still no filename, set it to the one used as a parameter.
                prnt.set_filename(filename)

            if os.path.exists(filename + "bonds.npz"):
                prnt.set_bonds(numpy.load(filename + "bonds.npz")['arr_0'])

            if os.path.exists(filename + "coordinates_undo_point.npz"):
                prnt.set_coordinates_undo_point(
                    numpy.load(filename + "coordinates_undo_point.npz")['arr_0']
                )


    def load_pdbqt_trajectory_into(self, filename, bonds_by_distance = False,
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

        self.load_pdbqt_into(
            filename, bonds_by_distance, serial_reindex, resseq_reindex,
            is_trajectory = True
        )

    def load_pdbqt_trajectory_into_using_file_object(self, file_obj,
                                                     bonds_by_distance = False,
                                                     serial_reindex = True,
                                                     resseq_reindex = False):
        """
        Loads molecular data from a python file object (pdbqt trajectory
        formatted) into the current scoria.Molecule object. Note that most
        users will want to use the load_pdbqt_trajectory_into() function
        instead, which is identical except that it accepts a filename string
        instead of a python file object.

        Should be called via the wrapper function
        :meth:`scoria.Molecule.Molecule.load_pdbqt_trajectory_into_using_file_object`

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

        self.load_pdbqt_into_using_file_object(
            file_obj, bonds_by_distance, serial_reindex, resseq_reindex,
            is_trajectory = True
        )

    def load_pdbqt_into(self, filename, bonds_by_distance = False,
                      serial_reindex = True, resseq_reindex = False,
                      is_trajectory = False):
        """
        Loads the molecular data contained in a pdbqt file into the current
        scoria.Molecule object. Note that this implementation is
        incomplete. It doesn't save atomic charges, for example. The atom
        types are stored in the "element_padded" and "element" columns.

        Should be called via the wrapper function 
        :meth:`~scoria.Molecule.Molecule.load_pdbqt_into`

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

        self.__parent_molecule.set_filename(filename)

        # open/read the file
        # Treat PDBQT just like PDB, but merge last two columns.
        afile = open(filename, "r")
        self.load_pdbqt_into_using_file_object(
            afile, bonds_by_distance, serial_reindex, resseq_reindex,
            is_trajectory = is_trajectory
        )
        afile.close()

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

        Should be called via the wrapper function
        :meth:`~scoria.Molecule.Molecule.load_pdbqt_into_using_file_object`

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
                    is multi-frame.
        """

        self.load_pdb_into_using_file_object(file_obj, bonds_by_distance,
                                             serial_reindex, resseq_reindex,
                                             is_trajectory = is_trajectory)

        # Now merge the last two columns.
        atom_inf = self.__parent_molecule.get_atom_information()

        # In some instances, the atom information value can be None
        if atom_inf is not None:

            atom_types = numpy.defchararray_add(
                atom_inf["element"], atom_inf["charge"]
            )

            atom_inf["element"] = numpy.defchararray_strip(atom_types)

            atom_inf["charge"] = "\n"

            atom_inf["element_padded"] = numpy.defchararray_rjust(
                atom_inf["element"], 2
            )

        self.__parent_molecule.set_atom_information(atom_inf)

    def load_pdb_trajectory_into(self, filename, bonds_by_distance = False,
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

        self.__parent_molecule.set_filename(filename)

        # open/read the file
        afile = open(filename, "r")
        self.load_pdb_trajectory_into_using_file_object(
            afile, bonds_by_distance, serial_reindex, resseq_reindex
        )
        afile.close()

    def load_pdb_trajectory_into_using_file_object(self, file_obj,
                                                   bonds_by_distance = False,
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

        # Frames are separated by "ENDMDL or "END".
        # A generator to return one frame at a time.
        def get_next_frame(file_obj):
            lines = []
            line = " "

            # Loop through until EOF
            while len(line) != 0:

                # Get the current line
                line = file_obj.readline()

                # If the current line if ENDMDL or END, it's likely a new frame
                # now.
                if line.strip() in ["ENDMDL", "END"]:
                    # print "len(lines)", len(lines)
                    # if len(lines) == 1: print lines

                    # Make sure the list isn't just " "
                    while len(lines) > 0 and lines[0] == " ":
                        lines = lines[1:]

                    if len(lines) > 0:
                        # == 0 if, for example, when ENDMDL then END on next
                        # line.

                        # So turn list into string and yield that value.
                        to_yield = "".join(lines)
                        yield to_yield

                    # Clear list (since new frame) and continue loop.
                    lines = []
                    line = " "
                    continue

                lines.append(line)
            
            # In rare cases, it could be that the PDB file doesn't end in END
            # or ENDMDL. so consider any remaining lines in the lines list.
            while len(lines) > 0 and (lines[0] == " " or lines[0] == ""):
                lines = lines[1:]
            if len(lines) > 0:
                yield "".join(lines)
        
        # Get the information from each frame and place it in the current Molecule object.
        # Note that this could be parallelized in the future...
        first_line = True
        trajectoryList = []
        for pdb_frame in get_next_frame(file_obj):
            if numpy.python_version == 2:
                str_file_obj = StringIO.StringIO(pdb_frame)
            else:
                str_file_obj = StringIO(pdb_frame)

            if first_line == True:
                # First frame, load it into the current molecule
                first_line = False
                self.load_pdb_into_using_file_object(
                    str_file_obj, bonds_by_distance, serial_reindex,
                    resseq_reindex, is_trajectory = False
                )
                trajectoryList.append(self.__parent_molecule.get_coordinates())
            else:
                # subsequent frames, load it into a tmp molecule and copy over
                # the coordinates.
                tmp_mol = scoria.Molecule()
                tmp_mol.load_pdb_into_using_file_object(
                    str_file_obj, bonds_by_distance = False, 
                    serial_reindex = False, resseq_reindex = False,
                    is_trajectory = False
                )
                
                trajectoryList.append(tmp_mol.get_coordinates())

        self.__parent_molecule.set_trajectory_coordinates(trajectoryList)        

    def load_pdb_into(self, filename, bonds_by_distance = False,
                      serial_reindex = True, resseq_reindex = False,
                      is_trajectory = False):
        """
        Loads the molecular data contained in a pdb file into the current
        scoria.Molecule object.

        Should be called via the wrapper function :meth:`~scoria.Molecule.Molecule.load_pdb_into`

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

        self.__parent_molecule.set_filename(filename)

        # open/read the file
        afile = open(filename, "r")
        self.load_pdb_into_using_file_object(afile, bonds_by_distance,
                                             serial_reindex, resseq_reindex,
                                             is_trajectory)
        afile.close()

    def load_pdb_into_using_file_object(self, file_obj,
                                        bonds_by_distance = False,
                                        serial_reindex = True,
                                        resseq_reindex = False,
                                        is_trajectory = False):
        """
        Loads molecular data from a python file object (pdb formatted) into
        the current scoria.Molecule object. Note that most users will want
        to use the load_pdb_into() function instead, which is identical except
        that it accepts a filename string instead of a python file object.

        Should be called via the wrapper function 
        :meth:`~scoria.Molecule.Molecule.load_pdb_into_using_file_object`

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

        if is_trajectory == True:
            self.load_pdb_trajectory_into_using_file_object(
                file_obj, bonds_by_distance = bonds_by_distance,
                serial_reindex = serial_reindex, resseq_reindex =
                resseq_reindex
            )
            return

        # source_data = numpy.genfromtxt(file_obj,
        # dtype="S6,S5,S5,S4,S2,S4,S4,S8,S8,S8,S6,S6,S10,S2,S2",
        # names=['record_name', 'serial', 'name_padded', 'resname_padded', 'chainid_padded',
        # 'resseq', 'empty', 'x', 'y', 'z', 'occupancy', 'tempfactor',
        # 'empty2', 'element_padded', 'charge'], delimiter=[6, 5, 5, 4, 2, 4, 4, 8, 8,
        # 8, 6, 6, 10, 2, 2])


        if numpy.python_version == 2:
            source_data = numpy.genfromtxt(
                file_obj,
                dtype = "S6,S5,S5,S5,S1,S4,S4,S8,S8,S8,S6,S6,S10,S2,S3",
                names = ['record_name', 'serial', 'name_padded', 'resname_padded', 'chainid_padded',
                    'resseq', 'empty', 'x', 'y', 'z', 'occupancy',
                    'tempfactor', 'empty2', 'element_padded', 'charge'],
                delimiter = [6, 5, 5, 5, 1, 4, 4, 8, 8, 8, 6, 6, 10, 2, 3]
            )
        elif numpy.python_version == 3:
            data = []
            for line in file_obj:
                toappend = (line[:6], line[6:11], line[11:16], line[16:21], line[21:22], line[22:26], line[26:30], line[30:38], line[38:46], line[46:54], line[54:60], line[60:66], line[66:76], line[76:78], line[78:81])
                #toappend = [s.decode("utf-8") for s in toappend]
                data.append(toappend)

            #s1 = ["|S6","|S5","|S5","|S5","|S1","|S4","|S4","|S8","|S8","|S8","|S6","|S6","|S10","|S2","|S3"]
            s1 = ["|U6","|U5","|U5","|U5","|U1","|U4","|U4","|U8","|U8","|U8","|U6","|U6","|U10","|U2","|U3"]
            names = ['record_name', 'serial', 'name_padded', 'resname_padded', 'chainid_padded', 'resseq', 'empty', 'x', 'y', 'z', 'occupancy', 'tempfactor', 'empty2', 'element_padded', 'charge']
            dtype = [l for l in zip(names, s1)]
            source_data = numpy.array(data, dtype=dtype)

        # get the remarks, if any. good to hold on to this because some of my
        # programs might retain info via remarks
        remark_indices = numpy.nonzero(
            source_data['record_name'] == "REMARK"
        )[0]

        remarks = []
        for index in remark_indices:
            astr = ""

            for name in source_data.dtype.names[1:]:
                astr = astr + source_data[name][index]
            remarks.append(astr.rstrip())

        self.__parent_molecule.set_remarks(remarks)

        # in case the pdb file has only one line
        if source_data.ndim == 0:
            source_data = source_data.reshape(1, -1)

        # get the ones that are ATOM or HETATOM in the record_name
        or_matrix = numpy.logical_or((source_data['record_name'] == "ATOM  "),
                                     (source_data['record_name'] == "HETATM"))
        indices_of_atom_or_hetatom = numpy.nonzero(or_matrix)[0]

        self.__parent_molecule.set_atom_information(
            source_data[indices_of_atom_or_hetatom]
        )

        # now, some of the data needs to change types
        # first, fields that should be numbers cannot be empty strings
        atom_inf = self.__parent_molecule.get_atom_information()
        for field in (self.__parent_molecule.get_constants()['i8_fields'] +
                      self.__parent_molecule.get_constants()['f8_fields']):
            check_fields = atom_inf[field]
            check_fields = numpy.defchararray_strip(check_fields)
            indices_of_empty = numpy.nonzero(check_fields == '')[0]
            atom_inf[field][indices_of_empty] = '0'

        # now actually change the type
        old_types = atom_inf.dtype

        descr = old_types.descr

        for field in self.__parent_molecule.get_constants()['i8_fields']:
            index = atom_inf.dtype.names.index(field)
            descr[index] = (descr[index][0], 'i8')
        for field in self.__parent_molecule.get_constants()['f8_fields']:
            index = atom_inf.dtype.names.index(field)
            descr[index] = (descr[index][0], 'f8')
        # You need to create this descr object. strings are prefixed with |,
        # and int and float with <
        new_types = numpy.dtype(descr)
        self.__parent_molecule.set_atom_information(atom_inf.astype(new_types))

        # remove some of the fields that just contain empty data
        atom_inf = self.__parent_molecule.get_atom_information()
        self.__parent_molecule.set_atom_information(
            self.__parent_molecule.numpy_structured_array_remove_field(
                atom_inf, ['empty', 'empty2']
            )
        )

        # the coordinates need to be placed in their own special numpy array to
        # facilitate later manipulation
        atom_inf = self.__parent_molecule.get_atom_information()

        self.__parent_molecule.set_coordinates(
            numpy.vstack([atom_inf['x'], atom_inf['y'], atom_inf['z']]).T
            )

        # now remove the coordinates from the atom_information object to save
        # memory
        self.__parent_molecule.set_atom_information(
            self.__parent_molecule.numpy_structured_array_remove_field(
                atom_inf, ['x', 'y', 'z']
            )
        )

        # string values in
        # self.__parent_molecule.information.get_atom_information() should also
        # be provided in stripped format for easier comparison
        fields_to_strip = ['name', 'resname', 'chainid', 'element']
        for f in fields_to_strip:
            self.__parent_molecule.set_atom_information(
                numpy.append_fields(
                    self.__parent_molecule.get_atom_information().copy(),
                    f,
                    data = numpy.defchararray_strip(atom_inf[f + '_padded']),
                    usemask=False
                )
            )

        # now determine element from atom name for those entries where it's not
        # given note that the
        # molecule.information.assign_elements_from_atom_names function can be
        # used to overwrite this and assign elements based on the atom name
        # only.
        indicies_where_element_is_not_defined = numpy.nonzero(
            numpy.defchararray_strip(atom_inf['element_padded']) == ''
        )[0]

        self.__parent_molecule.assign_elements_from_atom_names(
            indicies_where_element_is_not_defined
        )

        # now, if there's conect data, load it. this part of the code is not
        # that "numpyic"
        conect_indices = numpy.nonzero(
            source_data['record_name'] == "CONECT"
        )[0]

        if len(conect_indices) > 0:

            self.__parent_molecule.set_bonds(numpy.zeros((len(atom_inf),
                                                          len(atom_inf))))

            # build serial to index mapping
            serial_to_index = {}
            # is there a faster way?
            for index, inf in enumerate(atom_inf['serial']):
                serial_to_index[inf] = index

            # get the connect data
            for index in conect_indices:
                astr = ""
                for name in source_data.dtype.names[1:]:
                    astr = astr + source_data[name][index]
                astr = astr.rstrip()

                indices = []
                for i in range(0, len(astr), 5):
                    indices.append(serial_to_index[int(astr[i:i + 5])])

                for partner_index in indices[1:]:
                    self.__parent_molecule.add_bond(indices[0], partner_index)

        # else: # create empty bond array
        #    self.__parent_molecule.information.get_bonds() =
        #    numpy.zeros((len(self.__parent_molecule.information.
        #        get_atom_information()),
        #    len(self.__parent_molecule.information.get_atom_information())))

        if bonds_by_distance == True:
            self.__parent_molecule.create_bonds_by_distance(True)

        if serial_reindex == True:
            self.__parent_molecule.serial_reindex()

        if resseq_reindex == True:
            self.__parent_molecule.resseq_reindex()

    def save_pym(self, filename, save_bonds = False, save_filename = False,
                 save_remarks = False, save_hierarchy = False,
                 save_coordinates_undo_point = False):
        """
        Saves the molecular data contained in a scoria.Molecule object
        to a pym file.

        Requires the :any:`numpy` library.

        Should be called via the wrapper function 
        :meth:`~scoria.Molecule.Molecule.save_pym`

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

        if not numpy.class_dependency("save pym files", "NUMPY"):
            return

        # Why not just pickle self.parent.information? Because it's a huge
        # file, can't selectively not save bonds, for example, and numpy.save
        # is faster than cPickle protocol 2 on numpy arrays

        # if the directory already exists, first delete it
        if os.path.exists(filename):
            try:
                shutil.rmtree(filename)
            except:
                pass

            # it could be a file, not a directory
            try:
                os.remove(filename)
            except:
                pass

        # filename is actually a directory, so append separator if needed
        if filename[-1:] != os.sep:
            filename = filename + os.sep

        # make directory
        os.mkdir(filename)

        # save components

        # python objects must be pickled
        if save_hierarchy == True:
            # note this is a combo of python objects and numpy arrays, so must
            # be pickled.
            pickle.dump(self.__parent_molecule.get_hierarchy(),
                        open(filename + 'hierarchy', 'wb'), -1)

        if save_remarks == True:
            # using the latest protocol
            pickle.dump(self.__parent_molecule.get_remarks(),
                        open(filename + 'remarks', 'wb'), -1)

        if save_filename == True:
            pickle.dump(self.__parent_molecule.get_filename(),
                        open(filename + 'filename', 'wb'), -1)

        # unfortunately, the speedy numpy.save doesn't work on masked arrays
        # masked arrays have a dump method, but it just uses cPickle so we're
        # just going to cPickle masked arrays. Could be so much faster if numpy
        # were up to speed... :( not clear that numpy.ma.dump accepts protocol
        # parameter, so let's just use cPickle directly
        pickle.dump(self.__parent_molecule.get_atom_information(),
                    open(filename + 'atom_information', 'wb'), -1)

        # fortunately, coordinates and bonds are regular numpy arrays they can
        # be saved with numpy's speedy numpy.save function note that I'm
        # compressing them here. benchmarking suggests this takes longer to
        # save, but is much faster to load. so I'm prioritizing load times over
        # save times note also that numpy.savez can save multiple arrays to a
        # single file, probably speeding up load.

        numpy.savez(filename + "coordinates.npz",
                    self.__parent_molecule.get_coordinates())

        if save_bonds == True:
            numpy.savez(filename + "bonds.npz",
                        self.__parent_molecule.get_bonds())

        if save_coordinates_undo_point == True:
            numpy.savez(filename + "coordinates_undo_point.npz",
                        self.__parent_molecule.get_coordinates_undo_point())

    def save_pdb(self, filename = "", serial_reindex = True, resseq_reindex = False,
                 return_text = False, frame = None):
        """
        Saves the molecular data contained in a scoria.Molecule object
        to a pdb file.

        Should be called via the wrapper function :meth:`~scoria.Molecule.Molecule.save_pdb`

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
        :rtype: *str* or *None*
        """

        if len(self.__parent_molecule.get_atom_information()) > 0:
            # so the pdb is not empty (if it is empty, don't save)

            # Quick check if it has more than one frame. If so, switch to a
            # multi-frame output function.
            if (self.__parent_molecule.get_trajectory_frame_count() > 1 and 
                frame is None):
                # It has multiple frames and the user hasn't specified a
                # specific frame, so generate the whole trajectory.
                return self._save_pdb_trajectory(filename, serial_reindex, 
                                                 resseq_reindex, return_text)

            # If you get to this point, set the frame to 0.
            if frame is None:
                frame = self.__parent_molecule.get_default_trajectory_frame()

            if serial_reindex == True: self.__parent_molecule.serial_reindex()
            if resseq_reindex == True: self.__parent_molecule.resseq_reindex()

            if return_text == False:
                afile = open(filename, "w")
            else:
                return_string = ""

            # print out remarks
            for line in self.__parent_molecule.get_remarks():
                remark = "REMARK" + line + "\n"

                if return_text == False:
                    afile.write(remark)
                else:
                    return_string += remark

            # print out coordinates
            atom_information = self.__parent_molecule.get_atom_information()
            coordinates = self.__parent_molecule.get_coordinates(frame)

            #if numpy.python_version == 2: 
            dtype_to_use = '|S5'
            #else: 
            #    dtype_to_use = '|U5'  # python3 needs this instead

            printout = numpy.defchararray_add(
                atom_information['record_name'].astype('|S6'),
                numpy.defchararray_rjust(
                    atom_information['serial'].astype('|S5'), 5
                )
            )

            printout = numpy.defchararray_add(printout,
                                              atom_information['name_padded'].astype('|S5'))

            printout = numpy.defchararray_add(printout,
                                              atom_information['resname_padded'].astype('|S5'))

            printout = numpy.defchararray_add(printout,
                                              atom_information['chainid_padded'].astype('|S1'))

            #if numpy.python_version == 2: 
            dtype_to_use = '|S4'
            #else: dtype_to_use = '|U4'  # python3 needs this instead

            printout = numpy.defchararray_add(
                printout, numpy.defchararray_rjust(
                    atom_information['resseq'].astype(dtype_to_use), 4
                )
            )

            printout = numpy.defchararray_add(printout, '    '.encode())

            dtype_to_use = '|S7'

            printout = numpy.defchararray_add(
                printout, numpy.defchararray_rjust(
                    numpy.array(
                        ["%.3f" % t for t in numpy.get_col(coordinates, 0)]
                        ).astype(dtype_to_use), 8
                )
            )

            printout = numpy.defchararray_add(
                printout, numpy.defchararray_rjust(
                    numpy.array(
                        ["%.3f" % t for t in numpy.get_col(coordinates, 1)]
                        ).astype(dtype_to_use), 8
                )
            )

            printout = numpy.defchararray_add(
                printout, numpy.defchararray_rjust(
                    numpy.array(["%.3f" % t
                        for t in numpy.get_col(coordinates, 2)]
                        ).astype(dtype_to_use), 8
                )
            )

           #dtype_to_use = '|S5'

            printout = numpy.defchararray_add(
                printout, numpy.defchararray_rjust(
                    numpy.array(["%.2f" % t
                        for t in atom_information['occupancy']]
                        ).astype(dtype_to_use), 6
                )
            )

            printout = numpy.defchararray_add(
                printout, numpy.defchararray_rjust(
                    numpy.array(["%.2f" % t
                        for t in atom_information['tempfactor']]
                        ).astype(dtype_to_use), 6
                )
            )

            printout = numpy.defchararray_add(printout, '          '.encode())

            dtype_to_use = '|S2'

            printout = numpy.defchararray_add(
                printout, atom_information['element_padded'].astype(dtype_to_use)
            )

            dtype_to_use = '|S3'

            printout = numpy.defchararray_add(
                printout, atom_information['charge'].astype(dtype_to_use)
            )

            printout_string = []
            for i in printout:
                printout_string.append(i.decode('utf-8'))

            if return_text == False:
                if printout_string[0][-1:] == "\n":
                    afile.write("".join(printout_string) + "\n")
                else:
                    afile.write("\n".join(printout_string) + "\n")
            else:
                if printout_string[0][-1:] == "\n":
                    return_string += "".join(printout) + "\n"
                else:
                    return_string += "\n".join(printout) + "\n"

            # print out connect
            prnt = self.__parent_molecule
            atm_inf = prnt.get_atom_information()
            sel_atms_bnd_to_sel = prnt.select_all_atoms_bound_to_selection
            if not prnt.get_bonds() is None:
                for indx in range(len(prnt.get_bonds())):
                    indices_of_bond_partners = sel_atms_bnd_to_sel(
                        numpy.array([indx])
                    )

                    if len(indices_of_bond_partners) > 0:

                        if return_text == False:
                            afile.write(
                                "CONECT" +
                                str(atm_inf["serial"][indx]).rjust(5) +
                                "".join([str(atm_inf["serial"][t]).rjust(5)
                                         for t in indices_of_bond_partners]) +
                                "\n"
                            )
                        else:
                            return_string += (
                                "CONECT" +
                                str(atm_inf["serial"][indx]).rjust(5) +
                                "".join([str(atm_inf["serial"][t]).rjust(5)
                                         for t in indices_of_bond_partners]) +
                                "\n"
                            )

            if return_text == False:
                afile.close()
            else:
                return return_string

        else:
            print(("ERROR: Cannot save a Molecule with no atoms " +
                  "(file name \"" + filename + "\")"))

    def _save_pdb_trajectory(self, filename = "", serial_reindex = True,
                 resseq_reindex = False, return_text = False):
        """
        Saves the molecular trajectory data contained in a scoria.Molecule
        object to a pdb file.

        Should be called via the wrapper function :meth:`scoria.Molecule.Molecule.save_pdb`

        :param str filename: An string, the filename to use for saving.
        :param bool serial_reindex: An optional boolean, whether or not to
                    reindex the pdb serial field. True by default.
        :param bool resseq_reindex: An optional boolean, whether or not to
                    reindex the pdb resseq field. False by default.
        :param bool return_text: An optional boolean, whether or not to return
                    text instead of writing to a file. If True, the filename
                    variable is ignored.

        :returns: If return_text is True, a PDB-formatted string. Otherwise,
                returns nothing.
        :rtype: *str* or *None*
        """

        num_frames = self.__parent_molecule.get_trajectory_frame_count()
        if num_frames > 1:
            # Only proceed if it's a multi-frame trajectory.

            if return_text == True:
                all_txt = ""
            else:
                out = open(filename, 'w')
    

            for frame in range(num_frames):
                pdb_txt = self.save_pdb(
                    filename, serial_reindex, resseq_reindex, 
                    True, frame = frame
                )

                # Remove the CONECT and REMARK information. I don't know how
                # to deal with that when it's multiple frames.
                pdb_txt = "\n".join(l for l in pdb_txt.split("\n") if not l.startswith("CONECT") and not l.startswith("REMARK"))
                
                # Add start and end tags
                pdb_txt = "MODEL" + str(frame + 1).rjust(9) + "\n" + pdb_txt.strip() + "\nENDMDL\n"

                if return_text == True:
                    all_txt = all_txt + pdb_txt
                else:
                    out.write(pdb_txt)

            if return_text == True:
                all_txt = all_txt + "END"
                return all_txt
            else:
                out.write("END")
                out.close()
                return

    def load_MDAnalysis_into(self, *args):
        """
        ***Note this function is only functional in Scoria_MDA***

        Allows import of molecular structure with MDAnalysis.
    
        Requires the :any:`MDAnalysis <MDAnalysis.core.AtomGroup>` library.
    
        Should be called via the wrapper function
        :meth:`~scoria.Molecule.Molecule.load_MDAnalysis_into`
    
        :params \*args: Filename, filenames, or list of file names. Used to
            inizalize a MDAnalysis.Universe object.
    
        """
        print("Compatibility with MDAnalysis objects is not supported in this version.")
    
    def load_MDAnalysis_into_using_universe_object(self, universe):
        """
        ***Note this function is only functional in Scoria_MDA***
        
        Allows import of molecular structure from an MDAnalysis object.
    
        Requires the :any:`MDAnalysis <MDAnalysis.core.AtomGroup>` library.
    
        Should be called via the wrapper function
        :meth:`~scoria.Molecule.Molecule.load_MDAnalysis_into`
    
        :param mdanalysis.universe universe: An MDAnalysis universe object to
            import.
        """
        print("Compatibility with MDAnalysis objects is not supported in this version.")
