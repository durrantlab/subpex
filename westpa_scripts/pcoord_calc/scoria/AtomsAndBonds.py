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
from .six.moves import range


class AtomsAndBonds():
    """
    A class for adding and deleting atoms and bonds. Subclass to the
    :py:class:`scoria.Molecule` class.
    """

    def __init__(self, parent_molecule_object):
        """
        Initializes the scoria.AtomsAndBonds class.

        :param scoria.Molecule parent_molecule_object: The scoria.Molecule object
                    associated with this class.

        """

        self.__parent_molecule = parent_molecule_object

    def create_bonds_by_distance(self, remove_old_bond_data = True,
                                 delete_excessive_bonds = True):
        """
        Determines which atoms are bound to each other based on their
        proximity.

        Requires the :any:`numpy` and :any:`scipy<scipy.spatial>` libraries.

        Should be called via the wrapper function
        :meth:`~scoria.Molecule.Molecule.create_bonds_by_distance`.

        :param bool remove_old_bond_data: An optional boolean, whether or not to
                    discard old bond data before adding in bonds determined by
                    distance. True by default.
        :param bool delete_excessive_bonds: An optional boolean, whether or not
                    to check for and delete excessive bonds. True by default.
        """

        if not numpy.class_dependency("calculate bonds by distance", "NUMPY"):
            return

        if not numpy.class_dependency("calculate bonds by distance", "SCIPY"):
            return

        atom_inf = self.__parent_molecule.information.get_atom_information()
        consts = self.__parent_molecule.get_constants()

        # create/recreate the bond array if needed
        if (remove_old_bond_data == True or
            self.__parent_molecule.get_bonds() is None):

            self.__parent_molecule.set_bonds(
                numpy.zeros((len(atom_inf), len(atom_inf)))
            )

        # get the longest bond length on record
        max_bond_length = numpy.max([
            consts['bond_length_dict'][key]
            for key in consts['bond_length_dict'].keys()
        ])

        # which ones could possibly be bound (less than the max_bond_length)
        distances = numpy.squareform(numpy.pdist(self.__parent_molecule.get_coordinates()))
        ones_to_consider = numpy.nonzero(distances < max_bond_length * 1.2)

        for index in range(len(ones_to_consider[0])):
            index1 = ones_to_consider[0][index]
            index2 = ones_to_consider[1][index]

            atom_inf = self.__parent_molecule.get_atom_information()
            consts = self.__parent_molecule.get_constants()

            if index1 != index2:
                # so an atom is not bound to itself.__parent_molecule
                key = (atom_inf['element'][index1] + '-' +
                    atom_inf['element'][index2])

                try: bond_dist = consts['bond_length_dict'][key]
                except:

                    print("ERROR: Unknown bond distance between elements " +
                        atom_inf['element'][index1] + ' and ' +
                        atom_inf['element'][index2] +
                        '. Assuming ' + str(max_bond_length) + '.')
                    bond_dist = max_bond_length

                if (distances[index1][index2] < bond_dist * 1.2 and
                    distances[index1][index2] > bond_dist * 0.5):

                    # so they should be bonded
                    self.__parent_molecule.add_bond(index1, index2)

        if delete_excessive_bonds == True:
            # now do a sanity check. C cannot have more than 4 bonds, O cannot
            # have more than 2, and N cannot have more than 2 if more, than use
            # ones closest to ideal bond length
            parnt = self.__parent_molecule
            atom_inf = parnt.get_atom_information()
            sel_func = parnt.select_all_atoms_bound_to_selection
            consts = self.__parent_molecule.get_constants()
            for index in range(len(atom_inf)):
                # get the info of the index atom
                element = atom_inf['element'][index]

                bond_partner_indices = sel_func(numpy.array([index]))
                number_of_bonds = len(bond_partner_indices)

                try:
                    if (number_of_bonds >
                        consts['max_number_of_bonds_permitted'][element]):

                        # so this atom has too many bonds
                        # get the distances of this atoms bonds
                        dists = distances[index][bond_partner_indices]

                        # get the ideal distances of those bonds
                        # initialize the vector
                        ideal_dists = numpy.empty(len(dists))

                        for t in range(len(bond_partner_indices)):
                            # populate the ideal-bond-length vector
                            index_partner = bond_partner_indices[t]
                            element_partner = atom_inf['element'][
                                index_partner]
                            ideal_dists[t] = consts['bond_length_dict'][
                                element + '-' + element_partner]

                        # get the distance
                        diff = numpy.absolute(dists - ideal_dists)

                        # identify the bonds to discard
                        indices_in_order = diff.argsort()
                        indx_2_thrw_out = indices_in_order[
                            consts['max_number_of_bonds_permitted'][element]:
                        ]
                        indx_2_thrw_out = bond_partner_indices[indx_2_thrw_out]

                        # discard the extra bonds
                        for throw_out_index in indx_2_thrw_out:
                            self.__parent_molecule.delete_bond(index,
                                                               throw_out_index)

                except:
                    pass # element probably wasn't in the dictionary

    def get_number_of_bond_partners_of_element(self, atom_index, the_element):
        """
        Counts the number of atoms of a given element bonded to a specified
        atom of interest.

        Requires the :any:`numpy` library.

        Should be called via the wrapper function
        :meth:`~scoria.Molecule.Molecule.get_number_of_bond_partners_of_element`.

        :param int atom_index: An int, the index of the atom of interest.
        :param str the_element: A string describing the element of the neighbors
                    to be counted.

        :returns: An int, the number of neighboring atoms of the specified
                element.
        :rtype: *int*
        """

        if not numpy.class_dependency("count the number of bond partners of a given element", "NUMPY"):
            return

        # this function is really here for historical reasons. it's similar to
        # the old number_of_neighors_of_element function. it could be done
        # pretty easily with numpy

        sel_func = self.__parent_molecule.select_all_atoms_bound_to_selection
        atom_inf = self.__parent_molecule.get_atom_information()

        the_element = the_element.strip()
        bond_partners_selection = sel_func(numpy.array([atom_index]))
        elements = atom_inf['element'][bond_partners_selection]
        return len(numpy.nonzero(elements == the_element)[0])

    def get_index_of_first_bond_partner_of_element(self, atom_index,
                                                   the_element):
        """
        For a given atom of interest, returns the index of the first
        neighbor of a specified element.

        Requires the :mod:`numpy` and :mod:`scipy.spatial` libraries.

        Should be called via the wrapper function
        :meth:`~scoria.Molecule.Molecule.get_index_of_first_bond_partner_of_element`.

        :param int atom_index: An int, the index of the atom of interest.
        :param str the_element: A string specifying the desired element of the
                    neighbor.

        :returns: An int, the index of the first neighbor atom of the specified
                element. If no such neighbor exists, returns -1.
        :rtype: *int*
        """

        if not numpy.class_dependency("get the index of a bond partner of element X", "NUMPY"):
            return

        if not numpy.class_dependency("get the index of a bond partner of element X", "SCIPY"):
            return

        # this function is really here for historical reasons. it's similar to
        # the old index_of_neighbor_of_element function. it could be done
        # pretty easily with numpy

        sel_func = self.__parent_molecule.select_all_atoms_bound_to_selection
        atom_inf = self.__parent_molecule.get_atom_information()

        the_element = the_element.strip()
        bond_partners_selection = sel_func(numpy.array([atom_index]))
        elements = atom_inf['element'][bond_partners_selection]

        return bond_partners_selection[numpy.nonzero(elements == the_element)[0]][0]

    def delete_bond(self, index1, index2):
        """
        Deletes a bond.

        Should be called via the wrapper function :meth:`~scoria.Molecule.Molecule.delete_bond`.
        
        :param int index1: An int, the index of the first atom of the bonded
                    pair.
        :param int index2: An int, the index of the second atom of the bonded
                    pair.
        """

        if not numpy.class_dependency("delete bonds", "NUMPY"):
            return

        bonds = self.__parent_molecule.get_bonds()
        try:
            bonds[index1][index2] = 0
            bonds[index2][index1] = 0
        except: print(("Could not delete bond between " + str(index1) +
                       " and " + str(index2) + "."))

    def add_bond(self, index1, index2, order = 1):
        """
        Adds a bond.

        Should be called via the wrapper function :meth:`~scoria.Molecule.Molecule.add_bond`.

        :param int index1: An int, the index of the first atom of the bonded
                    pair.
        :param int index2: An int, the index of the second atom of the bonded
                    pair.
        :param int order: An optional int, the order of the bond. 1 by default.
        """

        if not numpy.class_dependency("add bonds", "NUMPY"):
            return

        bonds = self.__parent_molecule.get_bonds()
        bonds[index1][index2] = order
        bonds[index2][index1] = order
        self.__parent_molecule.set_bonds(bonds)

    def delete_atom(self, index):
        """
        Deletes an atom.

        Should be called via the wrapper function :meth:`~scoria.Molecule.Molecule.delete_atom`.

        :param int index: An int, the index of the atom to delete.
        """

        if not numpy.class_dependency("delete atoms", "NUMPY"):
            return

        # remove the atom information
        self.__parent_molecule.set_atom_information(
            numpy.delete(
                self.__parent_molecule.get_atom_information(),
                index
            )
        )

        # remove the coordinates from all timesteps
        for frame in range(0, self.__parent_molecule.get_trajectory_frame_count() ):
            self.__parent_molecule.set_coordinates(
                numpy.delete(
                    self.__parent_molecule.get_coordinates(frame),
                    index, axis = 0
                ), frame
            )

        try:
            self.__parent_molecule.set_coordinates_undo_point(
                numpy.delete(
                    self.__parent_molecule.get_coordinates_undo_point(),
                    index, axis = 0
                )
            )
        except: pass

        # remove the relevant bonds
        self.__parent_molecule.set_bonds(
            numpy.delete(
                self.__parent_molecule.information.get_bonds(),
                index, 0
            )
        )

        self.__parent_molecule.set_bonds(
            numpy.delete(
                self.__parent_molecule.information.get_bonds(),
                index, 1
            )
        )

        # the hierarchy will have to be recomputed
        self.__hierarchy = {}

    def add_atom(self, record_name = "ATOM", serial = 1, name = "X",
                 resname = "XXX", chainid = "X", resseq = 1, occupancy = 0.0,
                 tempfactor = 0.0, charge = '', element = "X",
                 coordinates = numpy.array([0.0, 0.0, 0.0]), autoindex = True):
        """
        Adds an atom.

        Should be called via the wrapper function :meth:`~scoria.Molecule.Molecule.add_atom`.

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

        if not numpy.class_dependency("add atoms", "NUMPY"):
            return

        # add the atom information

        if len(record_name) < 6:
            record_name = record_name.ljust(6)

        if len(name) < 5:
            if len(name) < 4:
                name = name.rjust(4) + ' '
            else:
                name = name.rjust(5)

        if len(resname) < 4:
            resname = resname.rjust(4)

        if len(chainid) < 2:
            chainid = chainid.rjust(2)

        if len(charge) < 2:
            charge = charge.ljust(2)

        if len(element) < 2:
            element = element.rjust(2)

        name_stripped = name.strip()
        resname_stripped = resname.strip()
        chainid_stripped = chainid.strip()
        element_stripped = element.strip()

        consts = self.__parent_molecule.get_constants()

        try:
            mass = consts['mass_dict'][element_stripped]
        except:
            mass = 0.0

        # if there is no atom_information, you need to create it.
        if self.__parent_molecule.get_atom_information() is None:
            self.__parent_molecule.information.set_atom_information(
                numpy.zeros(
                    (1,),
                    dtype = [('record_name', '|S6'), ('serial', '<i8'),
                             ('name_padded', '|S5'), ('resname_padded', '|S4'),
                             ('chainid_padded', '|S2'), ('resseq', '<i8'),
                             ('occupancy', '<f8'), ('tempfactor', '<f8'),
                             ('element_padded', '|S2'), ('charge', '|S2'),
                             ('name', '|S5'),
                             ('resname', '|S4'),
                             ('chainid', '|S2'),
                             ('element', '|S2')]
                )
            )

        # ********

        atom_information = numpy.ma.resize(
            self.__parent_molecule.get_atom_information(),
            self.__parent_molecule.information.get_total_number_of_atoms() + 1
        )

        atom_information['record_name'][-1] = record_name
        atom_information['name_padded'][-1] = name
        atom_information['resname_padded'][-1] = resname
        atom_information['chainid_padded'][-1] = chainid
        atom_information['charge'][-1] = charge
        atom_information['element_padded'][-1] = element
        atom_information['name'][-1] = name_stripped
        atom_information['resname'][-1] = resname_stripped
        atom_information['chainid'][-1] = chainid_stripped
        atom_information['element'][-1] = element_stripped
        atom_information['serial'][-1] = serial
        atom_information['resseq'][-1] = resseq
        atom_information['occupancy'][-1] = occupancy
        atom_information['tempfactor'][-1] = tempfactor

        self.__parent_molecule.set_atom_information(atom_information)
        #self.__parent_molecule.assign_masses()

        if 'mass' in atom_information.dtype.names:
            atom_information['mass'][-1] = mass

        # now add the coordinates
        if self.__parent_molecule.get_coordinates() is None:
            self.__parent_molecule.set_coordinates(numpy.array([coordinates]))
        else:
            for frame in range(0, len(self.__parent_molecule.get_trajectory_coordinates())):
                self.__parent_molecule.set_coordinates(
                    numpy.vstack((
                        self.__parent_molecule.get_coordinates(frame), coordinates
                    )),
                    frame
                )

        tot_atoms = self.__parent_molecule.get_total_number_of_atoms()

        # now add places for bonds, though bonds will only be added if done
        # explicitly, not here
        if self.__parent_molecule.get_bonds() is None:
            self.__parent_molecule.set_bonds(numpy.array([[0]]))
        else:
            self.__parent_molecule.set_bonds(numpy.vstack((
                self.__parent_molecule.information.get_bonds(),
                numpy.zeros(tot_atoms - 1)
            )))

            self.__parent_molecule.set_bonds(numpy.hstack((
                self.__parent_molecule.get_bonds(),
                numpy.zeros((1, tot_atoms)).T
            )))
