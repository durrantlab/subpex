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

class Manipulation():
    """
    A class for translating and rotating the atomic coordinates of a
    scoria.Molecule object
    """

    def __init__(self, parent_molecule_object):
        """
        Initializes the scoria.Manipulation class.

        :param scoria.Molecule parent_molecule_object: The scoria.Molecule object
                    associated with this class.

            """

        self.__parent_molecule = parent_molecule_object

    def set_coordinate_undo_point(self):
        """
        Sets ("saves") the undo point of the atom coordinates. Any
        subsequent manipulations of atomic coordinates can be "undone" by
        reseting to this configuration via the coordinate_undo function.
        
        """

        self.__parent_molecule.set_coordinates_undo_point(
            copy.deepcopy(self.__parent_molecule.get_trajectory_coordinates())
        )

    def coordinate_undo(self):
        """
        Resets the coordinates of all atoms to those saved using the
        set_coordinate_undo_point function.
        """

        self.__parent_molecule.set_trajectory_coordinates(
            copy.deepcopy(self.__parent_molecule.get_coordinates_undo_point())
        )

    def set_atom_location(self, atom_index, new_location):
        """
        Translates the entire molecular model (without rotating) so that the
        atom with the specified index is located at the specified coordinate.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.set_atom_location`
        
        :param int atom_index: An int, the index of the target atom.
        :param numpy.array new_location: A numpy.array specifying the new (x, y, z)
                    coordinate of the specified atom.

        :returns: A numpy.array specifying the (delta_x, delta_y, delta_z) vector
                by which the pmolecule.Molecule was translated.
        """

        if new_location.shape == (3,):
            new_location = numpy.array([new_location])

        currentloc = self.__parent_molecule.get_coordinates()[atom_index]
        delta = new_location - currentloc

        self.translate_molecule(delta)

        return delta

    def translate_molecule(self, delta):
        """
        Translate all the atoms of the molecular model by a specified
        vector.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.translate_molecule`

        :param numpy.array delta: A numpy.array (delta_x, delta_y, delta_z) specifying the
            amount to move each atom along the x, y, and z coordinates.
        """

        if delta.shape == (3,):
            delta = numpy.array([delta])

        self.__parent_molecule.set_coordinates(
            self.__parent_molecule.get_coordinates() + delta
        )

        hierarchy = self.__parent_molecule.get_hierarchy()
        if hierarchy is not None and 'spheres' in list(hierarchy.keys()):
            gt_hrchy_sph = self.__parent_molecule.get_hierarchy()['spheres']

            # so update location of hierarchical elements
            gt_hrchy_sph['molecule']['center'] = (
                gt_hrchy_sph['molecule']['center'] + delta
            )

            gt_hrchy_sph['chains']['centers'] = (
                gt_hrchy_sph['chains']['centers'] + delta
            )

            gt_hrchy_sph['residues']['centers'] = (
                gt_hrchy_sph['residues']['centers'] + delta
            )

    def rotate_molecule_around_a_line_between_points(self, line_point1,
                                                     line_point2, rotate):
        """
        Rotate the molecular model about a line segment. The end points of
        the line segment are explicitly specified coordinates.

        Wrapper function for 
        :meth:`~scoria.Molecule.Molecule.rotate_molecule_around_a_line_between_points`

        :param numpy.array line_point1: A numpy.array (x, y, z) corresponding to one end
                    of the line segment.
        :param numpy.array line_point2: A numpy.array (x, y, z) corresponding to the
                    other end of the line segment.
        :param float rotate: A float, the angle of rotation, in radians.
        """

        if line_point1.shape == (1, 3):
            line_point1 = line_point1[0]

        if line_point2.shape == (1, 3):
            line_point2 = line_point2[0]

        a = line_point1[0]
        b = line_point1[1]
        c = line_point1[2]
        #d = line_point2[0]
        #e = line_point2[1]
        #f = line_point2[2]

        delta = line_point2 - line_point1

        u = delta[0] #d - a
        v = delta[1] #e - b
        w = delta[2] #f - c

        v_2_plus_w_2 = numpy.power(v, 2) + numpy.power(w, 2)
        u_2_plus_w_2 = numpy.power(u, 2) + numpy.power(w, 2)
        u_2_plus_v_2 = numpy.power(u, 2) + numpy.power(v, 2)
        u_2_plus_v_2_plus_w_2 = u_2_plus_v_2 + numpy.power(w, 2)

        cos = numpy.cos(rotate)
        sin = numpy.sin(rotate)

        coordinates = self.__parent_molecule.get_coordinates()

        ux_plus_vy_plus_wz = numpy.sum(
            coordinates * delta, 1
        )

        # Now rotate molecule. In a perfect world, I'd have an awesome better
        # numpified version of this, perhaps with tensor or matrix
        # multiplication

        for t in range(len(coordinates)):
            # so t is an atom index
            x_not, y_not, z_not = coordinates[t]
            ux_plus_vy_plus_wz = u * x_not + v * y_not + w * z_not

            coordinates[t][0] = (
                a * v_2_plus_w_2 + u * (-b * v - c * w + ux_plus_vy_plus_wz) +
                (
                    -a * v_2_plus_w_2 +
                    u * (b * v + c * w - v * y_not - w * z_not) +
                    v_2_plus_w_2 * x_not
                ) * cos +
                numpy.sqrt(u_2_plus_v_2_plus_w_2) *
                (-c * v + b * w - w * y_not + v * z_not) * sin
            ) # /u_2_plus_v_2_plus_w_2

            coordinates[t][1] = (
                b * u_2_plus_w_2 + v * (-a * u - c * w + ux_plus_vy_plus_wz) +
                (
                    -b * u_2_plus_w_2 +
                    v * (a * u + c * w - u * x_not - w * z_not) +
                    u_2_plus_w_2 * y_not
                ) * cos +
                numpy.sqrt(u_2_plus_v_2_plus_w_2) *
                (c * u - a * w + w * x_not - u * z_not) * sin
            ) # /u_2_plus_v_2_plus_w_2

            coordinates[t][2] = (
                c * u_2_plus_v_2 + w * (-a * u - b * v + ux_plus_vy_plus_wz) +
                (
                    -c * u_2_plus_v_2 +
                    w * (a * u + b * v - u * x_not - v * y_not) +
                    u_2_plus_v_2 * z_not
                ) * cos +
                numpy.sqrt(u_2_plus_v_2_plus_w_2) *
                (-b * u + a * v - v * x_not + u * y_not) * sin
            ) # /u_2_plus_v_2_plus_w_2

        coordinates = coordinates * (1.0 / u_2_plus_v_2_plus_w_2)

        self.__parent_molecule.set_coordinates(
            coordinates
        )

        # here I'm going to just delete the hierarchical info because I'm lazy.
        try:
            # calculated bounding spheres, if any, are no longer valid.
            del self.__parent_molecule.get_hierarchy()['spheres']
        except:
            pass

    def rotate_molecule_around_a_line_between_atoms(self, line_point1_index,
                                                    line_point2_index, rotate):
        """
        Rotate the molecular model about a line segment. The end points of
        the line segment are atoms of specified indices.

        Wrapper function for 
        :meth:`~scoria.Molecule.Molecule.rotate_molecule_around_a_line_between_atoms`

        :param int line_point1_index: An int, the index of the first atom at one
                    end of the line segment.
        :param int line_point2_index: An int, the index of the second atom at
                    the other end of the line segment.
        :param float rotate: A float, the angle of rotation, in radians.
        """

        pt1 = self.__parent_molecule.get_coordinates()[line_point1_index]
        pt2 = self.__parent_molecule.get_coordinates()[line_point2_index]
        self.rotate_molecule_around_a_line_between_points(pt1, pt2, rotate)

        try:
            # calculated bounding spheres, if any, are no longer valid.
            del self.__parent_molecule.get_hierarchy()['spheres']
        except:
            pass

    def rotate_molecule_around_pivot_point(self, pivot, thetax,
                                           thetay, thetaz):
        """
        Rotate the molecular model around a specified atom.

        Requires the :any:`numpy` library.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.rotate_molecule_around_pivot_point`

        :param numpy.array pivot: A numpy.array, the (x, y, z) coordinate about which
                    the molecular model will be rotated.
        :param float thetax: A float, the angle to rotate relative to the x axis,
                    in radians.
        :param float thetay: A float, the angle to rotate relative to the y axis,
                    in radians.
        :param float thetaz: A float, the angle to rotate relative to the z axis,
                    in radians.
        """

        if not numpy.class_dependency("rotate the molecule about a point. Missing the dot-product function", "NUMPY"):
            return

        if pivot.shape == (3,): pivot = numpy.array([pivot])

        # First, move the Molecule so the pivot is at the origin
        self.__parent_molecule.set_coordinates(
            self.__parent_molecule.get_coordinates() - pivot
        )

        # do the rotation
        sinx = numpy.sin(thetax)
        siny = numpy.sin(thetay)
        sinz = numpy.sin(thetaz)
        cosx = numpy.cos(thetax)
        cosy = numpy.cos(thetay)
        cosz = numpy.cos(thetaz)

        rot_matrix = numpy.array([
            [(cosy * cosz),
             (sinx * siny * cosz + cosx * sinz),
             (sinx * sinz - cosx * siny * cosz)
            ],
            [-(cosy * sinz),
             (cosx * cosz - sinx * siny * sinz),
             (cosx * siny * sinz + sinx * cosz)
            ],
            [siny, -(sinx * cosy), (cosx * cosy)]
        ])

        self.__parent_molecule.set_coordinates(
            numpy.dot(rot_matrix, self.__parent_molecule.get_coordinates().T).T
        )

        # now move the pivot point back to it's old location
        self.__parent_molecule.set_coordinates(
            self.__parent_molecule.get_coordinates() + pivot
        )

        try:
            # calculated bounding spheres, if any, are no longer valid.
            del self.__parent_molecule.information.hierarchy['spheres']
        except:
            pass

    def rotate_molecule_around_pivot_atom(self, pivot_index, thetax,
                                          thetay, thetaz):
        """
        Rotate the molecular model around a specified atom.

        Requires the :any:`numpy` library.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.rotate_molecule_around_pivot_atom`
        
        :param int pivot_index: An int, the index of the atom about which the
                    molecular model will be rotated.
        :param float thetax: A float, the angle to rotate relative to the x axis,
                    in radians.
        :param float thetay: A float, the angle to rotate relative to the y axis,
                    in radians.
        :param float thetaz: A float, the angle to rotate relative to the z axis,
                    in radians.
        """

        if not numpy.class_dependency("rotate the molecule about an atom. Missing the dot-product function", "NUMPY"):
            return

        pivot = self.__parent_molecule.get_coordinates()[pivot_index]
        self.rotate_molecule_around_pivot_point(pivot, thetax, thetay, thetaz)

        try:
            # calculated bounding spheres, if any, are no longer valid.
            del self.__parent_molecule.get_hierarchy()['spheres']
        except:
            pass
