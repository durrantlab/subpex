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
from scoria.Quaternion import Quaternion


class OtherMolecules():
    """
    A class for characterizing the relationships between multiple
    scoria.Molecule objects.
    """

    def __init__(self, parent_molecule_object):
        """
        Initializes the scoria.OtherMolecules class.

        :param scoria.Molecule parent_molecule_object: The scoria.Molecule object
                associated with this class.
        """

        self.__parent_molecule = parent_molecule_object

    def get_other_molecules_aligned_to_this(self, other_mol, tethers,
                                           weight_mat = None):
        """
        Aligns a molecule to self (this scoria.Molecule object) using a
        quaternion RMSD alignment.

        Requires the :any:`numpy` library.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_other_molecules_aligned_to_this`
                
        :param scoria.Molecule other_mol: A scoria.Molecule that is to be aligned to
                    this one.
        :param list[list] tethers: A list of lists, where the inner list is
                          (tether1_index, tether2_index). That inner list can
                          also be a numpy array. Each inner array contains the
                          indices of self and other_mol, respectively, such
                          that equivalent atoms are listed in the same order.
                          So, for example, if (atom 1, self = atom 3, other)
                          and (atom2,  self = atom6, other) than the tethers
                          would be (numpy.array([1, 2]), numpy.array([3, 6])),
                          or [(1, 2), (3, 6)].

        :returns: The new molecule.
        """

        if not numpy.class_dependency("align molecules. Missing the dot-product function", "NUMPY"):
            return
        
        tethers = numpy.array(tethers).T

        # Adapted from Itzhack Y. Bar-Itzhack. New Method for Extracting the
        # Quaternion from a Rotation Matrix. Journal of Guidance, Control, and
        # Dynamics 2000
        if tethers is None:
            raise Exception('No tethers specified for RMSD alignment')
        elif tethers.shape[0] != 2:
            raise Exception('Tethers should have only 2 columns')  # Because see T above, which makes rows at this point.

        # If weight_matrix isn't specified, then treat all atoms equally
        if weight_mat is None: 
            weight_mat = numpy.identity(tethers.shape[1])
        
        # get the atoms corresponding to the tethers, in tether order
        self_static_atom_coordinates = (
            self.__parent_molecule.get_coordinates()[tethers[0]]
        )

        other_dynamic_atom_coordinates = (
            other_mol.get_coordinates()[tethers[1]]
        )

        # translate the tether atoms to the origin
        center_self = numpy.mean(self_static_atom_coordinates, 0)
        center_other = numpy.mean(other_dynamic_atom_coordinates, 0)

        self_static_atom_coordinates = (self_static_atom_coordinates -
                                        center_self)

        other_dynamic_atom_coordinates = (other_dynamic_atom_coordinates -
                                          center_other)

        # get optimal rotation
        M = numpy.dot(
            numpy.dot(
                numpy.transpose(self_static_atom_coordinates),
                weight_mat
            ),
            other_dynamic_atom_coordinates
        )

        #Create symmetric 4x4 matrix K from M
        K = numpy.array([[M[0, 0] + M[1, 1] + M[2, 2],
                          M[1, 2] - M[2, 1],
                          M[2, 0] - M[0, 2],
                          M[0, 1] - M[1, 0]
                         ],
                         [
                          M[1, 2] - M[2, 1],
                          M[0, 0] - M[1, 1] - M[2, 2],
                          M[1, 0] + M[0, 1],
                          M[2, 0] + M[0, 2]
                         ],
                         [
                          M[2, 0] - M[0, 2],
                          M[1, 0] + M[0, 1],
                          M[1, 1] - M[0, 0] - M[2, 2],
                          M[1, 2] + M[2, 1]
                         ],
                         [
                          M[0, 1] - M[1, 0],
                          M[2, 0] + M[0, 2],
                          M[1, 2] + M[2, 1],
                          M[2, 2] - M [0, 0] - M[1, 1]
                         ]
                        ]
        )

        # Find eigenvector associated with the most positive eigenvalue of K.
        # Multiple quaternions can
        E, V = numpy.linalg.eig(K)
        index = numpy.argmax(E)
        eigenvector = V[:, index]
        rot_quat = Quaternion(eigenvector[0], eigenvector[1],
                              eigenvector[2], eigenvector[3])

        rot_mat = rot_quat.to_matrix()

        #Apply translation and rotation to the other molecule
        new_mol = other_mol.copy()

        new_mol.set_coordinates(new_mol.information.get_coordinates() -
                                center_other)

        new_mol.set_coordinates(numpy.dot(
            new_mol.information.get_coordinates(),
            rot_mat
        ))

        new_mol.set_coordinates(new_mol.information.get_coordinates() +
                                center_self)

        return new_mol

    def steric_clash_with_another_molecule(self, other_mol, cutoff,
                                           pairwise_comparison = True):
        """
        Detects steric clashes between the scoria.Molecule (self) and
        another scoria.Molecule.

        Requires the :any:`numpy` and :any:`scipy<scipy.spatial>` libraries.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.steric_clash_with_another_molecule`

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

        if not numpy.class_dependency("calculate the steric clashes with another molecule", "NUMPY"):
            return

        if not numpy.class_dependency("calculate the steric clashes with another molecule", "SCIPY"):
            return

        prnt = self.__parent_molecule

        sel_cls_aatms_frm_diff_mols = (
            prnt.select_close_atoms_from_different_molecules
        )

        if pairwise_comparison == True:
            # so use a simple pairwise comparison to find close atoms
            indices1, indices2 = sel_cls_aatms_frm_diff_mols(other_mol,
                                                             cutoff, True)
        else: # so the more sophisticated heirarchical method
            # terminate early is true because you don't want all close ones
            indices1, indices2 = sel_cls_aatms_frm_diff_mols(other_mol, cutoff,
                                                             False, True)

        if len(indices1) == 0 and len(indices2) == 0:
            return False
        else:
            return True

    def merge_with_another_molecule(self, other_molecules):
        """
        Merges two molecular models into a single model.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.merge_with_another_molecule`
        
        :param scoria.Molecule other_molecules: A molecular model (scoria.Molecule
                    object).

        :returns: A single scoria.Molecule object containing the atoms of
                    this model combined with the atoms of other_molecules.
        """

        merged = self.__parent_molecule.copy()

        # if masses have been assigned to either molecule, they must be
        # assigned to both
        slf_atm_inf = self.__parent_molecule.get_atom_information()
        if ('mass' in merged.information.get_atom_information().dtype.names or
            'mass' in slf_atm_inf.dtype.names):
            self.__parent_molecule.assign_masses()
            merged.information.assign_masses()

        merged.filename = ""
        merged.get_remarks().extend(other_molecules.get_remarks())

        merged.set_atom_information(
            numpy.stack_arrays(
                (
                    merged.get_atom_information(),
                    other_molecules.get_atom_information()
                ),
                usemask = False
            )
        )

        merged.set_coordinates(numpy.vstack((
            merged.get_coordinates(), other_molecules.get_coordinates()
        )))

        merged.set_coordinates_undo_point(None)

        # merge the bonds, though note that bonds between the two molecules
        # will not be set
        if (not merged.get_bonds() is None and
            not other_molecules.get_bonds() is None):

            bonds1 = merged.get_bonds().copy()
            bonds2 = other_molecules.get_bonds().copy()

            bonds1_v2 = numpy.hstack((
                bonds1, numpy.zeros((len(bonds1), len(bonds2)))
            ))

            bonds2_v2 = numpy.hstack((
                numpy.zeros((len(bonds2), len(bonds1))), bonds2
            ))

            merged.set_bonds(numpy.vstack((bonds1_v2, bonds2_v2)))
        else:
            merged.set_bonds(None)

        # the molecule center will be redefined, so you might as well start the
        # hierarchy all over
        try:
            del merged.information.hierarchy['spheres']
        except:
            pass

        return merged

    def get_distance_to_another_molecule(self, other_molecules,
                                         pairwise_comparison = True):
        """
        Computes the minimum distance between any of the atoms of this
        molecular model and any of the atoms of a second specified model.

        Requires the :any:`numpy` and :any:`scipy<scipy.spatial>` libraries.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_distance_to_another_molecule`

        :param scoria.Molecule other_molecules: a scoria.Molecule, the other molecular
                    model.
        :param bool pairwise_comparison: An optional boolean, whether or not to
                    perform a simple pairwise distance comparison (if True) or
                    to use a more sophisitcated method (if False). True by
                    default.

        :returns: A float, the minimum distance between any two atoms of the two
                specified molecular models (self and other_molecules).
        """

        if not numpy.class_dependency("calculate the distance to another molecule", "NUMPY"):
            return

        if not numpy.class_dependency("calculate the distance to another molecule", "SCIPY"):
            return

        if pairwise_comparison == True:
            return numpy.amin(numpy.cdist(
                self.__parent_molecule.get_coordinates(),
                other_molecules.get_coordinates()
            ))
        else:
            # so use the more sophisticated methods for comparison
            # note that this is not the fastest way to do this, but it uses
            # existing functions
            # and is still pretty fast, so I'm going to stick with it.

            # first, get a cutoff distance. Let's just do a quick survey of the
            # two molecules to pick a good one.
            slf_gt_crs = self.__parent_molecule.get_coordinates()
            oth_gt_crs = other_molecules.get_coordinates()
            self_tmp = slf_gt_crs[numpy.arange(0, len(slf_gt_crs),
                                               len(slf_gt_crs) / 10.0,
                                               dtype = int)]

            other_tmp = oth_gt_crs[numpy.arange(0, len(oth_gt_crs),
                                                len(oth_gt_crs) / 10.0,
                                                dtype = int)]

            cutoff = numpy.amin(numpy.cdist(self_tmp, other_tmp))

            # now get all the indices that come within that cutoff
            prnt = self.__parent_molecule

            sel_cls_atms_dif_mols = (
                prnt.select_close_atoms_from_different_molecules
            )

            self_indices, other_indices = sel_cls_atms_dif_mols(other_molecules,
                                                                cutoff, False)

            self_coors = self.__parent_molecule.get_coordinates()[self_indices]
            self_other = other_molecules.get_coordinates()[other_indices]

            return numpy.amin(numpy.cdist(self_coors, self_other))

    def get_rmsd_equivalent_atoms_specified(self, other_mol, tethers):
        """
        Calculates the RMSD between this scoria.Molecle object and
        another, where equivalent atoms are explicitly specified.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_rmsd_equivalent_atoms_specified`
        
        :param scoria.Molecule other_mol: The other scoria.Molecule object.
        :param tuple tethers: A tuple of two numpy.array objects, where each array
                    contains the indices of self and other_mol, respectively,
                    such that equivalent atoms are listed in the same order.
                    So, for example, if (atom 1, self = atom 3, other) and
                    (atom2, self = atom6, other) than the tethers would be
                    (numpy.array([1, 2]), numpy.array([3, 6])).

        :returns: A float, the RMSD between self and other_mol.
        """

        tethers = numpy.transpose(tethers)

        slf_gt_crs = self.__parent_molecule.get_coordinates()
        if (len(slf_gt_crs) !=
            len(other_mol.get_coordinates())):

            print("Cannot calculate RMSD: number of atoms are not equal.")
            print(("\t" + (str(len(slf_gt_crs)) +
                          " vs. " + str(len(other_mol.get_coordinates())) +
                          " atoms.")))
            return 99999999.0

        self_coors_in_order = slf_gt_crs[tethers[0]]
        other_coors_in_order = other_mol.get_coordinates()[tethers[1]]

        delta = self_coors_in_order - other_coors_in_order
        norm_squared = numpy.sum(delta * delta, axis = -1)
        rmsd = numpy.power(numpy.sum(norm_squared) / len(norm_squared), 0.5)
        return rmsd

    def get_rmsd_order_dependent(self, other_mol):
        """
        Calculates the RMSD between two structures, where equivalent atoms
        are listed in the same order.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_rmsd_order_dependent`
        
        :param scoria.Molecule other_mol: The other scoria.Molecule object.

        :returns: A float, the RMSD between self and other_mol.
        """

        self_index_in_order = numpy.arange(
            0, len(self.__parent_molecule.get_coordinates()), 1, dtype = int
        )

        other_index_in_order = numpy.arange(
            0, len(other_mol.get_coordinates()), 1, dtype = int
        )

        return self.get_rmsd_equivalent_atoms_specified(
            other_mol,
            list(zip(self_index_in_order, other_index_in_order))
        )

    def get_rmsd_heuristic(self, other_mol):
        """
        Caluclates the RMSD between two identical molecules with different
        conformations, per the definition given in "AutoDock Vina: Improving
        the speed and accuracy of docking with a new scoring function,
        efficient optimization, and multithreading,"" by Oleg Trott and Arthur
        J. Olson. Note: Identical means the order of the atoms is the same as
        well.

        Requires the :any:`numpy` and :any:`scipy<scipy.spatial>` libraries.

        Wrapper function for :meth:`~scoria.Molecule.Molecule.get_rmsd_heuristic`

        :param scoria.Molecule other_mol: The other scoria.Molecule object.

        :returns: A float, the RMSD between self and other_mol.
        """

        if not numpy.class_dependency("calculate an RMSD using a heuristical algorithm", "NUMPY"):
            return

        # Group the other_mol atoms by element (atom type in pdbqt speak)
        atom_inf = self.__parent_molecule.get_atom_information()
        self_atom_coors = self.__parent_molecule.get_coordinates()
        other_atom_coors = other_mol.get_coordinates()

        self_atom_grps = {}
        other_atom_grps = {}

        for i, atm in enumerate(atom_inf):
            element_stripped = atm["element"]
            if not element_stripped in self_atom_grps.keys():
                self_atom_grps[element_stripped] = []
                other_atom_grps[element_stripped] = []
            self_atom_grps[element_stripped].append(self_atom_coors[i])
            other_atom_grps[element_stripped].append(other_atom_coors[i])
        
        for element in self_atom_grps.keys():
            self_atom_grps[element] = numpy.array(self_atom_grps[element])
            other_atom_grps[element] = numpy.array(other_atom_grps[element])
        
        # Calculate the rmsds going both ways
        rmsd1 = self._get_rmsd_heuristic_helper_func(self_atom_grps, other_atom_grps)
        rmsd2 = self._get_rmsd_heuristic_helper_func(other_atom_grps, self_atom_grps)
        
        return numpy.max((rmsd1, rmsd2))

    def _get_rmsd_heuristic_helper_func(self, atom_grp1, atom_grp2):
        """
        A helper function for calculating heuristic RMSD.

        Requires the :any:`scipy<scipy.spatial>` library.

        :param dict atom_grp1: A dictionary, where the keys are atom types
                    and the values are numpy arrays of the coordinates.
        :param dict atom_grp2: The same, but now given the atoms of the
                    other molecule.
            
        :returns: A float, the heuristic RMSD between the two molecules
                    (atom_grp1 to atom_grp2)
        """
        
        if not numpy.class_dependency("use this helper function", "SCIPY"):
            return

        # Go through each of the element types
        tot_sum_min_dists_sqrd = 0
        num_heavy_atms = 0
        for element in atom_grp1.keys():
            if element[:1] != "H":  # Because only heavy atoms
                coor1 = atom_grp1[element]
                coor2 = atom_grp2[element]
                
                dists = numpy.cdist(coor1, coor2)
                dists_sqr = dists * dists
                min_dists_sqr = numpy.min(dists_sqr, axis = 0)

                tot_sum_min_dists_sqrd = tot_sum_min_dists_sqrd + numpy.sum(min_dists_sqr) 
                num_heavy_atms = num_heavy_atms + len(coor1)
        rmsd = numpy.sqrt(tot_sum_min_dists_sqrd / num_heavy_atms)
        return rmsd

