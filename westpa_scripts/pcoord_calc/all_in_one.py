# all in one

import MDAnalysis
from MDAnalysis.analysis import align
import os
import sys
from align_traj import _align_trajectory_
import pars_pdb
from fix_settings import fix_settings
import multiprocessing
import numpy as np
import scipy as sp
import scoria
import SubPEX_tweaked

# Load Files                                                                                                # Load Files
psf             = os.path.abspath(sys.argv[1])                                                              # System Parameters
dcd_eq          = os.path.abspath(sys.argv[2])                                                              # Equilibrium Trajectory 
pdb_ref         = os.path.abspath(sys.argv[3])                                                              # Reference Structure

# Name Output                                                                                               # Name Output 
pdb_out         = "temp/"+sys.argv[2][:-4]+'_aligned.pdb'                                                   # Add _aligned suffix

# Get Trajectories                                                                                          # Get Trajectories  
traj_eq         = MDAnalysis.Universe(psf, dcd_eq)                                                          # Equilibrium Trajectory 

# Align Trajectories                                                                                        # Align Trajectories
traj_eq_aligned = _align_trajectory_(traj_eq,pdb_filename=pdb_ref)                                          # Equilibrium Trajectory 


# Save as PDB                                                                                               # Save as PDB
with MDAnalysis.Writer(pdb_out, multiframe=True, bonds=None, n_atoms=traj_eq_aligned.atoms.n_atoms) as PDB: # Writer "PDB"
    
    # Write Equilibrium Trajectory                                                                          # Write Equilibrium Trajectory
    eq_sel = traj_eq_aligned.select_atoms('protein')                                                        # Strip Waters
    traj_eq_aligned.trajectory[-1]                                                                          # Last Frame Only
    PDB.write(eq_sel.atoms)


pdb_path = os.path.abspath(sys.argv[4])                         # Get PDB
fix_settings(pdb_path)                                                                                # Write Atoms

# get inputs from user
#inputs = get_inputs()
settings = get_settings()
inputs = settings[0:3]


# initialize scoria Molecule by reading in the file from the inputs
mol = scoria.Molecule(inputs[-1])

# get the binding pocket center 
center_unsnapped = np.array(settings[3:6], float)
center = get_snapped_center(center_unsnapped, inputs[0])

# create a static field_of_points to be referenced but never modified
reference_fop = field_of_points()

# get all the alpha carbon indices
all_alphas = get_all_alphas()

# get the multi-frame coordinates of all the alpha carbons using aforementioned
# indices
all_alpha_trajectory = all_alphas.get_trajectory_coordinates()
print 'There are {} frames, each containing {} alpha carbons'.format(len(all_alpha_trajectory),len(all_alpha_trajectory[0]))

# only keep the alpha carbons coordinates in each frame that are in range of
# the binding pocket dimension parameter obtained from the user
kept_alphas = keep_pocket_alphas()

# summary outputs
for frame in kept_alphas:
    print "This frame has {} alpha carbons in the binding pocket".format(len(frame))

# get the coordinates for only those atoms whose coordinates in each frame are
# in range of the binding pocket dimension parameter obtained from the user
kept_atoms = keep_pocket_atoms()

# summary outputs
for frame in kept_atoms:
    print "This frame has {} atoms in the binding pocket".format(len(frame))

# now use the data to do the pocket analysis
for frame_index in xrange(len(kept_atoms)):                         # for each frame present grab
    frame_atoms = kept_atoms[frame_index]                           # coords for all atoms/alphas
    frame_alphas = kept_alphas[frame_index]                         # kept for this frame
    local_fop = reference_fop                                       # local dynamic copy of fop
    tree_atoms = sp.spatial.cKDTree(frame_atoms)                    # distances btw atoms
    tree_fop = sp.spatial.cKDTree(local_fop)                        # distances in field of points
    r = tree_atoms.query_ball_tree(tree_fop, 1.5, p=2, eps=0)       # points too close to atoms
    eliminate_points_that_clash_w_protein(local_fop)                # eliminate too close points
    hull = sp.spatial.ConvexHull(frame_alphas)                      # create reference hull

    def hull_checker_inputs():                                      # function to set up for hull checker
        hull_inputs = []                                        # because NAN points don't need to be
        for index in xrange(len(local_fop)):                    # checked... we already know they aren't
            if not np.isnan(index):                         # in the hull
                hull_inputs.append(index)
            else:
                pass
        return hull_inputs

    to_be_checked = hull_checker_inputs()                           # run the setup function

    def hull_checker(index):                                        # function to check whether a point
        point = reference_fop[index]                            # is within the convex hull made by
        if not np.isnan(point[0]):                              # kept alphas.  Essentially if adding
            points = list(frame_alphas)                     # that point to alphas changes the 
            points.append(point)                            # hull's vertices (see scipy.spatial
            points = np.array(points)                       # ConvexHull documentation), then it's 
            new_hull = sp.spatial.ConvexHull(points)        # outside the hull
            if list(hull.vertices) == list(new_hull.vertices):
                return index
            else:
                return None

        def chunk_data_for_multiprocess():               
                shape = []                                              # list to store indices of points in hull;
                if len(to_be_checked) <= 32766:
                        shape = multithreading(to_be_checked, 4, hull_checker)
                else:
                        chunk_indices = []
                        start = 0
                        end = 32767
                        while len(to_be_checked) > end - 1:
                                chunk_indices.append([start,end])
                                start += 32767
                                end += 32767
                        chunk_indices.append([start, len(to_be_checked)])
                        for i in chunk_indices:
                                in_hull = multithreading(a[i[0]:i[1]], 4, hull_checker)
                                for j in in_hull:
                                        shape.append(j)
                return shape

    pocket_shape = chunk_data_for_multiprocess()                    # if something other than None gets
    final_shape = []
    for index in pocket_shape:
            if index != None:
                    final_shape.append(reference_fop[index])
            else:
                    pass

    # summary outputs
    print len(final_shape)                                        
    print "Pocket volume: {}".format(len(final_shape)*(inputs[0]**3))

    # file output
    out_path = inputs[2][:-4] + '_frame{}'.format(frame_index) + '.xyz'
    g = open(out_path, 'w')
    g.write(str(len(final_shape)) + "\n")
    g.write("Generated in SubPEx.\tRes: {}\tDim: {} Ang\tVolume: {} Ang^3\n".format(inputs[0],inputs[1],len(final_shape)*(inputs[0]**3)))
    for point in final_shape:
        line = "CA\t{}\t{}\t{}\n".format(point[0],point[1],point[2])
        g.write(line)
    # summary output
    print "Done with frame {} of {}".format(frame_index+1, len(kept_atoms))