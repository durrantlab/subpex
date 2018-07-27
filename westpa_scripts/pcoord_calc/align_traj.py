#!/usr/bin/env python

### 
### This script is specifically for the last equilibrium trajectory.
### 
### What this script does:
###   - Aligns Trajectories
###   - Strips Waters
###   - Converts to pdb
### 
### How to use Script:
###   - open terminal 
###   - python align_traj.py mol.psf eq.dcd ref.pdb needs to have a temp/ folder in same directory
###   - [Enter]
### 


# Imports
import MDAnalysis
from MDAnalysis.analysis import align
import os
import sys


# Functions                                                                                                 # Functions
def _align_trajectory_(traj, sel_str="name CA", pdb_filename=None):                                             # Align Trajectory Function
    
    # Get Reference                                                                                             # Get Reference
    ref = MDAnalysis.Universe(pdb_filename)                                                                     # Load Reference Structure
    sel_ref = ref.select_atoms(sel_str)                                                                         # Select Alignment Basis

    # Align Trajectory                                                                                          # Align Trajectory
    alignment = align.AlignTraj(traj, ref, in_memory=True, select=sel_str)                                      # Set Alignment Settings
    alignment.run()                                                                                             # Align
    
    # Return Aligned Trajectory                                                                                 # Return Aligned Trajectory
    return traj                                                                                                 # Return Aligned Trajectory


# Main                                                                                                      # Main
if __name__ == "__main__":
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
        PDB.write(eq_sel.atoms)                                                                             # Write Atoms
