#!/usr/bin/python

### 
### Note: SubPEX Settings are Hard Coded here.
### 
### What this script does:
###   - Replaces Chain in SubPEX Settings
### 
### How to use Script:
###   - open terminal 
###   - python fix_settings.py mol.pdb
###   - [Enter]
### 


# Imports for Main
import pars_pdb
if (__name__ == "__main__"):
    import sys
    import os


# Functions                                                         # Functions
def fix_settings(pdb_path):
    
    # Initialize Files                                              # Initialize Files
    pdb_old    = open(pdb_path)                                     # Get  old PDB File
    
    for line in pdb_old:                                            # For Each Line
        fields = line.strip().split()                               # Split to Columns
        
        if(fields[0] == "ATOM"):                                    # At Model Start, Delete Model (Don't Write)
            chain = pars_pdb.line_to_pdb_atom(line)[5]              # Get Chain ID
            
            # Update Chain                                          # Update Chain in SubPEX Settings
            settings            = "SubPEX_settings_pre"                 # Default Settings Location
            subpex_settings     = open(settings)                        # Default Settings Open
            new_settings        = "SubPEX_settings"                     # New Settings Location
            subpex_new_settings = open(new_settings,"w")                # New Settings Open
            for line in subpex_settings:                                # Splits Rows
                fields = line.strip().split()                           # Splits Columns
                if(fields[0].lower() == "chain"):                       # Finds Chain
                    subpex_new_settings.write("chain           "+chain) # Rename Chain
                else:                                                   # Anything Else, Leave Alone
                    subpex_new_settings.write(line)                     # Leave Alone (Write to settings)
            break                                                   # Stop Looking at the PDB
        else:                                                       # Anything Else, Leave Alone
            continue                                                # Delete Model (Don't Write)
    pdb_old.close()                                                 # Close Old PDB


# Main                                                          # Main
if (__name__ == "__main__"):
    pdb_path = os.path.abspath(sys.argv[1])                         # Get PDB
    fix_settings(pdb_path)                                          # Convert to Trajectory