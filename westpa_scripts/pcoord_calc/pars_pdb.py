# Tools for Processing PDB Files.
# 
# PDB Format Functions (below) via: 
# https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
# Mostly followed that, mostly.
#
# ATOM
#1-4    "ATOM"                                  character   "{0:<4}  "
#7-11   Atom Serial Number              right   integer     "{0:>5} "
#13-16  Atom Name                       left*   character   " {0:<3}"
#17     Alternate Location Indicator            character   "{0:<1}"
#18-20  Residue Name                    right   character   "{0:>3} "
#22     Chain Identifier                        character   "{0:<1}"
#23-26  Residue Sequence Number         right   integer     "{0:>4}"
#27     Insertion Codes of Residues             character   "{0:<1}   "
#31-38  X Orthogonal A Coordinate       right   real        "{0:>8}"
#39-46  Y Orthogonal A Coordinate       right   real        "{0:>8}"
#47-54  Z Orthogonal A Coordinate       right   real        "{0:>8}"
#55-60  Occupancy                       right   real        "{0:>6}"
#61-66  Temperature Factor              right   real        "{0:>6}      "
#73-76  Segment Identifier              left    character   "{0:>4}"
#77-78  Element Symbol                  right   character   "{0:>2}"
#
# TER
#1-3    "TER"                                   character   "{0:<3}   "
#7-11   Serial Number                   right   integer     "{0:>5}      "
#18-20  Residue Name                    right   character   "{0:>3} "
#22     Chain Identifier                        character   "{0:<1}"
#23-26  Residue Sequence Number         right   integer     "{0:>4}"
#27     Insertion Codes of Residues             character   "{0:<1}"
#
# ATOM
def line_to_pdb_atom(line):
    fields = []
    fields.append(line[0:7].strip())        # 0  ATOM
    fields.append(line[7:11].strip())       # 1  Atom Serial Number
    fields.append(line[12:16].strip())      # 2  Atom Name
    fields.append(line[16].strip())         # 3  Alternate Location Indicator
    fields.append(line[17:20].strip())      # 4  Residue Name
    fields.append(line[21].strip())         # 5  Chain Identifier
    fields.append(line[22:26].strip())      # 6  Residue Sequence Number
    fields.append(line[26].strip())         # 7  Insertion Codes of Residues
    fields.append(line[30:38].strip())      # 8  X Orthogonal A Coordinate
    fields.append(line[38:46].strip())      # 9  Y Orthogonal A Coordinate
    fields.append(line[46:54].strip())      # 10 Z Orthogonal A Coordinate
    fields.append(line[54:60].strip())      # 11 Occupancy
    fields.append(line[60:66].strip())      # 12 Temperature Factor
    fields.append(line[72:76].strip())      # 13 Segment Identifier
    fields.append(line[76:78].strip())      # 14 Element Symbol
    return fields
def pdb_atom_to_line(fields):
    line = (
         "{0:<6.6}".format(fields[0])         # 0  ATOM
        +"{0:>5.5} ".format(fields[1])        # 1  Atom Serial Number
        +" {0:<3.3}".format(fields[2])        # 2  Atom Name
        +"{0:<1.1}".format(fields[3])         # 3  Alternate Location Indicator
        +"{0:>3.3} ".format(fields[4])        # 4  Residue Name
        +"{0:<1.1}".format(fields[5])         # 5  Chain Identifier
        +"{0:>4.4}".format(fields[6])         # 6  Residue Sequence Number
        +"{0:<1.1}   ".format(fields[7])      # 7  Insertion Codes of Residues
        +"{0:>8.7}".format(fields[8])         # 8  X Orthogonal A Coordinate
        +"{0:>8.7}".format(fields[9])         # 9  Y Orthogonal A Coordinate
        +"{0:>8.7}".format(fields[10])        # 10 Z Orthogonal A Coordinate
        +"{0:>6.6}".format(fields[11])        # 11 Occupancy
        +"{0:>6.6}      ".format(fields[12])  # 12 Temperature Factor
        +"{0:>4.4}".format(fields[13])        # 13 Segment Identifier
        +"{0:>2.2}".format(fields[14])        # 14 Element Symbol
        )
    return line
# TER
def line_to_pdb_ter(line):
    fields = []
    fields.append(line[0:3].strip())        # 0  TER
    fields.append(line[6:11].strip())       # 1  Atom Serial Number
    fields.append(line[17:20].strip())      # 2  Residue Name
    fields.append(line[21].strip())         # 3  Chain Identifier
    fields.append(line[22:26].strip())      # 4  Residue Sequence Number
    fields.append(line[26].strip())         # 5  Insertion Codes of Residues
    return fields
def pdb_ter_to_line(fields):
    line = (
         "{0:<3.3}   ".format(fields[0])      # 0  ATOM
        +"{0:>5.5}      ".format(fields[1])   # 1  Atom Serial Number
        +"{0:>3.3} ".format(fields[2])        # 2  Residue Name
        +"{0:<1.1}".format(fields[3])         # 3  Chain Identifier
        +"{0:>4.4}".format(fields[4])         # 4  Residue Sequence Number
        +"{0:<1.1}   ".format(fields[5])      # 5  Insertion Codes of Residues
        )
    return line