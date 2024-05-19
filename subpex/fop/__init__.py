"""
The fop module provides functions for reading, writing, and parsing fields of points (FOP)
in XYZ and PDB file formats. These functions support both single-frame and multi-frame FOP data,
making it easy to visualize molecular structures and fields of points using software such as VMD.

Functions:
    read_fop: Determines the file type (.xyz or .pdb) and parses the FOP accordingly.
    write_fop: Writes the FOP to a file in either XYZ or PDB format, handling single-frame and multi-frame data.
"""
