"""Reading and writing field of points."""

import os
from collections.abc import Sequence

from loguru import logger


def parse_pdb_fop(file_path: str) -> Sequence[Sequence[float]]:
    """Parse and return the field of points (FOP) stored in a PDB file.

    This function reads a PDB file, extracts the XYZ coordinates from lines that start
    with 'ATOM', and returns them as a list of lists of floats.

    Args:
        file_path: Path to the PDB file containing the field of points.

    Returns:
        A list of [X, Y, Z] coordinates.

    Raises:
        FileNotFoundError: If the file at the given path does not exist.
        ValueError: If the file contains improperly formatted lines.
    """
    logger.debug(f"Parsing {file_path} for FOP")
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
    except FileNotFoundError as e:
        logger.error(f"File not found: {file_path}")
        raise e

    fop = []
    for line in lines:
        if line.startswith("ATOM"):
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                fop.append([x, y, z])
            except ValueError as e:
                logger.error(f"Error parsing line: {line.strip()}")
                raise e

    logger.debug(f"Successfully parsed {len(fop)} points from {file_path}")
    return fop


def parse_xyz_fop(file_path: str) -> Sequence[Sequence[float]]:
    """Reads a filename and parses the field of points (FOP) from an XYZ file format.

    This function reads an XYZ file, extracts the XYZ coordinates from the lines, and
    returns them as a list of lists of floats. The first two lines are typically headers
    and are skipped.

    Args:
        file_path: Filename of the field of points to open and parse.

    Returns:
        A list of [X, Y, Z] coordinates.

    Raises:
        FileNotFoundError: If the file at the given path does not exist.
        ValueError: If the file contains improperly formatted lines.
    """
    logger.debug(f"Parsing {file_path} for FOP")
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
    except FileNotFoundError as e:
        logger.error(f"File not found: {file_path}")
        raise e

    try:
        fop = [
            [float(line.split()[1]), float(line.split()[2]), float(line.split()[3])]
            for line in lines[2:]
            if len(line.split()) == 4
        ]
    except ValueError as e:
        logger.error(f"Error parsing file: {file_path}")
        raise e

    logger.debug(f"Successfully parsed {len(fop)} points from {file_path}")
    return fop


def read_fop(file_path: str) -> Sequence[Sequence[float]]:
    """Load field of points (FOP) from supported file types (.xyz or .pdb).

    This function determines the file type based on the file extension and
    calls the appropriate parsing function to load the field of points.

    Args:
        file_path: Path to the .xyz or .pdb file containing the field of points.

    Returns:
        A list of [X, Y, Z] coordinates.

    Raises:
        ValueError: If the file extension is not supported.
        FileNotFoundError: If the file at the given path does not exist.
        ValueError: If the file contains improperly formatted lines.
    """
    fop_type = file_path.split(".")[-1].lower()
    logger.debug(f"Reading FOP from {file_path}, detected file type: {fop_type}")

    if fop_type == "xyz":
        fop_ref = parse_xyz_fop(file_path)
    elif fop_type == "pdb":
        fop_ref = parse_pdb_fop(file_path)
    else:
        logger.error("Unsupported file extension: must be .xyz or .pdb")
        raise ValueError("File must end in `.xyz` or `.pdb`.")

    logger.debug(f"Successfully loaded FOP from {file_path}")
    return fop_ref


def points_to_pdb(
    file_path: str,
    fop_input: Sequence[Sequence[float]] | Sequence[Sequence[Sequence[float]]],
) -> None:
    """Writes a PDB file full of C-alphas so users can visualize the field of
    points in software such as VMD.

    This function creates a PDB file with each point represented as a C-alpha
    (CA) atom, allowing for visualization in molecular visualization software.
    The function supports creating a multiframe PDB file if the input sequence
    contains multiple frames of points.

    Args:
        file_path: Name of the PDB file to create.
        fop_input: Field of points XYZ coordinates. Can contain multiple frames of points.
    """
    logger.debug(f"Writing FOP to {file_path}")

    # Ensure fop is a 3D list
    if all(isinstance(i, (float, int)) for i in fop_input[0]):
        fop = [fop_input]
    else:
        fop = fop_input

    with open(file_path, "w", encoding="utf-8") as f:
        frame_number = 1
        for frame in fop:
            atom_number = 1
            f.write(f"MODEL     {frame_number}\n")
            for point in frame:
                text = f"ATOM  {atom_number:>5}  CA  ALA A {atom_number:>4}    "
                text += f"{point[0]:>8.3f}{point[1]:>8.3f}{point[2]:>8.3f}  1.00  0.00           C\n"
                f.write(text)
                atom_number += 1
            f.write("ENDMDL\n")
            frame_number += 1

    logger.debug(f"Successfully wrote FOP to {file_path}")


def points_to_xyz(
    file_path: str,
    fop_input: Sequence[Sequence[float]] | Sequence[Sequence[Sequence[float]]],
) -> None:
    """Writes an XYZ file that users can load into visualization software such as VMD.

    This function creates an XYZ file with each point represented as a CA atom,
    allowing for visualization in molecular visualization software. The function
    supports creating a multiframe XYZ file if the input sequence contains multiple
    frames of points.

    Args:
        file_path: Name for the XYZ file to be created.
        fop: Contains XYZ coordinates for each atom. Can be a 2D list (single frame)
             or a 3D list (multiple frames).
    """
    logger.debug(f"Writing FOP to {file_path}")

    # Ensure fop is a 3D list
    if all(isinstance(i, (float, int)) for i in fop_input[0]):
        fop = [fop_input]  # Wrap in an outer list to make it a 3D list
    else:
        fop = fop_input

    with open(file_path, "w", encoding="utf-8") as f:
        for frame in fop:
            f.write(f"{len(frame)}\n\n")  # Number of atoms and a blank line
            for point in frame:
                f.write(f"CA    {point[0]:.3f}    {point[1]:.3f}    {point[2]:.3f}\n")

    logger.debug(f"Successfully wrote FOP to {file_path}")


def write_fop(
    fop: Sequence[Sequence[float]] | Sequence[Sequence[Sequence[float]]],
    f_path: str,
    data_dir: str | None = None,
) -> None:
    """Writes the field of points (FOP) to a file in either XYZ or PDB format.

    This function saves the given field of points to a specified file path. It supports
    both XYZ and PDB file formats and can handle both single-frame and multi-frame data.

    Args:
        fop: Field of points XYZ coordinates. Can be a 2D list (single frame) or a 3D list (multiple frames).
        f_path: File path where the FOP should be saved. The file extension must be either `.xyz` or `.pdb`.
        data_dir: Directory to save the file in. If None, the file will be saved in the current directory.

    Raises:
        ValueError: If the file extension is not `.xyz` or `.pdb`.
    """
    if data_dir is None:
        data_dir = ""

    full_path = os.path.join(data_dir, f_path)
    fop_type = full_path.split(".")[-1].lower()

    if fop_type == "xyz":
        points_to_xyz(full_path, fop)
    elif fop_type == "pdb":
        points_to_pdb(full_path, fop)
    else:
        raise ValueError("f_path must end in `.xyz` or `.pdb`.")
