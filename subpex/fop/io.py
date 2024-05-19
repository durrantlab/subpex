"""Reading and writing field of points."""

import os
from collections.abc import Sequence

from loguru import logger

from ..configs import SubpexConfig
from .prop import get_fop_volume


def parse_pdb_fop(file_path: str) -> Sequence[Sequence[float]]:
    """Parse and return the field of points that is stored in a pdb file.

    Args:
        file_path: Path to the pdb file with the field of points.

    Returns:
        Field of points XYZ coordinates.
    """
    logger.debug(f"Parsing {file_path} for FOP")
    with open(file_path, "r", encoding="utf-8") as f:
        text = f.readlines()
    fop = []
    for line in text:
        if line.startswith("ATOM"):
            split = line.split()
            x = float(split[4])
            y = float(split[5])
            z = float(split[6])
            fop.append([x, y, z])
        else:
            pass
    return fop


def parse_xyz_fop(file_path: str) -> Sequence[Sequence[float]]:
    """Reads a filename and parses the field of points to convert it into a
    list of lists.

    Args:
        file_path: filename of the field of points to open and parse

    Returns:
        Field of points XYZ coordinates.
    """
    logger.debug(f"Parsing {file_path} for FOP")
    # open reference fop xyz file
    with open(file_path, "r", encoding="utf-8") as f:
        text = f.readlines()

    # parse xyz file
    fop = []
    for i in text[2:]:
        line = i.split()
        if len(line) == 4:
            point = [float(line[1]), float(line[2]), float(line[3])]
            fop.append(point)
        else:
            pass

    return fop


def read_fop(file_path: str) -> Sequence[Sequence[float]]:
    """Load field of points from supported file types.

    Args:
        file_path: Path to xyz or pdb file.

    Returns:
        Loaded field of points.
    """
    fop_type = file_path.split(".")[-1]
    if fop_type.lower() == "xyz":
        fop_ref = parse_xyz_fop(file_path)
    elif fop_type.lower() == "pdb":
        fop_ref = parse_pdb_fop(file_path)
    else:
        raise ValueError("ref_fop_write must end in `.xyz` or `.pdb`.")
    return fop_ref


def points_to_pdb(file_path: str, fop: Sequence[Sequence[float]]) -> None:
    """Writes a pdb file full of C-alphas so users can visualize the field of
    points in software such as VMD.

    Args:
        filename: name of the pdb file to create.
        fop: Field of points XYZ coordinates.
    """
    logger.debug(f"Writing FOP to {file_path}")
    # TODO: Erich - modify to be able to do multiframe pdb
    with open(file_path, "w", encoding="utf-8") as f:
        # f.write(header)
        atom_number = 1
        for i in fop:
            text = "ATOM" + "{:>7}".format(atom_number)
            text += "{:^6}".format("CA")
            text += "{:>3}".format("ALA")
            text += "{:>6}".format("A")
            text += "{:>12}".format(i[0])
            text += "{:>8}".format(i[1])
            text += "{:>8}".format(i[2])
            text += "{:>6}".format("1.0")
            text += "{:>6}".format("1.0")
            text += "\n"
            f.write(text)
            atom_number += 1
        f.write("\n")


def points_to_xyz(
    file_path: str, fop: Sequence[Sequence[float]], resolution: float, radius: float
) -> None:
    """Takes the coordinates and the resolution, writes write an xyz file that
    users can load into visualization software such as VMD.

    Args:
        file_path: name for the xyz file to be created.
        fop: contains XYZ coordinates for each atom.
        resolution: resolution in Angstroms.
        radius: radius in Angstrom for the field of points.
    """
    logger.debug(f"Writing FOP to {file_path}")

    volume = get_fop_volume(fop, resolution)

    # TODO: Erich - modify to be able to do multiframe xyz
    with open(file_path, "w", encoding="utf-8") as f:
        f.write(str(len(fop)) + "\n")
        f.write(
            "Generated in SubPEx. Res: {res}    Rad: {dim} Ang    Volume: {vol} Ang^3 \n".format(
                res=resolution, dim=radius, vol=volume
            )
        )
        # Adding each coordinate as a C-alpha for visualization purposes
        for i in fop:
            f.write(f"CA    {i[0]}    {i[1]}    {i[2]} \n")


def write_fop(
    fop: Sequence[Sequence[float]],
    f_path: str,
    spx_config: SubpexConfig,
    data_dir: str | None = None,
) -> None:
    if data_dir is None:
        data_dir = ""
    f_path = os.path.join(data_dir, f_path)

    fop_type = f_path.split(".")[-1]
    if fop_type.lower() == "xyz":
        points_to_xyz(
            f_path,
            fop,
            spx_config.pocket.resolution,
            spx_config.pocket.radius,
        )
    elif fop_type.lower() == "pdb":
        points_to_pdb(f_path, fop)
    else:
        raise ValueError("f_path must end in `.xyz` or `.pdb`.")
