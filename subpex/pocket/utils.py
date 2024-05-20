from typing import Sequence


def link_resnames(resnames: Sequence[str], linker: str = "or") -> str:
    """Convert a list of residue names into a formatted selection string.

    Args:
        resnames: List of residue names to string together.
        linker: String to link (e.g., put between) residue names.

    Returns:
        str: Formatted selection string.

    Example:
        >>> resnames = ["WAT", "HOH", "H2O"]
        >>> formatted_string = format_resnames(resnames)
        >>> print(formatted_string)
        resname WAT or resname HOH or resname H2O
        >>> formatted_string = format_resnames(resnames, linker="and")
        >>> print(formatted_string)
        resname WAT and resname HOH and resname H2O
    """
    return f" {linker} ".join(f"resname {resname}" for resname in resnames)
