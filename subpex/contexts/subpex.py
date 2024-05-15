from typing import Any

from collections.abc import Iterable, MutableSequence

from loguru import logger

from .base import BaseContextManager, BaseContextValidator


class SubpexContextManager(BaseContextManager):
    """Contexts for SubPEx."""

    # pylint: disable-next=too-many-statements
    def __init__(
        self, yaml_paths: str | Iterable[str] | None = None, **kwargs: dict[str, Any]
    ) -> None:
        """
        Args:
            yaml_paths: Path(s) to YAML file(s) to load into the context.
        """
        self.restart: bool = False
        """Restart a subpex simulation if directly already exists.

        TODO: run the sp_restart script if westpa file exists and this is True.
        """
        self.md_engine: str = "Amber"
        """Molecular dynamics engine to use. Both `Amber` and `NAMD` are implemented.
        """
        self.pocket_center: Iterable[float] | None = None
        """The center of the pocket."""
        self.pocket_radius: float = 10.9
        """Pocket radius to consider."""
        self.progress_coord: MutableSequence[str] = ["composite"]
        """Progress coordinates to calculate.

        -   `jd` for Jaccard distance,
        -   `prmsd` for pocket heavy atoms RMSD,
        -   `bb` for backbone RMSD,
        -   `composite` for composite RMSD.
        """
        self.clustering_engine: str = "cpptraj"
        """Clustering engine to use. `cpptraj` is the only option at the moment.
        """
        self.calculated_points: int = -1
        """Number of point to calculate per trajectory segment. If `-1`, it will
        calculate all.
        """
        self.n_clusters: int = 25
        """Number of independent clusters to identify."""
        self.min_n_clusters_gen_bin: int = 3
        """TODO:"""
        self.max_n_clusters_gen_bin: int = 25
        """TODO:"""
        self.clustering_region: str = "pocket"
        """Region to separate into clusters. This can be `pocket` or `backbone`."""
        self.clustering_method: str = "hierarchical"
        """TODO:"""
        self.cluster_generation: bool = True
        """TODO:"""
        self.cluster_bin: bool = False
        """TODO:"""
        self.plot_n_walkers: bool = False
        """TODO:"""
        self.fop_filetype: str = "xyz"
        """TODO:"""

        super().__init__(yaml_paths, **kwargs)


# pylint: disable-next=too-few-public-methods
class SubpexContextValidator(BaseContextValidator):
    """Base class for validating Subpex contexts."""

    # pylint: disable=unused-argument

    @staticmethod
    def write(value: Any, context: dict[str, Any]) -> bool:
        """Validate `write`"""
        return True

    @staticmethod
    def pocket_center(value: Any, context: dict[str, Any]) -> bool:
        """Validate `pocket_center`"""
        if len(value) != 3:
            logger.error("pocket_center must be a list of x, y, and z.")
            logger.error(f"pocket_center length should be 3, but is {len(value)}")
            return False
        for ele in value:
            if not isinstance(ele, float):
                logger.error(
                    f"pocket_center values must be `float` type; found {type(ele)}"
                )
                return False
        return True

    @staticmethod
    def pocket_radius(value: Any, context: dict[str, Any]) -> bool:
        """Validate `pocket_radius`"""
        if not isinstance(value, float):
            logger.error(f"pocket_radius value must be `float`; found {type(ele)}")
            return False
        return True


def check_input_settings_file(filename: str) -> dict:
    """Loads the file containing the SubPEx settings and checks that all the
    parameters are available. Otherwise, sets defaults if possible.

    Args:
        filename (str): filename of the settings yaml file (e.g., west.cfg).

    Raises:
        IOError: raises error if files do not exist or can't be read.

    Returns:
        settings (dict): settings for running SubPEx and some of the analysis tools.
    """

    # Checking that we specify which progress coordinates we will print in the pcoord file
    available_pcoords = ["jd", "prmsd", "bb", "composite", "pvol", "rog"]
    if "pcoord" in settings and len(settings["pcoord"]) < 1 or "pcoord" not in settings:
        logging.critical(
            "There is an error with the progress coordinate (pcoord) setting. Need to have one of "
            + (", ".join(available_pcoords))
        )
        sys.exit("Error with setting pcoord")
    else:
        for i in settings["pcoord"]:
            if i not in available_pcoords:
                logging.critical(
                    "There is an error with the progress coordinate (pcoord) setting. Need to have one of "
                    + (", ".join(available_pcoords))
                )
                sys.exit(f"Error with setting pcoord: {i} is not a valid pcoord")

    # Cheking resolution setting
    if "resolution" in settings and type(settings["resolution"]) == float:
        pass
    else:
        settings["resolution"] = 0.5
        logging.warning("Using a default resolution value of 0.5 Angstrom.")

    # Cheking calculated_points setting
    if "calculated_points" in settings and type(settings["calculated_points"]) == int:
        pass
    else:
        settings["calculated_points"] = -1
        logging.warning(
            "Using a default calculated_points value of -1. This will calculate progress coordinates for all the frames in the trajectory file"
        )

    # Checking west_home setting
    if "west_home" not in settings or len(glob.glob(settings["west_home"])) != 1:
        logging.critical("There is an error with the west home directory")
        sys.exit("There is an error with the west home directory")
    elif len(glob.glob(settings["west_home"] + "/traj_segs")) != 1:
        logging.warning(
            "Could not locate the traj_segs directory. Please be sure it exists and the trajectories are there before running the clustering script. Note: This should not exist if you have not run SubPEx"
        )
    else:
        pass
    if settings["west_home"][-1] != "/":
        settings["west_home"] = settings["west_home"] + "/"

    # checking reference and topology files
    check_file_exists(settings, "reference")
    check_file_exists(settings, "topology")

    # checking that the selection_file exists, it is different for get_reference_fop.py script
    if __file__ == "jdistance.py" or __file__ == "clustering.py":
        check_file_exists(settings, "selection_file")
    else:
        if "selection_file" not in settings:
            settings["selection_file"] = settings["west_home"] + "/selection_string.txt"
            logging.warning(
                "Selection file was not determined used the default of {}/selection_string.txt".format(
                    settings["west_home"]
                )
            )

    # making sure the reference_fop exists and the format is in the settings file
    accepted_fop_files = ["xyz", "pdb"]
    if (
        "fop_filetype" not in settings
        and settings["fop_filetype"] not in accepted_fop_files
    ):
        settings["fop_filetype"] = "xyz"

    if __file__ == "jdistance.py":
        check_file_exists(settings, "reference_fop")
    else:
        if "reference_fop" not in settings:
            settings["reference_fop"] = settings[
                "west_home"
            ] + "/reference_fop.{}".format(settings["fop_filetype"])
            logging.warning(
                "Selection file was not determined used the default of {}/reference_fop.{}".format(
                    settings["west_home"]
                ),
                settings["fop_filetype"],
            )

    logging.info("The parameters used for {} are {}".format(__file__, settings))
    return settings
