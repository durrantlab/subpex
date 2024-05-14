from typing import Any

from collections.abc import Iterable, MutableMapping

from loguru import logger
from ruamel.yaml import YAML


class WestpaContextManager:
    """Contexts for Westpa.

    More information on the configuration file can be found
    [here](https://github.com/westpa/westpa/wiki/Configuration-File).
    """

    # pylint: disable-next=too-many-statements
    def __init__(
        self, yaml_paths: str | Iterable[str] | None = None, **kwargs: dict[str, Any]
    ) -> None:
        """
        Args:
            yaml_paths: Path(s) to YAML file(s) to load into the context.
        """
        self.system_driver: str = "adaptive.System"
        """The `driver` parameter must be set to a subclass of `WESTSystem`, and given
        in the form `module.class`.
        """
        self.system_module_path: str = "$WEST_SIM_ROOT/adaptive_binning/"
        """The `module_path` parameter is appended to the system path and indicates
        where the class is defined.
        """
        self.propagation_gen_istates: bool = True
        """Boolean specifying whether to generate initial states from the basis states.
        The executable propagator defines a specific configuration block
        (add internal link to other section), and custom propagators should override
        the `WESTPropagator.gen_istate()` method.
        """
        self.propagation_block_size: int = 1
        """An integer defining how many segments should be passed to a worker at a time.
        When using the serial work manager, this value should be set to the maximum
        number of segments per iteration to avoid significant overhead incurred by the
        locking mechanism in the WMFutures framework. Parallel work managers might
        benefit from setting this value greater than one in some instances to
        decrease network communication load.
        """
        self.propagation_save_transition_matrices: bool = True
        """No description was provided.
        """
        self.propagation_max_run_wallclock: str | None = None
        """A time in `dd:hh:mm:ss` or `hh:mm:ss` specifying the maximum wallclock time
        of a particular WESTPA run. If running on a batch queuing system, this time
        should be set to less than the job allocation time to ensure that WESTPA shuts
        down cleanly.
        """
        self.propagation_max_total_iterations: int | None = None
        """An integer value specifying the number of iterations to run. This parameter is checked against the last completed iteration stored in the HDF5 file, not the number of iterations completed for a specific run. The default value of None only stops upon external termination of the code.
        """
        self.data_west_data_file: str = "west.h5"
        """The name of the main HDF5 data storage file for the WESTPA simulation.
        """
        self.data_aux_compression_threshold: int = 1048576
        """The threshold in bytes for compressing the auxiliary data in a dataset on an iteration-by-iteration basis.
        """
        self.data_iter_prec: int = 8
        """The length of the iteration index with zero-padding. For the default value, iteration 1 would be specified as iter_00000001.
        """
        self.data_load_auxdata: bool = False
        """All datasets without an explicit load: specification are loaded.
        """
        self.data_datasets: Iterable[MutableMapping[str, Any]] | None = None
        """A list of all datasets to load in. Each item must have the following
        information.

        - `h5path`: Custom path to dataset within the H5 file.
        - `store`: Store auxdata into H5 File.
        - `dtype`: Data type of auxdata to be saved.
        - `scaleoffset`: HDF5 dataset scale offset filter. Defaults to None.
        - `compression`: HDF5 dataset compression. Defaults to None.
        - `chunks`: HDF5 dataset chunk size. Defaults to None.
        """
        self.data_data_refs_segment: str = (
            r"$WEST_SIM_ROOT/traj_segs/{segment.n_iter:06d}/{segment.seg_id:06d}"
        )
        """TODO:"""
        self.data_data_refs_basis_state: str = (
            r"$WEST_SIM_ROOT/bstates/{basis_state.auxref}"
        )
        """TODO:"""
        self.data_data_refs_initial_state: str = (
            r"$WEST_SIM_ROOT/istates/{initial_state.iter_created}/{initial_state.state_id}.xml"
        )
        """TODO:"""

        if isinstance(yaml_paths, str):
            yaml_paths = [yaml_paths]
        if yaml_paths is not None:
            for yaml_path in yaml_paths:
                self.from_yaml(yaml_path)
        self.update(kwargs)

    def from_yaml(self, yaml_path: str | None) -> None:
        """Load context information from a YAML file. This will only update data
        contained in the YAML file.

        Args:
            yaml_path: Path to YAML file to load.
        """
        if yaml_path is not None:
            logger.info("Loading YAML context from {}", yaml_path)
            yaml = YAML(typ="safe")
            with open(yaml_path, "r", encoding="utf-8") as f:
                yaml_data = yaml.load(f)
            logger.debug("YAML data:\n{}", yaml_data)
            self.update(yaml_data)
        self.yaml_path = yaml_path

    def update(self, attr_dict: MutableMapping[str, Any]) -> dict[str, Any]:
        """Update attributes with values from the provided dictionary.

        Args:
            attr_dict: Dictionary containing attribute names and their
            corresponding values.
        """
        logger.debug("Updating context:\n{}", attr_dict)
        for key, value in attr_dict.items():
            setattr(self, key, value)
        return self.get()

    def get(self) -> dict[str, Any]:
        """Retrieve the context.

        Returns:
            A dictionary representing the current context.
        """
        # The following line filters methods and attributes like __dict__.
        context = {
            k: v for k, v in vars(self).items() if not callable(v) and "__" not in k
        }
        logger.debug("Retrieved context:\n{}", context)
        return context

    def __enter__(self) -> dict[str, Any]:
        """Enter the context and return the current context as a dictionary."""
        return self.get()

    def __exit__(self, exc_type, exc_value, exc_tb):
        """Exit the context.

        Args:
            exc_type: Type of the exception.
            exc_value: Value of the exception.
            exc_tb: Traceback information.
        """


# pylint: disable-next=too-few-public-methods
class ContextValidator:
    """Base class for validating contexts."""

    # pylint: disable=unused-argument

    @classmethod
    def validate(cls, context_manager: WestpaContextManager) -> bool:
        """Validate contexts for westpa.

        Args:
            context_manager: A westpa context manager to validate.

        Returns:
            If the context is valid.
        """
        logger.info("Validating with: {}", cls.__name__)
        is_valid = True
        context = context_manager.get()
        for key, value in context.items():
            try:
                checker = getattr(cls, key)
            except AttributeError:
                logger.debug("Cannot check {}", key)
                continue
            logger.debug("Checking {}", key)
            if value is not None:
                if isinstance(checker, tuple):
                    if value not in checker:
                        logger.error("  Invalid: {}", value)
                        is_valid = False
                    else:
                        logger.debug("  Valid: {}", value)
                if callable(checker):
                    is_value_valid = checker(value, context)
                    if not is_value_valid:
                        logger.error("  Invalid: {}", value)
                        is_valid = False
                    else:
                        logger.debug("  Valid: {}", value)
            else:
                logger.debug("  Skipping: None")
        return is_valid

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
