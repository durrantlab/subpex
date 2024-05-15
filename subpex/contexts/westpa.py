from typing import Any

from collections.abc import Iterable, MutableMapping

from loguru import logger

from .base import BaseContextManager, BaseContextValidator


class WestpaContextManager(BaseContextManager):
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

        super().__init__(yaml_paths, **kwargs)


# pylint: disable-next=too-few-public-methods
class WestpaContextValidator(BaseContextValidator):
    """Base class for validating contexts."""

    # pylint: disable=unused-argument
