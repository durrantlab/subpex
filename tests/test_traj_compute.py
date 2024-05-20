import os
import shutil

import MDAnalysis as mda
import numpy as np

from subpex.data.compute import run_compute_data
from subpex.pocket.detect import get_fop_pocket_convenience


def test_traj_data_m7g(path_m7g_paths, m7g_config, tmp_dir):
    """
    Test the computation of auxiliary data during a simulation for the binding pocket descriptors.

    This function tests the framework's ability to compute various descriptors for a binding pocket
    in a molecular dynamics simulation. The test uses a predefined PDB reference structure,
    topology, and trajectory files, and verifies the computed auxiliary data against known values.
    The reference structure is taken from the last frame of the same simulation, so properties
    compared to the reference should be close to zero.

    Args:
        path_m7g_paths (dict): Dictionary containing paths to the input files required for the test:
            - "ref_pdb": Path to the reference PDB file.
            - "topo": Path to the topology file.
            - "traj": Path to the trajectory file.
        m7g_config (object): Configuration object for the simulation which includes settings and
            parameters for computing the auxiliary data.
        tmp_dir (str): Path to a temporary directory for writing output data during the test.

    Steps:
        1. Load the reference PDB file using MDAnalysis.
        2. Detect the binding pocket using the `get_fop_pocket_convenience` function.
        3. Prepare the directory for writing computed data.
        4. Set the number of calculated points in the configuration.
        5. Run the data computation using `run_compute_data`.
        6. Verify the computed auxiliary data (backbone RMSD, pocket JD, pocket RMSD,
           pocket radius of gyration, and pocket volume) against expected values using
           `np.allclose`.

    Assertions:
        - Check if the computed backbone RMSD values match the expected values.
        - Check if the computed pocket JD values match the expected values.
        - Check if the computed pocket RMSD values match the expected values.
        - Check if the computed pocket radius of gyration values match the expected values.
        - Check if the computed pocket volume values match the expected values.

    Raises:
        AssertionError: If any of the computed values do not match the expected values.

    Notes:
        This test ensures that the framework accurately computes auxiliary data for the
        binding pocket during a molecular dynamics simulation, which is crucial for analyzing
        the stability and characteristics of the pocket over time. Since the reference structure
        is the last frame of the same simulation, properties compared to the reference should
        be close to zero.
    """
    u = mda.Universe(path_m7g_paths["ref_pdb"])
    fop_ref = get_fop_pocket_convenience(u.atoms, m7g_config)

    write_dir = os.path.join(tmp_dir, "m7g/calc-data/")
    if os.path.exists(write_dir):
        shutil.rmtree(write_dir)
    os.makedirs(write_dir, exist_ok=False)

    m7g_config.calculated_points = 3

    run_compute_data(
        topo_path=path_m7g_paths["topo"],
        traj_path=path_m7g_paths["traj"],
        ref_path=path_m7g_paths["ref_pdb"],
        fop_ref=fop_ref,
        subpex_config=m7g_config,
        selection_align_suffix=" and backbone",
        write=False,
        write_dir=write_dir,
    )

    assert np.allclose(
        np.array([1.1303496360778809, 0.9975303411483765, 0.00028853819821961224]),
        np.array(m7g_config.data.aux.backbone_rmsd.values),
    )
    assert np.allclose(
        np.array([0.16315789473684206, 0.2279226240538268, 0.0]),
        np.array(m7g_config.data.aux.pocket_jd.values),
    )
    assert np.allclose(
        np.array([1.5243936777114868, 1.6841306686401367, 0.0003748516319319606]),
        np.array(m7g_config.data.aux.pocket_rmsd.values),
    )
    assert np.allclose(
        np.array([3.1923352274061543, 3.1671394736443172, 3.2815136334233626]),
        np.array(m7g_config.data.aux.pocket_rog.values),
    )
    assert np.allclose(
        np.array([129.625, 131.0, 131.125]),
        np.array(m7g_config.data.aux.pocket_volume.values),
    )
