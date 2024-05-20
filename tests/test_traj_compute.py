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
        np.array(
            [
                1.1303496360778809,
                1.1373416185379028,
                1.0807933807373047,
                0.9600823521614075,
                0.8374161720275879,
            ]
        ),
        np.array(m7g_config.data.aux.backbone_rmsd.values),
    )
    assert np.allclose(
        np.array(
            [
                0.16315789473684206,
                0.28441127694859036,
                0.2343059239610964,
                0.22831858407079642,
                0.24537037037037035,
            ]
        ),
        np.array(m7g_config.data.aux.pocket_jd.values),
    )
    assert np.allclose(
        np.array(
            [
                1.5243936777114868,
                1.889655351638794,
                2.2993032932281494,
                1.395134687423706,
                1.8448500633239746,
            ]
        ),
        np.array(m7g_config.data.aux.pocket_rmsd.values),
    )
    assert np.allclose(
        np.array(
            [
                3.1923352274061543,
                3.1152869440628796,
                3.1127554945582774,
                3.10684407247731,
                3.051189269768271,
            ]
        ),
        np.array(m7g_config.data.aux.pocket_rog.values),
    )
    assert np.allclose(
        np.array([129.625, 126.25, 117.875, 121.5, 109.125]),
        np.array(m7g_config.data.aux.pocket_volume.values),
    )
