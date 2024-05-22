import os
import shutil

import MDAnalysis as mda
import numpy as np

from subpex.data.compute import run_compute_data
from subpex.pocket.detect import get_fop_pocket_convenience


def test_traj_data_larp1_dm15(path_larp1_dm15_paths, larp1_dm15_config, tmp_dir):
    """
    Test the computation of auxiliary data during a simulation for the binding pocket descriptors.

    This function tests the framework's ability to compute various descriptors for a binding pocket
    in a molecular dynamics simulation. The test uses a predefined PDB reference structure,
    topology, and trajectory files, and verifies the computed auxiliary data against known values.

    Args:
        path_larp1_dm15_paths (dict): Dictionary containing paths to the input files required for the test:
            - "ref_pdb": Path to the reference PDB file.
            - "topo": Path to the topology file.
            - "traj": Path to the trajectory file.
        larp1_dm15_config (object): Configuration object for the simulation which includes settings and
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
    u = mda.Universe(path_larp1_dm15_paths["ref_pdb"])
    fop_ref = get_fop_pocket_convenience(u.atoms, larp1_dm15_config)

    write_dir = os.path.join(tmp_dir, "larp1_dm15/calc-data/")
    if os.path.exists(write_dir):
        shutil.rmtree(write_dir)
    os.makedirs(write_dir, exist_ok=False)

    run_compute_data(
        topo_path=path_larp1_dm15_paths["topo"],
        traj_path=path_larp1_dm15_paths["traj"],
        ref_path=path_larp1_dm15_paths["ref_pdb"],
        fop_ref=fop_ref,
        subpex_config=larp1_dm15_config,
        selection_align_suffix=" and backbone",
        write=False,
        write_dir=write_dir,
    )

    assert np.allclose(
        np.array(
            [
                1.1303496360778809,
                1.0551434755325317,
                1.0171629190444946,
                0.9678191542625427,
                0.0002885486464947462,
            ]
        ),
        np.array(larp1_dm15_config.data.aux.backbone_rmsd.values),
    )
    assert np.allclose(
        np.array(
            [
                0.16315789473684206,
                0.2139737991266376,
                0.25634517766497467,
                0.26182432432432434,
                0.0,
            ]
        ),
        np.array(larp1_dm15_config.data.aux.pocket_jd.values),
    )
    assert np.allclose(
        np.array(
            [
                1.5243936777114868,
                1.9039777517318726,
                1.5977983474731445,
                2.2219390869140625,
                0.0003747374867089093,
            ]
        ),
        np.array(larp1_dm15_config.data.aux.pocket_rmsd.values),
    )
    assert np.allclose(
        np.array(
            [
                3.1923352274061543,
                3.12051673064247,
                3.127790000499829,
                3.05571329907486,
                3.2499845795814117,
            ]
        ),
        np.array(larp1_dm15_config.data.aux.pocket_rog.values),
    )
    assert np.allclose(
        np.array([129.625, 123.25, 126.25, 124.625, 129.0]),
        np.array(larp1_dm15_config.data.aux.pocket_volume.values),
    )
