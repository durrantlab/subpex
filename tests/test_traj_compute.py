import os
import shutil

import MDAnalysis as mda

from subpex.data.compute import run_compute_data
from subpex.pocket.detect import get_fop_pocket_convenience


def test_traj_data_m7g(path_m7g_paths, m7g_config, tmp_dir):
    u = mda.Universe(path_m7g_paths["ref_pdb"])
    fop_ref = get_fop_pocket_convenience(u.atoms, m7g_config)

    write_dir = os.path.join(tmp_dir, "m7g/calc-data/")
    if os.path.exists(write_dir):
        shutil.rmtree(write_dir)
    os.makedirs(write_dir, exist_ok=False)

    m7g_config.calculated_points = 2

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
