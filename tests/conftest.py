import os

import pytest

from subpex import enable_logging
from subpex.configs import SubpexConfig

TEST_DIR = os.path.dirname(__file__)


@pytest.fixture(scope="session", autouse=True)
def turn_on_logging():
    enable_logging(10)


@pytest.fixture
def tmp_dir():
    return os.path.join(TEST_DIR, "tmp")


@pytest.fixture
def path_m7g_paths():
    paths = {
        "topo": os.path.join(TEST_DIR, "files/m7g/mol.prmtop"),
        "ref_pdb": os.path.join(TEST_DIR, "files/m7g/equil_npt_last_frame.pdb"),
        "traj": os.path.join(TEST_DIR, "files/m7g/equil_npt.nc"),
    }
    return paths


@pytest.fixture
def m7g_config():
    subpex_config = SubpexConfig()
    subpex_config.pocket.center = [30.037991, 41.463818, 30.3776]
    subpex_config.pocket.radius = 6.5
    subpex_config.pocket.selection_dist = 6.7
    subpex_config.pocket.selection_str = "resid 124 or resid 125 or resid 128 or resid 88 or resid 89 or resid 121 or resid 92 or resid 129 and (not name H*)"
    for aux_data in subpex_config.data.aux:
        aux_data.active = True
    return subpex_config
