import os

import pytest

from subpex import enable_logging

TEST_DIR = os.path.dirname(__file__)
TMP_dir = os.path.join(TEST_DIR, "tmp")


@pytest.fixture(scope="session", autouse=True)
def turn_on_logging():
    enable_logging(10)


@pytest.fixture
def path_m7g_paths():
    paths = {
        "topo": os.path.join(TEST_DIR, "files/m7g/mol.prmtop"),
        "ref_pdb": os.path.join(TEST_DIR, "files/m7g/equil_npt_last_frame.pdb"),
        "traj": os.path.join(TEST_DIR, "files/m7g/equil_npt.nc"),
    }
    return paths
