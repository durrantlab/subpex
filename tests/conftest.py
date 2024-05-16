import os

import pytest

from subpex import enable_logging

TEST_DIR = os.path.dirname(__file__)


@pytest.fixture(scope="session", autouse=True)
def turn_on_logging():
    enable_logging(10)


@pytest.fixture
def path_m7g_ref_pdb():
    return os.path.join(TEST_DIR, "files/m7g/equil_npt_last_frame.pdb")
