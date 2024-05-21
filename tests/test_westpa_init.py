import os

from subpex.configs import WestpaConfig
from subpex.setup.westpa.main import run_westpa_setup


def test_westpa_init_m7g(path_m7g_paths, m7g_config, tmp_dir):
    westpa_config = WestpaConfig()
    write_dir = os.path.join(tmp_dir, "westpa_init_m7g")
    run_westpa_setup(
        subpex_config=m7g_config,
        westpa_config=westpa_config,
        write_dir=write_dir,
        overwrite=True,
    )
