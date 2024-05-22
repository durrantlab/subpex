import os
import shutil

from subpex.configs import WestpaConfig
from subpex.setup.westpa.main import run_westpa_setup


def test_westpa_init_larp1_dm15(path_larp1_dm15_paths, larp1_dm15_config, tmp_dir):
    westpa_config = WestpaConfig()
    write_dir = os.path.join(tmp_dir, "westpa_init_larp1_dm15")
    if os.path.exists(write_dir):
        shutil.rmtree(write_dir)

    # Setup WESTPA
    westpa_config.env.load_modules = ["gcc/10.2.0", "openmpi/4.1.1", "amber/22"]

    run_westpa_setup(
        subpex_config=larp1_dm15_config,
        westpa_config=westpa_config,
        write_dir=write_dir,
        overwrite=True,
    )
