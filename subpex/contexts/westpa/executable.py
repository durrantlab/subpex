from typing import Any

from collections.abc import MutableMapping


class Executable:
    """Executable settings for WESTPA."""

    def __init__(self):
        self.pre_iteration: MutableMapping[str, Any] = {
            "enabled": False,
            "executable": "$WEST_SIM_ROOT/westpa_scripts/pre_iter.sh",
            "stderr": "$WEST_SIM_ROOT/job_logs/pre_iter.log",
        }
        """Pre-iteration executable settings."""
        self.post_iteration: MutableMapping[str, Any] = {
            "enabled": True,
            "executable": "$WEST_SIM_ROOT/westpa_scripts/post_iter.sh",
            "stderr": "$WEST_SIM_ROOT/job_logs/post_iter.log",
        }
        """Post-iteration executable settings."""
