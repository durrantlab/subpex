from collections.abc import Sequence

from pydantic import BaseModel, Field


class ExeChildConfig(BaseModel):
    enabled: bool = Field(default=True)
    executable: str
    """Path to executable to run for data."""
    stdout: str | None = Field(default=None)
    stderr: str | None = Field(default=None)


class DatasetsConfig(BaseModel):
    name: str


class ExecutableConfig(BaseModel):
    datasets: Sequence[DatasetsConfig] = Field(
        default=[
            DatasetsConfig(name="pocket_rmsd"),
            DatasetsConfig(name="pocket_volume"),
            DatasetsConfig(name="backbone_rmsd"),
            DatasetsConfig(name="pocket_rog"),
            DatasetsConfig(name="pocket_jd"),
        ]
    )
    pre_iteration: ExeChildConfig = Field(
        default=ExeChildConfig(
            enabled=False,
            executable=r"$WEST_SIM_ROOT/westpa_scripts/pre_iter.sh",
            stderr=r"$WEST_SIM_ROOT/job_logs/pre_iter.log",
        )
    )
    """Pre-iteration executable settings."""
    post_iteration: ExeChildConfig = Field(
        default=ExeChildConfig(
            enabled=True,
            executable=r"$WEST_SIM_ROOT/westpa_scripts/post_iter.sh",
            stderr=r"$WEST_SIM_ROOT/job_logs/post_iter.log",
        )
    )
    """Post-iteration executable settings."""
    gen_istate: ExeChildConfig = Field(
        default=ExeChildConfig(
            enabled=True,
            executable=r"$WEST_SIM_ROOT/westpa_scripts/gen_istate.sh",
            stdout=r"/dev/null",
            stderr=r"$WEST_SIM_ROOT/job_logs/gen_istate.log",
        )
    )
    """Executable to generate initial state."""
    get_pcoord: ExeChildConfig = Field(
        default=ExeChildConfig(
            enabled=True,
            executable=r"$WEST_SIM_ROOT/westpa_scripts/get_pcoord.sh",
            stdout=r"/dev/null",
            stderr=r"$WEST_SIM_ROOT/job_logs/get_pcoord.log",
        )
    )
    """Executable to get progress coordinate."""
