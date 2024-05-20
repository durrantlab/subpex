import MDAnalysis as mda

from subpex.pocket.detect import get_fop_pocket_convenience, get_pocket_selection


def test_pocket_selection_m7g(path_m7g_paths, m7g_config):
    # Load structure
    u = mda.Universe(path_m7g_paths["ref_pdb"])

    pocket_selection = get_pocket_selection(
        u,
        center=m7g_config.pocket.center,
        radius=m7g_config.pocket.radius,
        distance_constraint=m7g_config.pocket.selection_dist,
    )
    assert pocket_selection == m7g_config.pocket.selection_str


def test_pocket_fop_gen_m7g(path_m7g_paths, m7g_config):
    u = mda.Universe(path_m7g_paths["ref_pdb"])
    fop_frame = get_fop_pocket_convenience(u.atoms, m7g_config)

    assert len(fop_frame) == 1057
    assert fop_frame[57] == [31.5, 38.5, 31.0]
