import MDAnalysis as mda

from subpex.pocket.detect import get_fop_pocket, get_pocket_selection


def test_pocket_selection_m7g(path_m7g_ref_pdb):
    # Manual parameters
    pocket_center = [30.037991, 41.463818, 30.3776]
    pocket_radius = 6.5
    distance_constraint = 6.7

    # Load structure
    u = mda.Universe(path_m7g_ref_pdb)

    pocket_selection = get_pocket_selection(
        u,
        center=pocket_center,
        radius=pocket_radius,
        distance_constraint=distance_constraint,
    )
    ref_string = "resid 124 or resid 125 or resid 128 or resid 88 or resid 89 or resid 121 or resid 92 or resid 129 and (not name H*)"
    assert pocket_selection == ref_string


def test_pocket_fop_gen_m7g(path_m7g_ref_pdb):
    # Manual parameters
    pocket_center = [30.037991, 41.463818, 30.3776]
    pocket_radius = 6.5
    pocket_resolution = 0.5
    pocket_selection_str = "resid 124 or resid 125 or resid 128 or resid 88 or resid 89 or resid 121 or resid 92 or resid 129 and (not name H*)"
    # Load structure
    u = mda.Universe(path_m7g_ref_pdb)
    protein_coords = u.select_atoms("protein").positions

    pocket_calpha = u.select_atoms(pocket_selection_str + " and name CA*").positions

    fop_frame = get_fop_pocket(
        protein_coords,
        pocket_calpha,
        pocket_center,
        pocket_resolution,
        pocket_radius,
    )

    assert len(fop_frame) == 1057
    assert fop_frame[57] == [31.5, 38.5, 31.0]
