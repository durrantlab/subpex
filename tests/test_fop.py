from subpex.fop.gen import gen_fop


def test_fop_gen():
    fop, fop_center = gen_fop(
        center=[0.0, 0.0, 0.0],
        resolution=0.5,
        radius=10.0,
        round_decimals=3,
    )
    assert len(fop) == 33371
    assert fop[5] == [-2.5, -9.5, -0.5]
    assert fop_center == [0.0, 0.0, 0.0]
