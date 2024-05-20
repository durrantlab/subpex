from subpex.fop.gen import gen_fop


def test_fop_gen():
    """
    Test the generation of a field of points (FOP) centered around a specified point.

    This function tests the `gen_fop` function's ability to generate a field of points with a
    given resolution and radius, centered at a specified point. The number of generated points,
    the coordinates of a specific point, and the center of the field are verified against expected
    values.

    Steps:
        1. Generate a field of points (FOP) using `gen_fop` with specified parameters:
           - center: [0.0, 0.0, 0.0]
           - resolution: 0.5
           - radius: 10.0
           - round_decimals: 3
        2. Verify the number of generated points.
        3. Verify the coordinates of a specific point in the field.
        4. Verify the center of the generated field.

    Raises:
        AssertionError: If the number of generated points, the coordinates of the specific point,
            or the center of the field do not match the expected values.
    """
    fop, fop_center = gen_fop(
        center=[0.0, 0.0, 0.0],
        resolution=0.5,
        radius=10.0,
        round_decimals=3,
    )
    assert len(fop) == 33371
    assert fop[5] == [-2.5, -9.5, -0.5]
    assert fop_center == [0.0, 0.0, 0.0]
