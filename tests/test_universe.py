import numpy.testing as npt


def test_universe(universe):
    npt.assert_array_almost_equal(
        universe.dimensions, [32.65, 32.65, 32.65, 90, 90, 90], decimal=2
    )
