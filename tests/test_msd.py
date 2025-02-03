import znmdakit


def test_msd(universe):
    msd = znmdakit.EinsteinMSD(
        universe=universe, select="name H", timestep=0.1, sampling_rate=10
    )
    msd.run()
    assert len(msd.results) == len(universe.trajectory)
