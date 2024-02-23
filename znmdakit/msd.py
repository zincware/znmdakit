import zntrack
import pandas as pd
import numpy as np
import tqdm
import matplotlib.pyplot as plt
from MDAnalysis.transformations import TransformationBase
import pint
import typing as t


class UnWrap(TransformationBase):
    def __init__(self, **kwargs):
        super().__init__(parallelizable=False, **kwargs)
        self.prev = None

    def _transform(self, ts):
        if ts.frame == 0:
            self.prev = ts.positions.copy()
            return ts

        assert np.all(ts.dimensions[3:] == [90, 90, 90])

        while np.any(ts.positions - self.prev > 0.5 * ts.dimensions[:3]):
            ts.positions = np.where(
                ts.positions - self.prev > 0.5 * ts.dimensions[:3],
                ts.positions - ts.dimensions[:3],
                ts.positions,
            )

        while np.any(ts.positions - self.prev < -0.5 * ts.dimensions[:3]):
            ts.positions = np.where(
                ts.positions - self.prev < -0.5 * ts.dimensions[:3],
                ts.positions + ts.dimensions[:3],
                ts.positions,
            )
        # assert that there is no jump larger than half the box size
        assert np.all(np.abs(ts.positions - self.prev) < 0.5 * ts.dimensions[:3])

        self.prev = ts.positions.copy()
        return ts


class EinsteinMSD(zntrack.Node):
    universe: zntrack.Node = zntrack.deps()
    select: str = zntrack.params()
    results: pd.DataFrame = zntrack.plots()
    timestep: float = zntrack.params(None)
    sampling_rate: int = zntrack.params(None)


    def run(self) -> None:
        import MDAnalysis.analysis.msd as msd
        import MDAnalysis as mda

        u = self.universe if isinstance(self.universe, mda.Universe) else self.universe.get_universe()
        
        u.trajectory.add_transformations(UnWrap())

        MSD = msd.EinsteinMSD(u, select=self.select, msd_type="xyz", fft=True)
        MSD.run(verbose=True)

        if self.timestep is not None and self.sampling_rate is not None:
            dt = self.timestep * self.sampling_rate
        else:
            dt = u.trajectory.dt

        lagtimes = np.arange(MSD.n_frames) * dt * u.trajectory.n_frames
        
        self.results = pd.DataFrame(
            {
                "lagtimes": lagtimes,
                "msd": MSD.results.timeseries,
            }
        )
        # set the index to lagtimes
        self.results.set_index("lagtimes", inplace=True)

