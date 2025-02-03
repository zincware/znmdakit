import MDAnalysis.analysis.msd as msd
import numpy as np
import pandas as pd
import zntrack
from MDAnalysis import Universe

from znmdakit.transformations import UnWrap, get_com_transform


class EinsteinMSD(zntrack.Node):
    universe: Universe = zntrack.deps()
    select: str = zntrack.params()
    results: pd.DataFrame = zntrack.plots(y="msd", x="lagtimes")
    apply_com_transform: bool = zntrack.params(False)

    # TODO: these two should be read from the h5md file
    timestep: float = zntrack.params(None)
    sampling_rate: int = zntrack.params(None)

    def run(self) -> None:
        universe = self.universe

        transformations = []

        if self.apply_com_transform:
            transformations.extend(get_com_transform(universe))
        universe.trajectory.add_transformations(*transformations, UnWrap()) # UnWrap is the last one!

        MSD = msd.EinsteinMSD(universe, select=self.select, msd_type="xyz", fft=True)
        MSD.run(verbose=True)

        if self.timestep is not None and self.sampling_rate is not None:
            dt = self.timestep * self.sampling_rate
        else:
            dt = universe.trajectory.dt

        lagtimes = np.arange(MSD.n_frames) * dt

        self.results = pd.DataFrame(
            {
                "lagtimes": lagtimes,
                "msd": MSD.results.timeseries,
            }
        )
        # set the index to lagtimes
        self.results.set_index("lagtimes", inplace=True)


class SelfDiffusionFromMSD(zntrack.Node):
    pass
