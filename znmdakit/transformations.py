import numpy as np
from MDAnalysis.transformations import TransformationBase


class UnWrap(TransformationBase):
    prev = None
    # def __init__(self, **kwargs):
    #     super().__init__(parallelizable=False, **kwargs)
    #     self.prev = None

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
