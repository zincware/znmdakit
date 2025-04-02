import ase.io
import zntrack
from tqdm import tqdm


class ReadData(zntrack.Node):
    file: str = zntrack.deps_path()

    def run(self):
        pass

    @property
    def frames(self) -> list[ase.Atoms]:
        if self.state.rev is not None or self.state.remote is not None:
            raise ValueError("Cannot read frames from a different rev/remote")
        return list(
            tqdm(
                ase.io.read(self.file, index=":"),
                desc="Reading frames",
                unit="frame",
                ncols=80,
            )
        )
