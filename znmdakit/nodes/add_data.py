import ase.io
import zntrack
from pathlib import Path
import znh5md
import h5py

class AddData(zntrack.Node):
    file: str = zntrack.deps_path()
    frames_path: Path = zntrack.outs_path(zntrack.nwd / "frames.h5")

    def run(self):
        frames = list(ase.io.read(self.file, index=":"))
        io = znh5md.IO(self.frames_path)
        io.extend(frames)
    
    @property
    def frames(self) -> list[ase.Atoms]:
        with self.state.fs.open(self.frames_path, "rb") as f:
            with h5py.File(f, "r") as h5:
                return znh5md.IO(file_handle=h5)[:]
