import zntrack
import MDAnalysis as mda
from ase.data import chemical_symbols
from MDAnalysis.coordinates.H5MD import H5MDReader
import os
import pathlib
import h5py


# TODO: rename to UniverseNode or so, and have UniverseNode.universe be the actual universe


class Universe(zntrack.Node):
    file: os.PathLike = zntrack.deps_path()


    def get_universe(self) -> mda.Universe:
        if pathlib.Path(self.file).suffix in [".h5", ".hdf5", ".h5md"]:
            with h5py.File(self.file) as f:
                # read the types
                keys = list(f["particles"].keys())
                types = f["particles"][keys[0]]["species/value"][0]
            u = mda.Universe.empty(len(types), trajectory=True)
            reader = H5MDReader(self.file, convert_units=False)
            u.trajectory = reader
            # set atoms types
            u.add_TopologyAttr("types", types)
            # use ase to convert atomic numbers to names
            u.add_TopologyAttr("names", [chemical_symbols[int(t)] for t in types])
        else:
            raise ValueError(f"Unsupported file format for {self.file}")

        return u
    
    # def get_universe(self) -> mda.Universe:
    #     if self.u_kwargs is None:
    #         self.u_kwargs = {}
    #     if self.use_ase:
    #         from .reader import ASEReader, ASETopologyReader

    #         import ase.io

    #         return mda.Universe(
    #             self.files,
    #             format=ASEReader,
    #             topology_format=ASETopologyReader,
    #             **self.u_kwargs,
    #         )

    #     u = mda.Universe(self.files, **self.u_kwargs)
    #     if self.cell_from_ase:
    #         import ase.io

    #         seed_atoms = ase.io.read(self.files, index=0)
    #         u.dimensions = seed_atoms.get_cell().cellpar()
    #     return u
