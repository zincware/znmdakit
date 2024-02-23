from MDAnalysis.coordinates.base import ReaderBase, Timestep
from MDAnalysis.topology.base import TopologyReaderBase

from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyattrs import (
    Atomnames,
    Atomids,
    Atomtypes,
    Masses,
    Resids,
    Resnums,
    Segids,
    Elements,
)
import ase.io
import numpy as np

# TODO: dont name it ASEReader, also support H5MD through ZnH5MD


class ASEReader(ReaderBase):
    units = {"time": "ps", "length": "Angstrom"}

    def __init__(self, filename, convert_units=True, **kwargs):
        super().__init__(filename, convert_units, **kwargs)
        self.ts = Timestep(self.n_atoms, **self._ts_kwargs)

        self._reader = iter(ase.io.iread(self.filename))
        self._read_next_timestep()
        # reset
        self._reader = iter(ase.io.iread(self.filename))

    def _read_next_timestep(self, ts=None) -> Timestep:
        if ts is not None:
            raise NotImplementedError("Parallel reading not implemented")
        ts = self.ts

        ts.frame += 1
        atoms = next(self._reader)
        ts.positions = atoms.get_positions()
        # set cell
        ts.dimensions = atoms.get_cell().cellpar()
        return ts

    def _reopen(self):
        pass  # abstract method

    @property
    def n_atoms(self):
        return ase.io.read(self.filename, index=0).get_number_of_atoms()

    @property
    def n_frames(self):
        return len(list(ase.io.iread(self.filename)))


class ASETopologyReader(TopologyReaderBase):
    def parse(self, **kwargs):
        """Read the file and return the structure.

        Returns
        -------
        MDAnalysis Topology object
        """
        atoms = ase.io.read(self.filename, index=0)
        natoms = atoms.get_number_of_atoms()
        names = atoms.get_chemical_symbols()
        atomtypes = atoms.get_chemical_symbols()
        masses = atoms.get_masses()

        attrs = [
            Atomnames(names),
            Atomids(np.arange(natoms) + 1),
            Atomtypes(atomtypes, guessed=True),
            Masses(masses, guessed=True),
            Resids(np.array([1])),
            Resnums(np.array([1])),
            Segids(np.array(["SYSTEM"], dtype=object)),
            Elements(names),
        ]

        top = Topology(natoms, 1, 1, attrs=attrs)

        return top
