from pathlib import Path
import typing as t

import plotly.graph_objects as go
import ase
import MDAnalysis as mda
import numpy as np
import znh5md
from ase.neighborlist import natural_cutoffs
from MDAnalysis.coordinates.H5MD import H5MDReader

FIGURES = t.Dict[str, go.Figure]
FRAMES = t.List[ase.Atoms]

class ComparisonResults(t.TypedDict):
    frames: FRAMES
    figures: FIGURES


def get_bonds(atoms: ase.Atoms, mult: float = 1.2) -> list[tuple[int, int]]:
    """Calculate the bonds in an ASE Atoms object.

    Parameters
    ----------
    atoms : ase.Atoms
        The ASE Atoms object for which to calculate bonds.
    mult : float
        The multiplier for the cutoff distance.

    Returns
    -------
    list of tuple of int
        A list of tuples, each containing two indices of bonded atoms.
    """
    if not all(atoms.pbc):
        raise ValueError("Periodic boundary conditions must be enabled for all axes.")

    # Compute distance matrix with minimum image convention
    distance_matrix = atoms.get_all_distances(mic=True)
    np.fill_diagonal(distance_matrix, np.inf)  # Ignore self-distances

    # Compute cutoff distances
    cutoffs = np.array(natural_cutoffs(atoms, mult=mult))
    cutoff_matrix = cutoffs[:, None] + cutoffs[None, :]

    # Find bonded pairs
    bonded_indices = np.argwhere(distance_matrix <= cutoff_matrix)

    # Ensure i < j to avoid duplicate bonds
    bonds = [(i, j) for i, j in bonded_indices if i < j]

    return bonds


def get_universe(file: Path) -> mda.Universe:
    """Create a MDAnalysis Universe from a H5MD file.

    Parameters
    ----------
    file : Path
        The path to the H5MD file.

    Returns
    -------
    mda.Universe
        The MDAnalysis Universe object.
    """

    if file.suffix not in [".h5", ".h5md", ".hdf5"]:
        raise ValueError("Currently, only HDF5 files are supported.")

    atoms: ase.Atoms = znh5md.IO(file)[0]
    # consider https://docs.mdanalysis.org/stable/documentation_pages/guesser_modules/default_guesser.html#MDAnalysis.guesser.default_guesser.DefaultGuesser.guess_bonds
    bonds = get_bonds(atoms)

    universe = mda.Universe.empty(n_atoms=len(atoms), trajectory=True)
    # https://userguide.mdanalysis.org/stable/examples/constructing_universe.html#Adding-topology-attributes
    universe.add_TopologyAttr("names", atoms.get_chemical_symbols())
    universe.add_TopologyAttr("type", atoms.get_chemical_symbols())
    universe.add_TopologyAttr("masses", atoms.get_masses())
    universe.add_TopologyAttr("bonds", bonds)

    # TODO: warning if the box changes in time! We do not support NPT yet!
    universe.dimensions = atoms.cell.cellpar()

    # add new unit translation
    # TODO: we assume these units, might not fit!
    H5MDReader._unit_translation["force"].update({"eV/Angstrom": "kJ/(mol*Angstrom)"})
    H5MDReader._unit_translation["velocity"].update({"Angstrom/fs": "Angstrom/fs"})

    reader = H5MDReader(file, convert_units=False)
    universe.trajectory = reader

    return universe
