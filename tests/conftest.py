import os
import pathlib
import shutil

import dvc.cli
import git
import h5py
import MDAnalysis as mda
import pytest
from ase.data import chemical_symbols
from MDAnalysis.coordinates.H5MD import H5MDReader

FILE = pathlib.Path(__file__).parent / "lammps_npt.h5"


def get_mda_universe():
    u = mda.Universe.empty(1000, trajectory=True)
    reader = H5MDReader(FILE, convert_units=False)
    u.trajectory = reader.trajectory

    with h5py.File(FILE) as f:
        # read the types
        types = f["particles/all/species/value"][0]

    # set atoms types
    u.add_TopologyAttr("types", types)
    u.add_TopologyAttr("names", [chemical_symbols[int(t)] for t in types])

    return u


def get_znmda_universe():
    import znmdakit

    return znmdakit.Universe(FILE).get_universe()


@pytest.fixture(params=[get_mda_universe, get_znmda_universe])
def universe(request):
    return request.param()


@pytest.fixture
def repo_path(tmp_path, request):
    shutil.copy(request.module.__file__, tmp_path)
    os.chdir(tmp_path)
    git.Repo.init()
    dvc.cli.main(["init"])

    return tmp_path


@pytest.fixture
def lammps_npt():
    return FILE
