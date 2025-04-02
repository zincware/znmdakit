import numpy as np
from MDAnalysis import AtomGroup
from MDAnalysis.transformations import TransformationBase
from tqdm import tqdm


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


# https://gist.github.com/orbeckst/b611626e7640f399580e233195c6429c
# Center of Mass Transformation

# TODO: consider using the "BeadGroup" suggested in the gist instead


class COMTransform(TransformationBase):
    def __init__(self, reference: AtomGroup, name: str):
        self.reference = reference
        self.com_atoms = reference.select_atoms(f"name {name}")

        # sanity check
        a = self.get_com().shape
        b = self.com_atoms.positions.shape
        if a != b:
            raise ValueError(f"Shape mismatch: {a} != {b}")

    def get_com(self):
        return self.reference.center_of_mass(unwrap=True, compound="fragments")

    def __call__(self, ts):
        self.com_atoms.positions = self.get_com()
        return ts


def get_com_transform(universe) -> list:
    # TODO: sanity check, write to XYZ, construct a base Node and one node that just writes the transformed trajectory to disk!
    # when writing to disk, allow for a list of selection strings
    transformations = {}
    for residue in tqdm(universe.residues, desc="Updating atom names"):
        if len(residue.atoms) > 1:
            # we will later update the positions of the first atom in each residue
            residue.atoms[0].name = "COM"
    for residue in tqdm(universe.residues, desc="Preparing COMTransform"):
        if len(residue.atoms) > 1:
            if residue.resname not in transformations:
                transformations[residue.resname] = COMTransform(
                    reference=universe.select_atoms(f"resname {residue.resname}"),
                    name="COM",
                )
    return list(transformations.values())
