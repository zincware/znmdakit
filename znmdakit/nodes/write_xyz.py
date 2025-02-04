import zntrack
from MDAnalysis import Universe
from pathlib import Path
from znmdakit.transformations import UnWrap, get_com_transform
from tqdm import tqdm


class WriteXYZ(zntrack.Node):
    universe: Universe = zntrack.deps()
    frames_path: Path = zntrack.outs_path(zntrack.nwd / "frames.xyz")

    select: str = zntrack.params("all")
    apply_com_transform: bool = zntrack.params(False)
    unwrap: bool = zntrack.params(False)

    def run(self):
        self.frames_path.parent.mkdir(parents=True, exist_ok=True)
        universe = self.universe

        transformations = []
        if self.apply_com_transform:
            transformations.extend(get_com_transform(universe))
        if self.unwrap:
            transformations.append(UnWrap())
        if len(transformations) > 0:
            universe.trajectory.add_transformations(*transformations)

        atoms_group = universe.select_atoms(self.select)
        print(f"Writing {len(atoms_group)} atoms to {self.frames_path}")

        with open(self.frames_path, "w") as f:
            # TODO: write znh5md
            for idx, ts in tqdm(enumerate(universe.trajectory)):
                f.write(f"{len(atoms_group)}\n")
                f.write(f"Frame {ts.frame}\n")
                # TODO: somehow define atom type per residue for ASE?
                for atom in atoms_group:
                    f.write(
                        f"X {atom.position[0]} {atom.position[1]} {atom.position[2]}\n"
                    )
