import logging
from pathlib import Path

import zntrack
from rdkit2ase import smiles2atoms

from znmdakit.utils import get_universe

log = logging.getLogger(__name__)


class Universe(zntrack.Node):
    data_file: str = zntrack.deps_path()
    residues: dict[str, str] = zntrack.params()  # dict[id, smiles]

    def run(self):
        print(self.universe)

    @property  # cached property needed?
    def universe(self):
        log.critical("Constructing universe")
        universe = get_universe(Path(self.data_file))

        residues = {k: smiles2atoms(smiles=v) for k, v in self.residues.items()}

        mda_residues = {}

        for idx, residue in enumerate(residues):
            mda_residues[residue] = universe.add_Residue(
                resid=idx,
                resname=residue,
            )

        if len(self.residues) > 0:
            for mol in universe.atoms.fragments:
                for residue in residues:
                    # We assign the residues by matching the chemical symbols of the atoms
                    # There is no graph comparison to ensure the atoms are connected in the same way!
                    # If the connectivity graph changes over time, this is not captured and thus not suitable
                    # for reactive systems!
                    if sorted(mol.names) == sorted(
                        residues[residue].get_chemical_symbols()
                    ):
                        # https://docs.mdanalysis.org/stable/documentation_pages/core/universe.html#MDAnalysis.core.universe.Universe.add_Residue
                        universe.atoms[mol.indices].residues = mda_residues[residue]
                        break
                else:
                    log.warning(
                        f"Could not find residue for molecule {sorted(mol.names)} in residues {list(sorted(x.get_chemical_symbols()) for x in residues.values())}"
                    )

        return universe
