[![ZnTrack](https://img.shields.io/badge/Powered%20by-ZnTrack-%23007CB0)](https://zntrack.readthedocs.io/en/latest/)
[![zincware](https://img.shields.io/badge/Powered%20by-zincware-darkcyan)](https://github.com/zincware)

# ZnMDAKit

[ZnTrack](https://github.com/zincware/ZnTrack) enabled [MDAnalysis Toolkits](https://mdakits.mdanalysis.org/) with full parameter and result tracking capabilities.

## Example

```python
import znmdakit

project = znmdakit.Project()

with project:

    system = ips.Smiles2Gromacs(
        smiles=[BMIM, BF4],
        density=0.5,
        count=[64, 64],
        labels=["Im", "BF"],
        config_files=mdp_files,
        fudgeLJ=1,
        fudgeQQ=1,
        itp_files=itps,
        pdb_files=pdbs,
        tolerance=1.8,
        cleanup=True,
        maxwarn=0,
    )

    universe = znmdakit.Universe(
        # data_file=znflow.resolve(system.traj_file), # hotfix for https://github.com/zincware/ZnTrack/pull/875
        data=system.frames,
        residues={"BF": BF4, "Im": BMIM},
    )

with project.group("BF4"):
    znmdakit.InterRDF(
        universe=universe.universe,
        g1="name COM and resname BF",
        g2="name COM and resname BF",
        nbins=1000,
        apply_com_transform=True,
    )
    msd = znmdakit.EinsteinMSD(
        universe=universe.universe,
        select="name COM and resname BF",
        timestep=0.001, # ps
        sampling_rate=1000, # ps
        apply_com_transform=True,
    )

    znmdakit.SelfDiffusionFromMSD(
        data=msd,
        start_time=2000,
        end_time=7000,
        always_changed=False,
    )

project.run()

rdf.load()
msd.load()

print(rdf.results)
print(msd.results)
```

## Available Nodes
- `znmdakit.RDF`
- `znmdakit.EinsteinMSD`
