[![ZnTrack](https://img.shields.io/badge/Powered%20by-ZnTrack-%23007CB0)](https://zntrack.readthedocs.io/en/latest/)
[![zincware](https://img.shields.io/badge/Powered%20by-zincware-darkcyan)](https://github.com/zincware)

# ZnMDAKit

[ZnTrack](https://github.com/zincware/ZnTrack) enabled [MDAnalysis Toolkits](https://mdakits.mdanalysis.org/) with full parameter and result tracking capabilities.

## Example

```python
import znmdakit

project = znmdakit.Project()

with project:
    u = znmdakit.Universe("trajectory.h5")
    rdf = znmdakit.RDF(universe=u, select="name H", nbins=100)
    msd = znmdakit.EinsteinMSD(universe=u, select="name H", timestep=0.1, sampling_rate=10)

project.run()

rdf.load()
msd.load()

print(rdf.results)
print(msd.results)
```

## Available Nodes
- `znmdakit.RDF`
- `znmdakit.EinsteinMSD`
