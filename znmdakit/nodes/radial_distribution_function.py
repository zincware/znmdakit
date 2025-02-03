import pandas as pd
import zntrack
from MDAnalysis import Universe
from MDAnalysis.analysis import rdf


class InterRDF(zntrack.Node):
    universe: Universe = zntrack.deps()
    g1: str = zntrack.params()
    g2: str = zntrack.params()
    nbins: int = zntrack.params()
    range: str | tuple = zntrack.params("auto")

    results: pd.DataFrame = zntrack.plots(y="g(r)", x="r")

    def run(self):
        universe = self.universe

        RDF = rdf.InterRDF(
            universe.select_atoms(self.g1),
            universe.select_atoms(self.g2),
            nbins=1000,
            range=(0, universe.dimensions[0] / 2)
            if self.range == "auto"
            else self.range,
        )
        RDF.run(verbose=True)

        self.results = pd.DataFrame(
            {"r": RDF.results.bins[1:], "g(r)": RDF.results.rdf[1:]}
        )

        self.results.set_index("r", inplace=True)
