import zntrack
import pandas as pd


class RDF(zntrack.Node):
    universe: zntrack.Node = zntrack.deps()
    select: str | list | tuple = zntrack.params()
    nbins: int = zntrack.params()
    results: pd.DataFrame = zntrack.plots()

    def run(self) -> None:
        # TODO: wrap coordinates into box
        import MDAnalysis.analysis.rdf as rdf
        import MDAnalysis as mda

        u = (
            self.universe
            if isinstance(self.universe, mda.Universe)
            else self.universe.get_universe()
        )
        if isinstance(self.select, str):
            a1 = u.select_atoms(self.select)
            a2 = u.select_atoms(self.select)
        else:
            a1 = u.select_atoms(self.select[0])
            a2 = u.select_atoms(self.select[1])

        RDF = rdf.InterRDF(a1, a2, nbins=self.nbins, range=(0, u.dimensions[0] / 2))
        RDF.run(verbose=True)

        self.results = pd.DataFrame(
            {"r": RDF.results.bins[1:], "g(r)": RDF.results.rdf[1:]}
        )
        # set the index to the r values
        self.results.set_index("r", inplace=True)
