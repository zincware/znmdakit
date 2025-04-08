from pathlib import Path

import ase
import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
import seaborn as sns
import zntrack
from MDAnalysis import Universe
from MDAnalysis.analysis import rdf

from znmdakit.transformations import get_com_transform
from znmdakit.utils import ComparisonResults


class InterRDF(zntrack.Node):
    universe: Universe = zntrack.deps()
    g1: str = zntrack.params()
    g2: str = zntrack.params()
    nbins: int = zntrack.params()
    range: str | tuple = zntrack.params("auto")

    apply_com_transform: bool = zntrack.params(False)
    # replace the position of the first atom in
    # each residuewith the center of mass and name it COM

    figure_path: Path = zntrack.outs_path(zntrack.nwd / "rdf.png")

    results: pd.DataFrame = zntrack.plots(y="g(r)", x="r")

    def run(self):
        universe = self.universe

        if self.apply_com_transform:
            transformations = get_com_transform(universe)
            universe.trajectory.add_transformations(*transformations)

        radial_dist_fn = rdf.InterRDF(
            universe.select_atoms(self.g1),
            universe.select_atoms(self.g2),
            nbins=1000,
            range=(0, universe.dimensions[0] / 2)
            if self.range == "auto"
            else self.range,
        )
        radial_dist_fn.run(verbose=True)

        self.results = pd.DataFrame(
            {
                "r": radial_dist_fn.results.bins[1:],
                "g(r)": radial_dist_fn.results.rdf[1:],
            }
        )

        self.results.set_index("r", inplace=True)
        self.figure_path.parent.mkdir(parents=True, exist_ok=True)
        figure = self.get_figure()
        figure.savefig(self.figure_path, bbox_inches="tight")

    def get_figure(self):
        sns.set_theme(style="whitegrid")  # Use a clean seaborn style

        fig, ax = plt.subplots(
            figsize=(6, 4), dpi=100
        )  # Set figure size and resolution
        ax.plot(self.results.index, self.results["g(r)"], color="tab:blue", lw=2)

        ax.set_ylabel("Radial Distribution Function, g(r)", fontsize=12)
        ax.set_xlabel("Distance, r (Å)", fontsize=12)
        ax.set_title(f"InterRDF: {self.g1} - {self.g2}", fontsize=14, fontweight="bold")

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        return fig

    @staticmethod
    def compare(*nodes: "InterRDF") -> ComparisonResults:
        frames = [ase.Atoms()]
        fig = go.Figure()
        for node in nodes:
            name = node.name.replace(f"_{node.__class__.__name__}", "")
            fig.add_trace(
                go.Scatter(
                    x=node.results.index,
                    y=node.results["g(r)"],
                    mode="lines",
                    name=name,
                )
            )

        title = f"InterRDF: {nodes[0].g1} - {nodes[0].g2}"
        fig.update_layout(
            title=title,
            xaxis_title="Distance r / Å",
            yaxis_title="g(r)",
            legend_title="Models",
            plot_bgcolor="rgba(0, 0, 0, 0)",
            paper_bgcolor="rgba(0, 0, 0, 0)",
        )
        fig.update_xaxes(
            showgrid=True,
            gridwidth=1,
            gridcolor="rgba(120, 120, 120, 0.3)",
            zeroline=False,
        )
        fig.update_yaxes(
            showgrid=True,
            gridwidth=1,
            gridcolor="rgba(120, 120, 120, 0.3)",
            zeroline=False,
        )

        return ComparisonResults(
            frames=frames,
            figures={title: fig},
        )
