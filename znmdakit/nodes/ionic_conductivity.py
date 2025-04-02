import typing as t

import pint
import zntrack
from MDAnalysis import Universe

ureg = pint.UnitRegistry()


class NEIons(t.TypedDict):
    diffusion: float
    charge: int
    select: str


class NEParams(t.TypedDict):
    temperature: float
    ions: t.Dict[str, NEIons]


class NernstEinsteinIonicConductivity(zntrack.Node):
    """Compute the ionic conductivity using the Nernst-Einstein equation.


    References
    ----------
    https://www.sciencedirect.com/science/article/pii/S2772422024000120
    https://jcheminf.biomedcentral.com/articles/10.1186/s13321-023-00687-y
    https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.122.136001

    """

    universe: Universe = zntrack.deps()
    params: NEParams = zntrack.params()

    metrics: dict = zntrack.metrics()

    def run(self):
        volume = (
            self.universe.dimensions[0]
            * self.universe.dimensions[1]
            * self.universe.dimensions[2]
        ) * ureg.angstrom**3

        # e^2 / (T * kB * V)
        prefactor = (ureg.elementary_charge**2) / (
            self.params["temperature"] * ureg.kelvin * ureg.boltzmann_constant * volume
        )
        # is it 3 * volume or just volume?

        diff_x_charge = []
        for ion_data in self.params["ions"].values():
            atoms_group = self.universe.select_atoms(ion_data["select"])

            # D * z^2 per ion (thus * n_ions)
            n_ions = len(atoms_group.fragments)
            diff = ion_data["diffusion"] * ureg.angstrom**2 / ureg.nanosecond
            value = diff * ion_data["charge"] ** 2 * n_ions

            diff_x_charge.append(value)

        sigma_NE = prefactor * sum(diff_x_charge)
        sigma_NE = sigma_NE.to("S/m")

        self.metrics = {"Nernst-Einstein ionic conductivity": sigma_NE.magnitude}
