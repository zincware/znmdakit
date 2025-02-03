from pathlib import Path

import pandas as pd
import zntrack


class BulkDensity(zntrack.Node):
    data_file: Path | str = zntrack.deps_path()

    results: pd.DataFrame = zntrack.plots()

    def run(self):
        pass
