import znmdakit


def test_workflow(repo_path, lammps_npt):
    project = znmdakit.Project()

    with project:
        u = znmdakit.Universe(lammps_npt)
        rdf = znmdakit.RDF(universe=u, select="name H", nbins=100)
        msd = znmdakit.EinsteinMSD(
            universe=u, select="name H", timestep=0.1, sampling_rate=10
        )

    project.run()

    rdf.load()
    msd.load()

    print(rdf.results)
    print(msd.results)
