import znmdakit


def test_rdf(universe):
    rdf = znmdakit.RDF(universe=universe, select="name H", nbins=100)
    rdf.run()
    assert len(rdf.results)
