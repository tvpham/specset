import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from specset.utils import _psm_stream_to_specset


@pytest.fixture
def psm_factory():
    """Return a factory building PSM dicts accepted by the specset writer."""

    def make_psm(pep, z=2, rt=1.0, im=0.1, matches="b2;y5", intensities="100;200", **extra):
        psm = {"pep": pep, "z": z, "rt": rt, "im": im, "matches": matches, "intensities": intensities}
        psm.update(extra)
        return psm

    return make_psm


@pytest.fixture
def specset_builder():
    """Return a helper writing a specset directory from a list of PSM dicts."""

    def build(path, psms, **kwargs):
        _psm_stream_to_specset(iter(psms), str(path), **kwargs)
        return Path(path)

    return build


@pytest.fixture
def read_index():
    """Return {peptide: location} parsed from a specset _index.tsv."""

    def read(specset_dir):
        text = (Path(specset_dir) / "_index.tsv").read_text()
        return dict(line.split("\t", 1) for line in text.splitlines() if line)

    return read


@pytest.fixture
def read_location_index():
    """Return the raw lines of a location's __index.tsv."""

    def read(specset_dir, location):
        text = (Path(specset_dir) / location / "__index.tsv").read_text()
        return [line for line in text.splitlines() if line]

    return read


@pytest.fixture
def read_psm_fields():
    """Return {key: value} for a single PSM .tsv file."""

    def read(psm_path):
        fields = {}
        for line in Path(psm_path).read_text().splitlines():
            key, _, value = line.partition("\t")
            fields[key] = value
        return fields

    return read
