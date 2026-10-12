import csv
from pathlib import Path

import pytest

from specset.utils import _psm_stream_to_specset, specset_to_speclib

EXPECTED_HEADER = [
    "ModifiedPeptide",
    "PrecursorCharge",
    "RT",
    "IM",
    "FragmentType",
    "RelativeIntensity",
    "FragmentCharge",
    "FragmentNumber",
    "FragmentLossType",
    "PrecursorMz",
    "FragmentMz",
]


def build_specset(tmp_path, psms, **kwargs):
    """Write a specset directory from a list of PSM dicts using the writer."""
    specset_dir = tmp_path / "specset"
    _psm_stream_to_specset(iter(psms), str(specset_dir), **kwargs)
    return specset_dir


def read_library(lib_path):
    """Return (fieldnames, rows) parsed from the speclib TSV."""
    with open(lib_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return reader.fieldnames, list(reader)


def make_psm(pep, z=2, rt=1.0, im=0.1, matches="b2;y5", intensities="100;200", **extra):
    psm = {"pep": pep, "z": z, "rt": rt, "im": im, "matches": matches, "intensities": intensities}
    psm.update(extra)
    return psm


# One PSM: checks the exact header, the two fragment rows, and every column of the first row.
def test_header_and_single_psm(tmp_path):
    specset_dir = build_specset(tmp_path, [make_psm("PEPTIDE")])
    lib = tmp_path / "lib.tsv"

    specset_to_speclib(str(specset_dir), str(lib))

    header, rows = read_library(lib)
    assert header == EXPECTED_HEADER
    assert len(rows) == 2

    first, second = rows
    assert first["ModifiedPeptide"] == "_PEPTIDE_"
    assert first["PrecursorCharge"] == "2"
    assert first["RT"] == "1.0"
    assert first["IM"] == "0.1"
    assert first["FragmentType"] == "b"
    assert first["RelativeIntensity"] == "100.0"
    assert first["FragmentCharge"] == "1"
    assert first["FragmentNumber"] == "2"
    assert first["FragmentLossType"] == "noloss"
    assert second["FragmentType"] == "y"
    assert second["FragmentNumber"] == "5"


# Fixed reference m/z values for PEPTIDE (charge 2): b2/y5, single charge, no loss.
def test_precursor_and_fragment_mz_values(tmp_path):
    specset_dir = build_specset(tmp_path, [make_psm("PEPTIDE")])
    lib = tmp_path / "lib.tsv"

    specset_to_speclib(str(specset_dir), str(lib))

    _, rows = read_library(lib)
    assert {row["PrecursorMz"] for row in rows} == {"400.68725848012497"}
    assert rows[0]["FragmentMz"] == "227.10263340359"
    assert rows[1]["FragmentMz"] == "574.27188355666"


# The same peptide at two charges produces one group per charge (_read_location groups by z).
def test_multiple_charges_are_separate_groups(tmp_path):
    specset_dir = build_specset(
        tmp_path,
        [
            make_psm("PEPTIDE", z=2, rt=1.0, im=0.1, matches="b2", intensities="100"),
            make_psm("PEPTIDE", z=3, rt=2.0, im=0.2, matches="b2", intensities="300"),
        ],
    )
    lib = tmp_path / "lib.tsv"

    specset_to_speclib(str(specset_dir), str(lib))

    _, rows = read_library(lib)
    assert len(rows) == 2
    assert {row["PrecursorCharge"] for row in rows} == {"2", "3"}


# Two distinct peptides must both appear in the library.
def test_multiple_peptides(tmp_path):
    specset_dir = build_specset(
        tmp_path,
        [make_psm("PEPADA", matches="b2", intensities="10"), make_psm("PEPADK", matches="b2", intensities="20")],
    )
    lib = tmp_path / "lib.tsv"

    specset_to_speclib(str(specset_dir), str(lib))

    _, rows = read_library(lib)
    assert {row["ModifiedPeptide"] for row in rows} == {"_PEPADA_", "_PEPADK_"}


# n=2 processes only the first two peptides, in _index.tsv order.
def test_n_limits_number_of_peptides(tmp_path):
    peps = ["PEPADA", "PEPADK", "PEPADR"]
    specset_dir = build_specset(
        tmp_path, [make_psm(p, matches="b2", intensities="10") for p in peps]
    )
    lib = tmp_path / "lib.tsv"

    specset_to_speclib(str(specset_dir), str(lib), n=2)

    _, rows = read_library(lib)
    assert {row["ModifiedPeptide"] for row in rows} == {"_PEPADA_", "_PEPADK_"}


# 'X' is not in aa_mass, so the peptide is skipped and only the header is written.
def test_unknown_amino_acid_is_skipped(tmp_path):
    specset_dir = build_specset(tmp_path, [make_psm("PEPXIDE", matches="b2;y2", intensities="10;20")])
    lib = tmp_path / "lib.tsv"

    specset_to_speclib(str(specset_dir), str(lib))

    header, rows = read_library(lib)
    assert header == EXPECTED_HEADER
    assert rows == []


# additional_columns appends the column to the header and carries the PSM value through.
def test_additional_columns_are_written(tmp_path):
    specset_dir = build_specset(
        tmp_path, [make_psm("PEPTIDE", matches="b2", intensities="100", genes="GENE1")]
    )
    lib = tmp_path / "lib.tsv"

    specset_to_speclib(str(specset_dir), str(lib), additional_columns={"genes": "genes"})

    header, rows = read_library(lib)
    assert header == EXPECTED_HEADER + ["genes"]
    assert rows[0]["genes"] == "GENE1"


# UniMod 35/21/1/4 round-trip to the _PEP[UniMod:n]..._ output form.
@pytest.mark.parametrize(
    "pep, expected",
    [
        ("PEPM(UniMod:35)TIDE", "_PEPM[UniMod:35]TIDE_"),
        ("PEPT(UniMod:21)DE", "_PEPT[UniMod:21]DE_"),
        ("(UniMod:1)PEPTIDE", "_[UniMod:1]PEPTIDE_"),
        ("PEPC(UniMod:4)TIDE", "_PEPC[UniMod:4]TIDE_"),
    ],
)
def test_modified_peptide_notation_conversion(tmp_path, pep, expected):
    specset_dir = build_specset(tmp_path, [make_psm(pep, matches="b2;y2", intensities="10;20")])
    lib = tmp_path / "lib.tsv"

    specset_to_speclib(str(specset_dir), str(lib))

    _, rows = read_library(lib)
    assert rows
    assert {row["ModifiedPeptide"] for row in rows} == {expected}


# consensus_method="score" picks the higher-score PSM even though it has fewer fragments
# (the column must also be listed in additional_columns to be read).
def test_consensus_by_column_selects_matching_psm(tmp_path):
    specset_dir = build_specset(
        tmp_path,
        [
            make_psm("PEPTIDE", matches="b2", intensities="100", score=10),
            make_psm("PEPTIDE", rt=2.0, matches="b2;y5", intensities="300;400", score=5),
        ],
    )
    lib = tmp_path / "lib.tsv"

    specset_to_speclib(
        str(specset_dir),
        str(lib),
        consensus_method="score",
        additional_columns={"score": "score"},
    )

    _, rows = read_library(lib)
    assert len(rows) == 1
    assert rows[0]["FragmentNumber"] == "2"
    assert rows[0]["RelativeIntensity"] == "100.0"
    assert rows[0]["RT"] == "1.0"
    assert rows[0]["score"] == "10"


# With no consensus method the PSM with more fragments wins; RT/intensities come from that PSM.
def test_consensus_default_selects_psm_with_most_fragments(tmp_path):
    specset_dir = build_specset(
        tmp_path,
        [
            make_psm("PEPTIDE", matches="b2", intensities="100", score=10),
            make_psm("PEPTIDE", rt=2.0, matches="b2;y5", intensities="300;400", score=5),
        ],
    )
    lib = tmp_path / "lib.tsv"

    specset_to_speclib(str(specset_dir), str(lib), additional_columns={"score": "score"})

    _, rows = read_library(lib)
    assert len(rows) == 2
    assert {row["RelativeIntensity"] for row in rows} == {"300.0", "400.0"}
    assert {row["RT"] for row in rows} == {"2.0"}


# Expects per-fragment means (b2 -> 200, y5 -> 300). Currently xfail: _consensus merges the
# per-fragment averages against a single metadata row, so all but one fragment get NaN metadata
# and specset_to_speclib crashes.
@pytest.mark.xfail(
    reason="mean/median consensus merges per-fragment averages with a single metadata row, "
    "producing NaN metadata for all but one fragment (bug in _consensus)",
)
def test_consensus_mean_averages_intensities(tmp_path):
    specset_dir = build_specset(
        tmp_path,
        [
            make_psm("PEPTIDE", rt=1.0, matches="b2;y5", intensities="100;200"),
            make_psm("PEPTIDE", rt=2.0, matches="b2;y5", intensities="300;400"),
        ],
    )
    lib = tmp_path / "lib.tsv"

    specset_to_speclib(str(specset_dir), str(lib), consensus_method="mean")

    _, rows = read_library(lib)
    intensities = {row["FragmentNumber"]: row["RelativeIntensity"] for row in rows}
    assert intensities == {"2": "200.0", "5": "300.0"}


# Expects per-fragment medians (b2 -> 300, y5 -> 400). Same xfail reason as the mean test.
@pytest.mark.xfail(
    reason="mean/median consensus merges per-fragment averages with a single metadata row, "
    "producing NaN metadata for all but one fragment (bug in _consensus)",
)
def test_consensus_median_intensities(tmp_path):
    specset_dir = build_specset(
        tmp_path,
        [
            make_psm("PEPTIDE", rt=1.0, matches="b2;y5", intensities="100;200"),
            make_psm("PEPTIDE", rt=2.0, matches="b2;y5", intensities="300;400"),
            make_psm("PEPTIDE", rt=3.0, matches="b2;y5", intensities="500;600"),
        ],
    )
    lib = tmp_path / "lib.tsv"

    specset_to_speclib(str(specset_dir), str(lib), consensus_method="median")

    _, rows = read_library(lib)
    intensities = {row["FragmentNumber"]: row["RelativeIntensity"] for row in rows}
    assert intensities == {"2": "300.0", "5": "400.0"}
