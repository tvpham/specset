import csv

import pytest

from specset.utils import _psm_stream_to_specset, specset_to_alphapeptdeep_string

EXPECTED_HEADER = [
    "FragmentLossType",
    "FragmentNumber",
    "FragmentType",
    "FragmentCharge",
    "RelativeIntensity",
    "ModifiedPeptide",
    "PrecursorCharge",
    "RT",
    "IonMobility",
]


def build_specset(tmp_path, psms, **kwargs):
    """Write a specset directory from a list of PSM dicts using the writer."""
    specset_dir = tmp_path / "specset"
    _psm_stream_to_specset(iter(psms), str(specset_dir), **kwargs)
    return specset_dir


def build_specset_raw(tmp_path, psms):
    """Write a minimal specset directory directly (bypasses the writer for edge cases)."""
    specset_dir = tmp_path / "specset"
    loc = specset_dir / "0-4999" / "0"
    loc.mkdir(parents=True)
    index2 = []
    pep = None
    for i, psm in enumerate(psms):
        with open(loc / f"{i}.tsv", "w", newline="") as f:
            for key, value in psm.items():
                f.write(f"{key}\t{value}\n")
        index2.append(f"{i}\t{psm['z']}\n")
        pep = psm["pep"]
    (loc / "__index.tsv").write_text("".join(index2))
    (specset_dir / "_index.tsv").write_text(f"{pep}\t0-4999/0\n")
    return specset_dir


def read_buffer(buffer):
    """Return (fieldnames, rows) parsed from a yielded StringIO buffer."""
    buffer.seek(0)
    reader = csv.DictReader(buffer, delimiter="\t")
    return reader.fieldnames, list(reader)


def make_psm(pep, z=2, rt=1.0, im=0.1, matches="b2;y5", intensities="100;200"):
    return {
        "pep": pep,
        "z": z,
        "rt": rt,
        "im": im,
        "matches": matches,
        "intensities": intensities,
    }


# One buffer is yielded per PSM file, regardless of the number of peptides.
def test_yields_one_buffer_per_psm(tmp_path):
    specset_dir = build_specset(
        tmp_path,
        [make_psm("PEPADA", matches="b2"), make_psm("PEPADK", matches="b2")],
    )

    buffers = list(specset_to_alphapeptdeep_string(str(specset_dir)))

    assert len(buffers) == 2


# The header of each yielded buffer lists the AlphaPeptDeep columns in the expected order.
def test_header_and_column_order(tmp_path):
    specset_dir = build_specset(tmp_path, [make_psm("PEPTIDE")])

    (buffer,) = list(specset_to_alphapeptdeep_string(str(specset_dir)))

    header, _ = read_buffer(buffer)
    assert header == EXPECTED_HEADER


# Fragment and metadata values in the buffer match the PSM and its ms2 fragments.
def test_fragment_and_metadata_values(tmp_path):
    specset_dir = build_specset(tmp_path, [make_psm("PEPTIDE")])

    (buffer,) = list(specset_to_alphapeptdeep_string(str(specset_dir)))

    _, rows = read_buffer(buffer)
    assert len(rows) == 2

    first = rows[0]
    assert first["FragmentLossType"] == "noloss"
    assert first["FragmentNumber"] == "2"
    assert first["FragmentType"] == "b"
    assert first["FragmentCharge"] == "1"
    assert first["RelativeIntensity"] == "100.0"
    assert first["ModifiedPeptide"] == "PEPTIDE"
    assert first["PrecursorCharge"] == "2"
    assert first["RT"] == "1.0"
    assert first["IonMobility"] == "0.1"

    assert rows[1]["FragmentType"] == "y"
    assert rows[1]["FragmentNumber"] == "5"


# The peptide is copied verbatim: no UniMod-to-square-bracket conversion or _..._ wrapping.
@pytest.mark.parametrize(
    "pep",
    ["PEPTIDE", "(UniMod:1)PEPTIDE", "PEPM(UniMod:35)TIDE", "_PEPTIDE_"],
)
def test_modified_peptide_written_as_is(tmp_path, pep):
    specset_dir = build_specset(tmp_path, [make_psm(pep, matches="b2", intensities="10")])

    (buffer,) = list(specset_to_alphapeptdeep_string(str(specset_dir)))

    _, rows = read_buffer(buffer)
    assert rows[0]["ModifiedPeptide"] == pep


# PSMs of the same peptide are yielded in __index.tsv order.
def test_multiple_psms_same_peptide_in_psm_order(tmp_path):
    specset_dir = build_specset(
        tmp_path,
        [
            make_psm("PEPTIDE", z=2, matches="b2", intensities="100"),
            make_psm("PEPTIDE", z=3, matches="b2", intensities="300"),
        ],
    )

    buffers = list(specset_to_alphapeptdeep_string(str(specset_dir)))

    assert [read_buffer(b)[1][0]["PrecursorCharge"] for b in buffers] == ["2", "3"]


# n limits the number of peptides processed, not the number of PSMs within one peptide.
def test_n_limits_number_of_peptides(tmp_path):
    specset_dir = build_specset(
        tmp_path,
        [make_psm(p, matches="b2", intensities="10") for p in ["PEPADA", "PEPADK", "PEPADR"]],
    )

    buffers = list(specset_to_alphapeptdeep_string(str(specset_dir), n=2))

    modified = [read_buffer(b)[1][0]["ModifiedPeptide"] for b in buffers]
    assert modified == ["PEPADA", "PEPADK"]
    assert "PEPADR" not in modified


# Missing rt/im fields fall back to 0 in the yielded buffer.
def test_missing_rt_and_im_fall_back_to_zero(tmp_path):
    specset_dir = build_specset_raw(tmp_path, [{"pep": "PEPTIDE", "z": "2", "ms2": "noloss,2,b,1,100.0;"}])

    (buffer,) = list(specset_to_alphapeptdeep_string(str(specset_dir)))

    header, rows = read_buffer(buffer)
    assert header == EXPECTED_HEADER
    assert rows[0]["RT"] == "0.0" or rows[0]["RT"] == "0"
    assert rows[0]["IonMobility"] == "0.0" or rows[0]["IonMobility"] == "0"


# A PSM without an ms2 field still yields a buffer, containing only the header.
def test_missing_ms2_yields_header_only_buffer(tmp_path):
    specset_dir = build_specset_raw(tmp_path, [{"pep": "PEPTIDE", "z": "2"}])

    (buffer,) = list(specset_to_alphapeptdeep_string(str(specset_dir)))

    header, rows = read_buffer(buffer)
    assert header == EXPECTED_HEADER
    assert rows == []