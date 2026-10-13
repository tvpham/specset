from pathlib import Path

import pytest

from specset.utils import maxquant_to_specset

MSMS_HEADER = [
    "Modified sequence",
    "Charge",
    "Retention time",
    "Matches",
    "Intensities",
    "Evidence ID",
]


def write_maxquant_dir(tmp_path, msms_rows, evidence_rows, msms_extra_header=()):
    """Write a minimal MaxQuant msms.txt/evidence.txt pair into tmp_path/mq."""
    txt_dir = tmp_path / "mq"
    txt_dir.mkdir()

    header = MSMS_HEADER + list(msms_extra_header)
    lines = ["\t".join(header)] + ["\t".join(row) for row in msms_rows]
    (txt_dir / "msms.txt").write_text("\n".join(lines) + "\n")

    ev_lines = ["\t".join(["id", "1/K0"])] + ["\t".join(row) for row in evidence_rows]
    (txt_dir / "evidence.txt").write_text("\n".join(ev_lines) + "\n")

    return txt_dir


def read_index(out_dir):
    """Return {peptide: location} from the main _index.tsv."""
    text = (Path(out_dir) / "_index.tsv").read_text()
    return dict(line.split("\t", 1) for line in text.splitlines() if line)


def read_psm_fields(psm_path):
    """Return {key: value} for a single PSM .tsv file."""
    fields = {}
    for line in Path(psm_path).read_text().splitlines():
        key, _, value = line.partition("\t")
        fields[key] = value
    return fields


def parse_ms2(ms2):
    """Split an ms2 field into a list of fragment column lists."""
    return [row.split(",") for row in ms2.rstrip(";").split(";") if row]


# One PSM: checks _index.tsv, __index.tsv, the PSM fields, and the parsed ms2 fragments.
def test_basic_single_psm(tmp_path):
    txt_dir = write_maxquant_dir(
        tmp_path,
        [["_PEPTIDE_", "2", "10.5", "b2;y2", "100;200", "1"]],
        [["1", "0.85"]],
    )
    out_dir = tmp_path / "out"

    maxquant_to_specset(str(txt_dir), str(out_dir))

    assert read_index(out_dir) == {"PEPTIDE": "0-4999/0"}

    psm_dir = out_dir / "0-4999" / "0"
    assert (psm_dir / "__index.tsv").read_text() == "0\t2\n"

    fields = read_psm_fields(psm_dir / "0.tsv")
    assert fields["pep"] == "PEPTIDE"
    assert fields["z"] == "2"
    assert fields["rt"] == "10.5"
    assert fields["im"] == "0.85"
    assert parse_ms2(fields["ms2"]) == [
        ["noloss", "2", "b", "1", "100.0"],
        ["noloss", "2", "y", "1", "200.0"],
    ]


# MaxQuant modification notation is converted to UniMod notation in the peptide.
@pytest.mark.parametrize(
    "modified_sequence, expected",
    [
        ("_PEPTIDE_", "PEPTIDE"),
        ("_T(Phospho (STY))PEPTIDE_", "T(UniMod:21)PEPTIDE"),
        ("_(Acetyl (Protein N-term))PEPTIDE_", "(UniMod:1)PEPTIDE"),
        ("_PEPM(Oxidation (M))TIDE_", "PEPM(UniMod:35)TIDE"),
        ("_PEPCTIDE_", "PEPC(UniMod:4)TIDE"),
    ],
)
def test_modified_sequence_conversion(tmp_path, modified_sequence, expected):
    txt_dir = write_maxquant_dir(
        tmp_path,
        [[modified_sequence, "2", "1.0", "b2;y2", "10;20", "1"]],
        [["1", "0.5"]],
    )
    out_dir = tmp_path / "out"

    maxquant_to_specset(str(txt_dir), str(out_dir))

    location = read_index(out_dir)[expected]
    fields = read_psm_fields(out_dir / location / "0.tsv")
    assert fields["pep"] == expected


# n=3 keeps only the first three PSMs in the output index.
def test_n_limits_number_of_psms(tmp_path):
    msms_rows = [
        [f"_PEP{i}_", "2", str(i), "b2;y2", "10;20", str(i)] for i in range(5)
    ]
    evidence_rows = [[str(i), "0.5"] for i in range(5)]
    txt_dir = write_maxquant_dir(tmp_path, msms_rows, evidence_rows)
    out_dir = tmp_path / "out"

    maxquant_to_specset(str(txt_dir), str(out_dir), n=3)

    assert set(read_index(out_dir)) == {"PEP0", "PEP1", "PEP2"}


# Multiple PSMs for the same peptide share one location dir; both entries are indexed.
def test_multiple_psms_same_peptide_share_location(tmp_path):
    txt_dir = write_maxquant_dir(
        tmp_path,
        [
            ["_PEPTIDE_", "2", "1.0", "b2;y2", "10;20", "1"],
            ["_PEPTIDE_", "3", "2.0", "b2;y2", "30;40", "2"],
        ],
        [["1", "0.5"], ["2", "0.6"]],
    )
    out_dir = tmp_path / "out"

    maxquant_to_specset(str(txt_dir), str(out_dir))

    assert read_index(out_dir) == {"PEPTIDE": "0-4999/0"}
    psm_dir = out_dir / "0-4999" / "0"
    assert (psm_dir / "__index.tsv").read_text() == "0\t2\n1\t3\n"

    second = read_psm_fields(psm_dir / "1.tsv")
    assert second["z"] == "3"
    assert second["rt"] == "2.0"
    assert second["im"] == "0.6"


# PSMs flagged as reverse hits (Reverse = "+") are filtered out of the output.
def test_reverse_hits_are_filtered(tmp_path):
    txt_dir = write_maxquant_dir(
        tmp_path,
        [
            ["_REVERSED_", "2", "1.0", "b2;y2", "10;20", "1", "+"],
            ["_PEPTIDEK_", "2", "2.0", "b2;y2", "10;20", "2", ""],
        ],
        [["1", "0.5"], ["2", "0.6"]],
        msms_extra_header=["Reverse"],
    )
    out_dir = tmp_path / "out"

    maxquant_to_specset(
        str(txt_dir), str(out_dir), additional_columns={"Reverse": "Reverse"}
    )

    assert set(read_index(out_dir)) == {"PEPTIDEK"}


# Additional columns passed via additional_columns are written into the PSM file.
def test_additional_columns_are_written(tmp_path):
    txt_dir = write_maxquant_dir(
        tmp_path,
        [["_PEPTIDE_", "2", "1.0", "b2;y2", "10;20", "1", "GENE1"]],
        [["1", "0.5"]],
        msms_extra_header=["Gene Names"],
    )
    out_dir = tmp_path / "out"

    maxquant_to_specset(
        str(txt_dir), str(out_dir), additional_columns={"Gene Names": "genes"}
    )

    location = read_index(out_dir)["PEPTIDE"]
    fields = read_psm_fields(out_dir / location / "0.tsv")
    assert fields["genes"] == "GENE1"


# PSMs with only invalid fragments are skipped, leaving an empty specset.
def test_invalid_fragments_are_skipped(tmp_path):
    txt_dir = write_maxquant_dir(
        tmp_path,
        [["_PEPTIDE_", "2", "1.0", "b2;acK", "10;20", "1"]],
        [["1", "0.5"]],
    )
    out_dir = tmp_path / "out"

    maxquant_to_specset(str(txt_dir), str(out_dir))

    assert read_index(out_dir) == {}
    assert (out_dir / "_index.tsv").read_text() == ""
    assert list((out_dir / "0-4999").iterdir()) == []


# New batch folders are created once max_files_in_directory is exceeded.
def test_batch_folders_are_created(tmp_path):
    msms_rows = [
        [f"_PEP{i}_", "2", str(i), "b2;y2", "10;20", str(i)] for i in range(3)
    ]
    evidence_rows = [[str(i), "0.5"] for i in range(3)]
    txt_dir = write_maxquant_dir(tmp_path, msms_rows, evidence_rows)
    out_dir = tmp_path / "out"

    maxquant_to_specset(str(txt_dir), str(out_dir), max_files_in_directory=2)

    assert (out_dir / "0-1").is_dir()
    assert (out_dir / "2-3").is_dir()
    index = read_index(out_dir)
    assert index["PEP0"].startswith("0-1/")
    assert index["PEP1"].startswith("0-1/")
    assert index["PEP2"].startswith("2-3/")


# Writing into a non-empty out_dir raises FileExistsError unless overwrite=True is given.
def test_overwrite_behaviour(tmp_path):
    txt_dir = write_maxquant_dir(
        tmp_path,
        [["_PEPTIDE_", "2", "1.0", "b2;y2", "10;20", "1"]],
        [["1", "0.5"]],
    )
    out_dir = tmp_path / "out"

    maxquant_to_specset(str(txt_dir), str(out_dir))

    with pytest.raises(FileExistsError):
        maxquant_to_specset(str(txt_dir), str(out_dir))

    maxquant_to_specset(str(txt_dir), str(out_dir), overwrite=True)
    assert read_index(out_dir) == {"PEPTIDE": "0-4999/0"}
