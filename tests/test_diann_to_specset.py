import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from specset.utils import diannlib_to_specset, diannreport_to_specset

LIB_COLUMNS = [
    "Precursor.Id",
    "Modified.Sequence",
    "Precursor.Charge",
    "RT",
    "IM",
    "Fragment.Series.Number",
    "Fragment.Type",
    "Fragment.Charge",
    "Fragment.Loss.Type",
    "Relative.Intensity",
    "Protein.Ids",
    "Genes",
]

REPORT_COLUMNS = [
    "Run",
    "Modified.Sequence",
    "Precursor.Charge",
    "RT",
    "IM",
    "Q.Value",
    "Decoy",
    "Ms2.Scan",
    "Protein.Ids",
    "Genes",
]


def write_parquet(path, columns, rows):
    table = pa.table({col: [r[i] for r in rows] for i, col in enumerate(columns)})
    pq.write_table(table, path)
    return path


def lib_row(pid, seq, z, num, ftype, fcharge=1, loss="noloss", intensity=100.0,
            rt=1.0, im=0.5, protein="P1", gene="G1"):
    return [pid, seq, z, rt, im, num, ftype, fcharge, loss, intensity, protein, gene]


def report_columns(n_frags=2):
    return REPORT_COLUMNS + [f"Fr.{i}.{field}" for i in range(n_frags) for field in ("Id", "Quantity")]


def report_row(run="run1", seq="PEPTIDE", z=2, rt=1.0, im=0.8, q=0.01, decoy=False, scan=10,
               protein="P1", gene="G1", frags=(("y8-unknown^1/816.4", 100.0), ("b2-noloss^1/200.5", 50.0))):
    return [run, seq, z, rt, im, q, decoy, scan, protein, gene] + [v for frag in frags for v in frag]


def parse_ms2(ms2):
    return [r.split(",") for r in ms2.rstrip(";").split(";") if r]


# Consecutive rows sharing a Precursor.Id become one PSM whose ms2 holds all fragments.
def test_library_basic(tmp_path, read_index, read_psm_fields):
    path = write_parquet(
        tmp_path / "lib.parquet",
        LIB_COLUMNS,
        [
            lib_row("1", "PEPTIDE", 2, 2, "b", intensity=100.0, protein="P1", gene="G1"),
            lib_row("1", "PEPTIDE", 2, 5, "y", intensity=200.0, protein="P1", gene="G1"),
            lib_row("2", "PEPADK", 3, 3, "b", intensity=50.0, rt=2.0, im=0.6),
        ],
    )
    out = tmp_path / "out"

    diannlib_to_specset(str(path), str(out))

    assert list(read_index(out)) == ["PEPTIDE", "PEPADK"]
    fields = read_psm_fields(out / read_index(out)["PEPTIDE"] / "0.tsv")
    assert fields["pep"] == "PEPTIDE"
    assert fields["z"] == "2"
    assert fields["rt"] == "1.0"
    assert fields["im"] == "0.5"
    assert fields["proteins"] == "P1"
    assert fields["genes"] == "G1"
    assert parse_ms2(fields["ms2"]) == [
        ["noloss", "2", "b", "1", "100.0"],
        ["noloss", "5", "y", "1", "200.0"],
    ]


# Library loss names and 2+ fragment charges map to the MaxQuant notation.
def test_library_loss_and_charge_conversion(tmp_path, read_index, read_psm_fields):
    path = write_parquet(
        tmp_path / "lib.parquet",
        LIB_COLUMNS,
        [
            lib_row("1", "PEPTIDE", 2, 2, "b", fcharge=2, loss="noloss"),
            lib_row("1", "PEPTIDE", 2, 3, "y", loss="H3PO4"),
            lib_row("1", "PEPTIDE", 2, 4, "b", loss="H2O"),
            lib_row("1", "PEPTIDE", 2, 5, "y", loss="NH3"),
        ],
    )
    out = tmp_path / "out"

    diannlib_to_specset(str(path), str(out))

    fields = read_psm_fields(out / read_index(out)["PEPTIDE"] / "0.tsv")
    assert parse_ms2(fields["ms2"]) == [
        ["noloss", "2", "b", "2", "100.0"],
        ["H3PO4", "3", "y", "1", "100.0"],
        ["H2O", "4", "b", "1", "100.0"],
        ["NH3", "5", "y", "1", "100.0"],
    ]


# A report row becomes one PSM carrying run/scan/protein metadata and its fragments.
def test_report_basic(tmp_path, read_index, read_psm_fields):
    path = write_parquet(tmp_path / "report.parquet", report_columns(2), [report_row()])
    out = tmp_path / "out"

    diannreport_to_specset(str(path), str(out))

    assert read_index(out) == {"PEPTIDE": "0-4999/0"}
    fields = read_psm_fields(out / read_index(out)["PEPTIDE"] / "0.tsv")
    assert fields["pep"] == "PEPTIDE"
    assert fields["z"] == "2"
    assert fields["rt"] == "1.0"
    assert fields["im"] == "0.8"
    assert fields["run"] == "run1"
    assert fields["scan"] == "10"
    assert fields["qvalue"] == "0.01"
    assert fields["proteins"] == "P1"
    assert fields["genes"] == "G1"
    assert parse_ms2(fields["ms2"]) == [
        ["noloss", "8", "y", "1", "100.0"],
        ["noloss", "2", "b", "1", "50.0"],
    ]


# Decoy rows are skipped.
def test_report_decoy_rows_skipped(tmp_path, read_index):
    rows = [report_row(seq="DECOY", decoy=True), report_row(seq="PEPADK")]
    path = write_parquet(tmp_path / "report.parquet", report_columns(2), rows)
    out = tmp_path / "out"

    diannreport_to_specset(str(path), str(out))

    assert read_index(out) == {"PEPADK": "0-4999/0"}


# string_equal keeps only the selected rows.
def test_report_column_filter(tmp_path, read_index):
    rows = [report_row(run="run1"), report_row(run="run2", seq="PEPADK")]
    path = write_parquet(tmp_path / "report.parquet", report_columns(2), rows)
    out = tmp_path / "out"

    diannreport_to_specset(str(path), str(out), string_equal={"Run": ["run1"]})

    assert read_index(out) == {"PEPTIDE": "0-4999/0"}


# Rows with no usable fragment (empty id, zero quantity) are skipped.
def test_report_rows_without_fragments_skipped(tmp_path, read_index):
    rows = [report_row(frags=(("", 0.0), ("", 0.0)))]
    path = write_parquet(tmp_path / "report.parquet", report_columns(2), rows)
    out = tmp_path / "out"

    diannreport_to_specset(str(path), str(out))

    assert read_index(out) == {}


# A fragment whose loss name is not recognized is dropped.
def test_report_unknown_loss_skipped(tmp_path, read_index, read_psm_fields):
    rows = [report_row(frags=(("y8-unknownloss^1/816.4", 100.0), ("b2-noloss^1/200.5", 50.0)))]
    path = write_parquet(tmp_path / "report.parquet", report_columns(2), rows)
    out = tmp_path / "out"

    diannreport_to_specset(str(path), str(out))

    fields = read_psm_fields(out / read_index(out)["PEPTIDE"] / "0.tsv")
    assert parse_ms2(fields["ms2"]) == [["noloss", "2", "b", "1", "50.0"]]


# Requesting a column that is not in the report raises.
def test_report_missing_column_raises(tmp_path):
    path = write_parquet(tmp_path / "report.parquet", report_columns(2), [report_row()])

    with pytest.raises(ValueError):
        diannreport_to_specset(str(path), str(tmp_path / "out"), string_equal={"NotAColumn": ["x"]})
