from specset.utils import easypqp_to_specset

COLUMNS = [
    "ModifiedPeptideSequence",
    "PrecursorCharge",
    "FragmentLossType",
    "FragmentSeriesNumber",
    "FragmentType",
    "FragmentCharge",
    "LibraryIntensity",
    "NormalizedRetentionTime",
    "PrecursorIonMobility",
    "ProteinId",
    "GeneName",
]


def write_easypqp(path, rows):
    lines = ["\t".join(COLUMNS)] + ["\t".join(str(v) for v in row) for row in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def row(seq="PEPTIDE", z=2, loss="noloss", num=2, ftype="b", fcharge=1, intensity=100.0,
        rt=1.5, im=0.8, protein="PROT1", gene="GENE1"):
    return [seq, z, loss, num, ftype, fcharge, intensity, rt, im, protein, gene]


def parse_ms2(ms2):
    return [r.split(",") for r in ms2.rstrip(";").split(";") if r]


# One peptide: fragment fields, PSM metadata and the extra ProteinId/GeneName columns.
def test_basic_single_peptide(tmp_path, read_index, read_psm_fields):
    path = write_easypqp(tmp_path / "easypqp.tsv", [row(), row(num=5, ftype="y", intensity=200.0)])
    out = tmp_path / "out"

    easypqp_to_specset(str(path), str(out))

    assert read_index(out) == {"PEPTIDE": "0-4999/0"}
    fields = read_psm_fields(out / read_index(out)["PEPTIDE"] / "0.tsv")
    assert fields["pep"] == "PEPTIDE"
    assert fields["z"] == "2"
    assert fields["rt"] == "1.5"
    assert fields["im"] == "0.8"
    assert fields["ProteinId"] == "PROT1"
    assert fields["GeneName"] == "GENE1"
    assert parse_ms2(fields["ms2"]) == [
        ["noloss", "2", "b", "1", "100.0"],
        ["noloss", "5", "y", "1", "200.0"],
    ]


# EasyPQP loss names and 2+ fragment charges are mapped to the MaxQuant notation.
def test_loss_and_charge_conversion(tmp_path, read_index, read_psm_fields):
    rows = [
        row(num=2, ftype="b", fcharge=2),
        row(loss="H3PO4", num=3, ftype="y"),
        row(loss="H2O", num=4, ftype="b"),
        row(loss="NH3", num=5, ftype="y"),
    ]
    path = write_easypqp(tmp_path / "easypqp.tsv", rows)
    out = tmp_path / "out"

    easypqp_to_specset(str(path), str(out))

    fields = read_psm_fields(out / read_index(out)["PEPTIDE"] / "0.tsv")
    assert parse_ms2(fields["ms2"]) == [
        ["noloss", "2", "b", "2", "100.0"],
        ["H3PO4", "3", "y", "1", "100.0"],
        ["H2O", "4", "b", "1", "100.0"],
        ["NH3", "5", "y", "1", "100.0"],
    ]


# The same peptide at two charges produces two PSMs under one location.
def test_same_peptide_two_charges(tmp_path, read_index, read_location_index):
    path = write_easypqp(tmp_path / "easypqp.tsv", [row(z=2), row(z=3, num=5, ftype="y")])
    out = tmp_path / "out"

    easypqp_to_specset(str(path), str(out))

    assert read_index(out) == {"PEPTIDE": "0-4999/0"}
    assert read_location_index(out, read_index(out)["PEPTIDE"]) == ["0\t2", "1\t3"]


# Distinct peptides are grouped in sequence order.
def test_multiple_peptides(tmp_path, read_index):
    rows = [row(seq="PEPADA"), row(seq="PEPADK")]
    path = write_easypqp(tmp_path / "easypqp.tsv", rows)
    out = tmp_path / "out"

    easypqp_to_specset(str(path), str(out))

    assert list(read_index(out)) == ["PEPADA", "PEPADK"]
