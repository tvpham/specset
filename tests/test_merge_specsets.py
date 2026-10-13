import pytest

from specset.utils import merge_specsets, merge_specsets_inplace


# Disjoint peptides from both specsets land in one output index, renumbered from 0.
def test_merge_disjoint_peptides(tmp_path, specset_builder, psm_factory, read_index, read_location_index):
    a = specset_builder(tmp_path / "a", [psm_factory("PEPA")])
    b = specset_builder(tmp_path / "b", [psm_factory("PEPB", z=3, matches="b2", intensities="10")])
    out = tmp_path / "out"

    merge_specsets(str(a), str(b), str(out))

    assert list(read_index(out)) == ["PEPA", "PEPB"]
    assert read_index(out)["PEPA"] == "0-4999/0"
    assert read_index(out)["PEPB"] == "0-4999/1"
    assert read_location_index(out, read_index(out)["PEPA"]) == ["0\t2"]
    assert read_location_index(out, read_index(out)["PEPB"]) == ["1\t3"]


# mode="ignore" keeps specset_a's PSMs and adds only the peptides new to it.
def test_merge_ignore_skips_overlapping_peptide(tmp_path, specset_builder, psm_factory, read_index, read_location_index):
    a = specset_builder(tmp_path / "a", [psm_factory("PEPTIDE", matches="b2;y5", intensities="100;200")])
    b = specset_builder(
        tmp_path / "b",
        [psm_factory("PEPTIDE", matches="b2", intensities="9"), psm_factory("PEPADK", matches="b3", intensities="7")],
    )
    out = tmp_path / "out"

    merge_specsets(str(a), str(b), str(out), mode="ignore")

    assert list(read_index(out)) == ["PEPTIDE", "PEPADK"]
    assert read_location_index(out, read_index(out)["PEPTIDE"]) == ["0\t2"]
    assert read_location_index(out, read_index(out)["PEPADK"]) == ["1\t2"]


# mode="both" keeps the overlapping peptide's PSMs from both specsets under one location.
def test_merge_both_keeps_overlapping_psms(tmp_path, specset_builder, psm_factory, read_index, read_location_index, read_psm_fields):
    a = specset_builder(tmp_path / "a", [psm_factory("PEPTIDE", matches="b2;y5", intensities="100;200")])
    b = specset_builder(tmp_path / "b", [psm_factory("PEPTIDE", matches="b2", intensities="9")])
    out = tmp_path / "out"

    merge_specsets(str(a), str(b), str(out), mode="both")

    assert list(read_index(out)) == ["PEPTIDE"]
    location = read_index(out)["PEPTIDE"]
    assert read_location_index(out, location) == ["0\t2", "1\t2"]
    assert read_psm_fields(out / location / "0.tsv")["ms2"] == "noloss,2,b,1,100.0;noloss,5,y,1,200.0;"
    assert read_psm_fields(out / location / "1.tsv")["ms2"] == "noloss,2,b,1,9.0;"


# An unknown mode is rejected before anything is written.
def test_merge_invalid_mode(tmp_path, specset_builder, psm_factory):
    a = specset_builder(tmp_path / "a", [psm_factory("PEPA")])
    b = specset_builder(tmp_path / "b", [psm_factory("PEPB")])

    with pytest.raises(ValueError):
        merge_specsets(str(a), str(b), str(tmp_path / "out"), mode="nope")


# The output must not be an input or a directory containing an input.
def test_merge_output_must_not_overlap_inputs(tmp_path, specset_builder, psm_factory):
    a = specset_builder(tmp_path / "a", [psm_factory("PEPA")])
    b = specset_builder(tmp_path / "b", [psm_factory("PEPB")])

    with pytest.raises(ValueError):
        merge_specsets(str(a), str(b), str(a))

    with pytest.raises(ValueError):
        merge_specsets(str(a), str(b), str(tmp_path))


# New peptides are appended to specset_a while the existing ones are preserved.
def test_merge_inplace_adds_new_peptides(tmp_path, specset_builder, psm_factory, read_index, read_location_index):
    a = specset_builder(tmp_path / "a", [psm_factory("PEPA")])
    b = specset_builder(tmp_path / "b", [psm_factory("PEPB", z=3, matches="b2", intensities="10")])

    merge_specsets_inplace(str(a), str(b))

    assert list(read_index(a)) == ["PEPA", "PEPB"]
    assert read_index(a)["PEPA"] == "0-4999/0"
    assert read_index(a)["PEPB"] == "0-4999/1"
    assert read_location_index(a, read_index(a)["PEPA"]) == ["0\t2"]
    assert read_location_index(a, read_index(a)["PEPB"]) == ["1\t3"]


# mode="ignore" adds only the peptides that are new to specset_a.
def test_merge_inplace_ignore_skips_overlapping_peptide(tmp_path, specset_builder, psm_factory, read_index, read_location_index):
    a = specset_builder(tmp_path / "a", [psm_factory("PEPTIDE")])
    b = specset_builder(
        tmp_path / "b",
        [psm_factory("PEPTIDE", matches="b2", intensities="9"), psm_factory("PEPADK", matches="b3", intensities="7")],
    )

    merge_specsets_inplace(str(a), str(b), mode="ignore")

    assert list(read_index(a)) == ["PEPTIDE", "PEPADK"]
    assert read_location_index(a, read_index(a)["PEPADK"]) == ["1\t2"]


# mode="both" appends the overlapping peptide's PSM to the existing location.
def test_merge_inplace_both_keeps_overlapping_psms(tmp_path, specset_builder, psm_factory, read_index, read_location_index):
    a = specset_builder(tmp_path / "a", [psm_factory("PEPTIDE", matches="b2;y5", intensities="100;200")])
    b = specset_builder(tmp_path / "b", [psm_factory("PEPTIDE", matches="b2", intensities="9")])

    merge_specsets_inplace(str(a), str(b), mode="both")

    assert list(read_index(a)) == ["PEPTIDE"]
    assert read_location_index(a, read_index(a)["PEPTIDE"]) == ["0\t2", "1\t2"]


def test_merge_inplace_missing_specset_a_index(tmp_path, specset_builder, psm_factory):
    a = tmp_path / "a"
    a.mkdir()
    b = specset_builder(tmp_path / "b", [psm_factory("PEPB")])

    with pytest.raises(FileNotFoundError):
        merge_specsets_inplace(str(a), str(b))


def test_merge_inplace_missing_specset_b_index(tmp_path, specset_builder, psm_factory):
    a = specset_builder(tmp_path / "a", [psm_factory("PEPA")])
    b = tmp_path / "b"
    b.mkdir()

    with pytest.raises(FileNotFoundError):
        merge_specsets_inplace(str(a), str(b))


def test_merge_inplace_missing_location_index(tmp_path, specset_builder, psm_factory):
    a = tmp_path / "a"
    a.mkdir()
    (a / "_index.tsv").write_text("PEPA\t0-4999/0\n")
    b = specset_builder(tmp_path / "b", [psm_factory("PEPB")])

    with pytest.raises(FileNotFoundError):
        merge_specsets_inplace(str(a), str(b))
