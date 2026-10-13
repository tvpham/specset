import pytest

from specset.utils import matches_intensities_to_psm


# MaxQuant loss suffixes map to their named FragmentLossType.
def test_loss_types_are_mapped():
    df = matches_intensities_to_psm("b2;y3*;b4-H2O;y5-NH3;b6-CO;y7-H4COS", "1;2;3;4;5;6")
    assert df["FragmentLossType"].tolist() == ["noloss", "H3PO4", "H2O", "NH3", "CO", "H4COS"]


# "(n+)" fragment charges are parsed; fragments without a suffix default to charge 1.
def test_fragment_charges_are_parsed():
    df = matches_intensities_to_psm("b2(2+);y3(3+);a1", "1;2;3")
    assert df["FragmentCharge"].tolist() == ["2", "3", "1"]


# Loss and charge suffixes can be combined on the same fragment.
def test_loss_and_charge_combined():
    df = matches_intensities_to_psm("y3*(2+);b4-H2O(3+)", "1;2")
    assert df["FragmentLossType"].tolist() == ["H3PO4", "H2O"]
    assert df["FragmentCharge"].tolist() == ["2", "3"]


def test_fragment_type_and_number_are_parsed():
    df = matches_intensities_to_psm("b2;y10;a1", "1;2;3")
    assert df[["FragmentType", "FragmentNumber"]].values.tolist() == [["b", 2], ["y", 10], ["a", 1]]


# A fragment from the "bad" set aborts the whole PSM.
def test_bad_fragment_returns_none():
    assert matches_intensities_to_psm("b2;acK", "1;2") is None
    assert matches_intensities_to_psm("", "1") is None


# pY and a2 are diagnostic/skip fragments and are dropped, not fatal.
def test_skip_fragments_are_dropped():
    df = matches_intensities_to_psm("b2;pY;y3;a2", "1;2;3;4")
    assert df["FragmentType"].tolist() == ["b", "y"]


# A fragment that does not match the pattern is skipped.
def test_unparseable_fragment_is_skipped():
    df = matches_intensities_to_psm("b2;??;y3", "1;2;3")
    assert df["FragmentType"].tolist() == ["b", "y"]


# If no fragment can be parsed, the PSM is rejected.
def test_no_parseable_fragments_returns_none():
    assert matches_intensities_to_psm("??;!!", "1;2") is None


# More fragments than intensities: the extra fragments are dropped.
def test_mismatched_lengths_are_truncated():
    df = matches_intensities_to_psm("b2;y3", "1")
    assert df["FragmentType"].tolist() == ["b"]


# Non-numeric or NaN intensities abort the whole PSM.
def test_invalid_intensity_returns_none():
    assert matches_intensities_to_psm("b2;y3", "1;oops") is None
    assert matches_intensities_to_psm("b2;y3", "1;nan") is None


# max_intensity > 0 rescales relative intensities so the maximum equals max_intensity.
def test_max_intensity_normalizes():
    df = matches_intensities_to_psm("b2;y3", "10;20", max_intensity=100.0)
    assert df["RelativeIntensity"].tolist() == [50.0, 100.0]


# max_intensity <= 0 leaves intensities untouched.
@pytest.mark.parametrize("max_intensity", [-1.0, 0.0])
def test_max_intensity_disabled_leaves_intensities(max_intensity):
    df = matches_intensities_to_psm("b2;y3", "10;20", max_intensity=max_intensity)
    assert df["RelativeIntensity"].tolist() == [10.0, 20.0]


# The writer threads max_intensity through to the stored ms2 field.
def test_writer_applies_max_intensity(tmp_path, specset_builder, psm_factory, read_index, read_psm_fields):
    specset = specset_builder(
        tmp_path / "s",
        [psm_factory("PEPTIDE", matches="b2;y5", intensities="10;20")],
        max_intensity=50.0,
    )

    location = read_index(specset)["PEPTIDE"]
    fields = read_psm_fields(specset / location / "0.tsv")
    assert fields["ms2"] == "noloss,2,b,1,25.0;noloss,5,y,1,50.0;"
