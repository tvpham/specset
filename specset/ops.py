"""Operations on specset directories."""

import io
import os
import json
import pandas as pd


def _read_psm_file(path):
    """Read a specset PSM .tsv file into a dict."""
    a = pd.read_csv(path, header=None, sep="\t", index_col=0)
    return {k: a.loc[k, 1] for k in a.index}


def _read_ms2(psm):
    """Parse the ms2 JSON string in a PSM dict into a DataFrame."""
    df = pd.read_json(io.StringIO(psm["ms2"]), orient="values")
    df.columns = ["FragmentLossType", "FragmentNumber", "FragmentType", "FragmentCharge", "RelativeIntensity"]
    df["FragmentNumber"] = df["FragmentNumber"].astype(int)
    df["RelativeIntensity"] = df["RelativeIntensity"].astype(float)
    return df


def _write_psm_file(path, psm):
    """Write a PSM dict to a specset .tsv file."""
    with open(path, "w", newline="") as f:
        for key, value in psm.items():
            f.write(f"{key}\t{value}\n")


def _aggregate_fragments(ms2_list, method):
    """Aggregate a list of fragment DataFrames by method: average, median, or max_score."""
    combined = pd.concat(ms2_list, ignore_index=True)
    group_cols = ["FragmentLossType", "FragmentNumber", "FragmentType", "FragmentCharge"]

    if method == "average":
        result = combined.groupby(group_cols)["RelativeIntensity"].mean().reset_index()
    elif method == "median":
        result = combined.groupby(group_cols)["RelativeIntensity"].median().reset_index()
    else:
        raise ValueError(f"Unknown method: {method!r}. Choose 'average', 'median', or 'max_score'.")

    max_inten = result["RelativeIntensity"].max()
    if max_inten > 0:
        result["RelativeIntensity"] = result["RelativeIntensity"] / max_inten

    return result


def consensus(in_dir, out_dir, method="max_score"):
    """Compute consensus spectra from a specset directory.

    Parameters
    ----------
    in_dir : str
        Path to the input specset directory.
    out_dir : str
        Path to write the output specset directory.
    method : str
        Aggregation method: 'average', 'median', or 'max_score'.
    """
    os.makedirs(out_dir, exist_ok=True)
    index_file = os.path.join(out_dir, "_index.tsv")

    if os.path.exists(index_file):
        os.remove(index_file)

    index_df = pd.read_csv(os.path.join(in_dir, "_index.tsv"),
                           header=None, names=["peptide", "location"], sep="\t")

    for i, row in index_df.iterrows():
        folder_name = os.path.dirname(row["location"])
        os.makedirs(os.path.join(out_dir, folder_name), exist_ok=True)

        precursor_dir = os.path.join(in_dir, row["location"])
        content_df = pd.read_csv(os.path.join(precursor_dir, "__index.tsv"),
                                 header=None, names=["file", "charge"], sep="\t")

        out_precursor_dir = os.path.join(out_dir, row["location"])
        os.makedirs(out_precursor_dir, exist_ok=True)

        # Write one entry per peptide to _index.tsv
        with open(index_file, "a") as f:
            f.write(f"{row['peptide']}\t{row['location']}\n")

        for charge, charge_group in content_df.groupby("charge"):
            psm_files = [_read_psm_file(os.path.join(precursor_dir, f"{f}.tsv"))
                         for f in charge_group["file"]]

            if method == "max_score":
                best = max(psm_files, key=lambda p: float(p["score"]))
                consensus_psm = dict(best)
            else:
                ms2_list = [_read_ms2(p) for p in psm_files]
                agg_df = _aggregate_fragments(ms2_list, method)
                best = max(psm_files, key=lambda p: float(p["score"]))
                consensus_psm = {
                    "ms2":   agg_df.to_json(orient="values"),
                    "ret":   float(best["ret"]),
                    "ccs":   float(best["ccs"]),
                    "pep":   best["pep"],
                    "z":     charge,
                    "score": float(best["score"]),
                    "id":    best["id"],
                }

            _write_psm_file(os.path.join(out_precursor_dir, f"{charge}.tsv"), consensus_psm)

            with open(os.path.join(out_precursor_dir, "__index.tsv"), "a") as f:
                f.write(f"{charge}\t{charge}\n")
