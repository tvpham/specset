import os
import io
import re
import shutil
import pandas as pd
from pathlib import Path
from collections import defaultdict
from typing import Iterator, Optional, Dict, List, Tuple, Any, Iterable


from .mz import get_ion_mz, get_precursor_mz, generate_aa_mass

_FRAG_BAD = {"", "acK", "acK*", "H13 C8 N O5"}
_FRAG_SKIP = {"pY", "a2"}
_FRAG_PATTERN = re.compile(r'^([a-z])([0-9]+)(\*|\-H2O|\-NH3|\-CO|\-H4COS)?(\([0-9]+\+\))?$')
_FRAG_LOSS_MAP = {
    None: "noloss", "": "noloss",
    "*": "H3PO4", "-NH3": "NH3", "-H2O": "H2O", "-CO": "CO", "-H4COS": "H4COS",
}


def iter_specset_locations(specset_dir: str, n: Optional[int] = None) -> Iterator[Tuple[str, str]]:
    """
    Yields (peptide, location) pairs from the main _index.tsv file.

    Parameters
    ----------
    specset_dir : str
        Path to the specset directory.
    n : int, optional
        Maximum number of entries to yield. If None, yields all entries.

    Returns
    -------
    Iterator[Tuple[str, str]]
        An iterator yielding tuples of (peptide, location).
    """
    i = 0
    index_path = Path(specset_dir) / "_index.tsv"

    with open(index_path) as index_file:
        for line in index_file:
            if n is not None and i >= n: 
                break
            line = line.strip()
            if not line: 
                continue
            parts = line.split(maxsplit=1)
            if len(parts) == 2:
                yield parts[0], parts[1]
                i += 1

def iter_specset_files(specset_dir: str, n: Optional[int] = None) -> Iterator[Tuple[str, str, Optional[str]]]:
    """
    Yields (location, file_name, charge) for each individual PSM file.

    Parameters
    ----------
    specset_dir : str
        Path to the specset directory.
    n : int, optional
        Maximum number of peptides to process from the main index.

    Returns
    -------
    Iterator[Tuple[str, str, Optional[str]]]
        An iterator yielding tuples of (location, file_name, charge).
    """
    specset_path = Path(specset_dir)

    for pep, location in iter_specset_locations(specset_dir, n):
        index2_path = specset_path / location / "__index.tsv"
        with open(index2_path) as index2:
            for line2 in index2:
                line2 = line2.strip()
                if not line2: 
                    continue
                parts = line2.split(maxsplit=1)
                file_name = parts[0]
                charge = parts[1] if len(parts) > 1 else None
                yield location, file_name, charge

def matches_intensities_to_psm(matches: str, intensities: str, sep: str = ";", max_intensity: float = -1.0) -> Optional[pd.DataFrame]:
    """
    Parse MaxQuant fragment matches and intensities into a DataFrame.

    losses such as *-NH3 are ignored.
    example:
        matches = "b3;y3-NH3;b4;y3;y4;b5;y5*;b6;b7;y5;b8;y17*(2+);b9;y18*(2+);y9*;y10*;y9;y10;y12*-NH3;y12*;y11;y13*-NH3;y13*;y12;y14*;y13;y15*-NH3;y15*;y14;b16;y15;y17*-NH3;y16;y18*-NH3;y17"
        intensities = "64.1662826538086;7.858124732971191;53.30729675292969;10.5747709274292;66.02859497070312;117.34032440185547;5.844237804412842;133.63040161132812;58.8831672668457;8.286745071411133;54.73167419433594;18.573741912841797;31.438446044921875;11.144254684448242;10.002944946289062;22.58058738708496;6.572253704071045;4.429126739501953;7.7240190505981445;34.01408767700195;6.000754356384277;31.58258056640625;25.00667953491211;19.00553321838379;18.720115661621094;27.294437408447266;7.71705436706543;15.690030097961426;9.860227584838867;3.571875810623169;28.576702117919922;6.858004570007324;14.86058235168457;4.143375396728516;6.000750541687012"

    Parameters
    ----------
    matches : str
        A delimited string of fragment matches (e.g., "b3;y3-NH3;b4").
    intensities : str
        A delimited string of corresponding fragment intensities.
    sep : str, optional
        The delimiter used in the strings. Default is ";".
    max_intensity : float, optional
        If greater than 0, normalizes relative intensities to this maximum value.

    Returns
    -------
    pd.DataFrame or None
        A DataFrame containing the parsed fragments and intensities, or None if parsing fails.

    """
    fragments = matches.split(sep)

    if _FRAG_BAD & set(fragments):
        return None

    try:
        inten_values = [float(x) for x in intensities.split(sep)]
    except ValueError:
        return None

    if any(x != x for x in inten_values):  # NaN check
        return None    

    pairs = [(f, i) for f, i in zip(fragments, inten_values) if f not in _FRAG_SKIP]
    if not pairs:
        return None
    fragments, inten_values = zip(*pairs)

    rows = []
    for frag, inten in zip(fragments, inten_values):
        m = _FRAG_PATTERN.match(frag)
        if not m:
            continue
        frag_type, frag_num, frag_loss, frag_charge = m.groups()
        frag_charge = frag_charge.lstrip("(").rstrip("+)") if frag_charge else "1"
        rows.append({
            "FragmentLossType": _FRAG_LOSS_MAP.get(frag_loss, frag_loss),
            "FragmentNumber": int(frag_num),
            "FragmentType": frag_type,
            "FragmentCharge": frag_charge,
            "RelativeIntensity": inten,
        })

    if not rows:
        return None

    df = pd.DataFrame(rows)

    if max_intensity > 0:
        df["RelativeIntensity"] = df["RelativeIntensity"] / df["RelativeIntensity"].max() * max_intensity

    return df

class _specset_writer:
    """
    Write PSMs one at a time to the specset directory format.

    Call add() for each PSM and close() at the end. See _psm_stream_to_specset for the parameters.
    """

    _CORE_KEYS = {"pep", "z", "rt", "im", "matches", "intensities"}

    def __init__(self, out_dir: str, max_files_in_directory: int = 5000, max_intensity: float = -1.0, flush_threshold: int = 1000, overwrite: bool = False) -> None:
        self.out_path = Path(out_dir)
        if self.out_path.exists() and any(self.out_path.iterdir()):
            if not overwrite:
                raise FileExistsError(f"out_dir is not empty: {self.out_path}")
            print(f"Output path exists. Removing {self.out_path}...")
            shutil.rmtree(self.out_path)
        print(f"Creating {self.out_path}...")
        self.out_path.mkdir(parents=True, exist_ok=True)

        self.max_files_in_directory = max_files_in_directory
        self.max_intensity = max_intensity
        self.flush_threshold = flush_threshold

        self.pep_dict = {}
        self.i = 0
        self.folder_name = f"0-{max_files_in_directory - 1}"
        (self.out_path / self.folder_name).mkdir(exist_ok=True)

        self.index_buffer = defaultdict(list)
        self.buffer_size = 0

    def _flush_index_buffer(self) -> None:
        for dir_path, lines in self.index_buffer.items():
            with open(dir_path / "__index.tsv", "a") as f:
                f.write("".join(lines))
        self.index_buffer.clear()
        self.buffer_size = 0

    def add(self, psm: Dict[str, Any]) -> None:
        """Write one PSM; PSMs without valid fragments are skipped."""
        d = matches_intensities_to_psm(psm["matches"], psm["intensities"], max_intensity=self.max_intensity)

        if d is None:
            return

        pep = psm["pep"]

        specset = {
            "pep": pep,
            "z":   psm["z"],
            "ms2": d.to_csv(sep=",", lineterminator=";", index=False, header=False),
            "rt": psm["rt"],
            "im": psm["im"],
        }

        for key, value in psm.items():
            if key not in self._CORE_KEYS:
                specset[key] = value

        if pep not in self.pep_dict:
            self.pep_dict[pep] = f"{self.folder_name}/{self.i}"
            dir_path = self.out_path / self.pep_dict[pep]
            dir_path.mkdir(parents=True, exist_ok=True)

            content_name = dir_path / "__index.tsv"
            if content_name.exists():
                content_name.unlink()

        dir_path = self.out_path / self.pep_dict[pep]

        with open(dir_path / f"{self.i}.tsv", "w", newline="") as f:
            for key, value in specset.items():
                f.write(f"{key}\t{value}\n")

        self.index_buffer[dir_path].append(f"{self.i}\t{psm['z']}\n")
        self.buffer_size += 1

        if self.buffer_size >= self.flush_threshold:
            self._flush_index_buffer()

        self.i += 1

        if self.i % self.max_files_in_directory == 0:
            self.folder_name = f"{self.i}-{self.i + self.max_files_in_directory - 1}"
            (self.out_path / self.folder_name).mkdir(exist_ok=True)

    def close(self) -> None:
        """Flush the remaining index lines and write _index.tsv."""
        if self.index_buffer:
            self._flush_index_buffer()

        with open(self.out_path / "_index.tsv", "w", newline="") as f:
            for key, value in self.pep_dict.items():
                f.write(f"{key}\t{value}\n")


def _psm_stream_to_specset(psm_iter: Iterable[Dict[str, Any]], out_dir: str, max_files_in_directory: int = 5000, max_intensity: float = -1.0, flush_threshold: int = 1000, overwrite: bool = False) -> None:
    """
    Write a stream of PSMs to the specset directory format.
    
    Parameters
    ----------
    psm_iter : Iterable[Dict[str, Any]]
        Each dict must have keys: pep, z, rt, im, matches, intensities.
    out_dir : str
        Output specset directory path.
    max_files_in_directory : int, optional
        Number of PSM files permitted per batch folder. Default is 5000.
    max_intensity : float, optional
        Max intensity for normalization. Default is -1.0 (no normalization).
    flush_threshold : int, optional
        Number of PSMs to buffer before flushing to disk. Default is 1000.
    overwrite : bool, optional
        If True, delete an existing out_dir before writing. If False, raise
        FileExistsError when out_dir exists and is not empty. Default is False.

    Returns
    -------
    None
        Writes directly to disk.
    """
    writer = _specset_writer(out_dir, max_files_in_directory, max_intensity, flush_threshold, overwrite)
    for psm in psm_iter:
        writer.add(psm)
    writer.close()


def psm_utils_to_specset(msms_file: str, out_dir: str, evidence_file: Optional[str] = None, max_files_in_directory: int = 5000, overwrite: bool = False) -> None:
    """
    Read a MaxQuant msms.txt file and write it to the specset directory format.

    Parameters
    ----------
    msms_file : str
        Path to the MaxQuant msms.txt file.
    out_dir : str
        Output directory path for the specset.
    evidence_file : str, optional
        Path to an optional evidence.txt file to extract ion mobility (1/K0).
    max_files_in_directory : int, optional
        Maximum number of PSMs per folder batch. Default is 5000.
    overwrite : bool, optional
        If True, delete an existing out_dir before writing. If False, raise
        FileExistsError when out_dir exists and is not empty. Default is False.

    Returns
    -------
    None
        Writes directly to disk.
    """

    from psm_utils.io.maxquant import MSMSReader
    
    psm_list = MSMSReader(msms_file).read_file()

    evidence_file_df = None
    if evidence_file is not None:
        evidence_file_df = pd.read_csv(evidence_file, sep="\t", index_col="id", dtype={"id": str}, low_memory=False)

    def psm_iter():
        for psm in psm_list:
            pep = (psm.peptidoform.modified_sequence
                   .replace("[Phospho (STY)]", "(UniMod:21)")
                   .replace("[Acetyl (Protein N-term)]-", "(UniMod:1)")
                   .replace("[Oxidation (M)]", "(UniMod:35)")
                   .replace("C", "C(UniMod:4)"))
            im = 0 if evidence_file_df is None else evidence_file_df.loc[psm.metadata["Evidence ID"]]["1/K0"]
            yield {
                "pep":         pep,
                "z":           psm.peptidoform.precursor_charge,
                "rt":          psm.retention_time,
                "im":          im,
                "score":       psm.score,
                "id":          psm.run,
                "matches":     psm.metadata["Matches"],
                "intensities": psm.metadata["Intensities"],
            }

    _psm_stream_to_specset(psm_iter(), out_dir, max_files_in_directory, overwrite=overwrite)


def maxquant_to_specset(txt_dir: str, out_dir: str, additional_columns: Dict[str, str] = {}, max_files_in_directory: int = 5000, n: int = -1, overwrite: bool = False) -> None:
    """
    Read MaxQuant msms.txt and evidence.txt files and write to specset format.

    Parameters
    ----------
    txt_dir : str
        Directory containing the msms.txt and evidence.txt files.
    out_dir : str
        Output directory for the specset.
    additional_columns : dict, optional
        Keys are column names in msms.txt, values are names to use in the specset dict.
        Example: {"Gene Names": "genes"}
    max_files_in_directory : int, optional
        Max number of PSMs per batch folder. Default is 5000.
    n : int, optional
        Maximum number of rows to process. Default is -1 (all rows).
    overwrite : bool, optional
        If True, delete an existing out_dir before writing. If False, raise
        FileExistsError when out_dir exists and is not empty. Default is False.

    Returns
    -------
    None
        Writes directly to disk.
    """

    print("Reading msms.txt...")
    msms_df = pd.read_csv(Path(txt_dir) / "msms.txt", sep="\t", usecols=["Modified sequence", "Charge", "Retention time", "Matches", "Intensities", "Evidence ID"] + list(additional_columns.keys()), dtype=str, keep_default_na=False)
    msms_df = msms_df.rename(columns=lambda c: c.replace(" ", "_"))

    print("Reading evidence.txt...")
    evidence_df = pd.read_csv(Path(txt_dir) / "evidence.txt", sep="\t", usecols=["id", "1/K0"], index_col="id", dtype={"id": str})

    has_im = "1/K0" in evidence_df.columns
    has_reverse = "Reverse" in msms_df.columns

    n_rows = 0

    def psm_iter():
        nonlocal n_rows
        for psm in msms_df.itertuples(index=False):

            if n > 0 and n_rows >= n:
                break

            if has_reverse and psm.Reverse == "+":
                continue

            n_rows += 1

            if n_rows % 1000 == 0:
                print(f"n = {n_rows}", end = '\r')

            pep = (psm.Modified_sequence[1:-1].replace("(Phospho (STY))", "(UniMod:21)")
                   .replace("(Acetyl (Protein N-term))", "(UniMod:1)")
                   .replace("(Oxidation (M))", "(UniMod:35)")
                   .replace("C", "C(UniMod:4)"))

            im = evidence_df.loc[str(psm.Evidence_ID)]["1/K0"] if has_im else 0
            
            d = {
                "pep":         pep,
                "z":           psm.Charge,
                "rt":          psm.Retention_time,
                "im":          im,
                "matches":     psm.Matches,
                "intensities": psm.Intensities
            }
            for msms_col, d_key in additional_columns.items():
                d[d_key] = getattr(psm, msms_col.replace(" ", "_"), None)
            yield d

    _psm_stream_to_specset(psm_iter(), out_dir, max_files_in_directory, overwrite=overwrite)

    print(f"n = {n_rows}")
    print("Done.")


def easypqp_to_specset(easypqp_file: str, out_dir: str, max_no_files_per_directory: int = 5000, overwrite: bool = False) -> None:
    """
    Convert an EasyPQP formatted file to the specset directory format.

    Parameters
    ----------
    easypqp_file : str
        Path to the EasyPQP TSV file.
    out_dir : str
        Output specset directory path.
    max_no_files_per_directory : int, optional
        Maximum number of files per batch directory. Default is 5000.
    overwrite : bool, optional
        If True, delete an existing out_dir before writing. If False, raise
        FileExistsError when out_dir exists and is not empty. Default is False.

    Returns
    -------
    None
        Writes directly to disk.
    """
    df = pd.read_csv(easypqp_file, sep="\t", usecols=["ModifiedPeptideSequence", "PrecursorCharge", "FragmentLossType", "FragmentSeriesNumber", "FragmentType", "FragmentCharge", "LibraryIntensity", "NormalizedRetentionTime", "PrecursorIonMobility", "ProteinId", "GeneName"])
    total_rows = len(df)
    print(f"Total rows: {total_rows}")

    loss_to_mq = {"noloss": "", "H3PO4": "*", "NH3": "-NH3", "H2O": "-H2O"}

    def psm_iter():
        rows_done = 0
        report_step = max(total_rows // 20, 1)
        next_report = report_step

        for (pep, z), group in df.groupby(["ModifiedPeptideSequence", "PrecursorCharge"]):
            matches = []
            intensities = []
            for row in group.itertuples(index=False):
                loss = loss_to_mq.get(row.FragmentLossType, "")
                charge_str = f"({row.FragmentCharge}+)" if int(row.FragmentCharge) > 1 else ""
                matches.append(f"{row.FragmentType}{row.FragmentSeriesNumber}{loss}{charge_str}")
                intensities.append(str(row.LibraryIntensity))
            yield {
                "pep":         pep,
                "z":           z,
                "rt":          group["NormalizedRetentionTime"].iloc[0],
                "im":          group["PrecursorIonMobility"].iloc[0],
                "matches":     ";".join(matches),
                "intensities": ";".join(intensities),
                "ProteinId":          group["ProteinId"].iloc[0],
                "GeneName":          group["GeneName"].iloc[0],
            }

            rows_done += len(group)
            if rows_done >= next_report or rows_done == total_rows:
                print(f"  {rows_done}/{total_rows} rows ({100*rows_done//total_rows}%)", end="\r")
                next_report += report_step

    _psm_stream_to_specset(psm_iter(), out_dir, max_no_files_per_directory, overwrite=overwrite)

    print("\nDone.")


def diannlib_to_specset(diannlib_file: str, out_dir: str, max_no_files_per_directory: int = 5000, overwrite: bool = False) -> None:
    """
    Convert a DIA-NN library parquet file to the specset directory format.

    Parameters
    ----------
    diannlib_file : str
        Path to the DIA-NN library parquet file.
    out_dir : str
        Output specset directory path.
    max_no_files_per_directory : int, optional
        Maximum number of files per batch directory. Default is 5000.
    overwrite : bool, optional
        If True, delete an existing out_dir before writing. If False, raise
        FileExistsError when out_dir exists and is not empty. Default is False.

    Returns
    -------
    None
        Writes directly to disk.
    """

    import pyarrow.parquet as pq

    pf = pq.ParquetFile(diannlib_file)
    total_rows = pf.metadata.num_rows
    print(f"Total rows: {total_rows}")

    loss_to_mq = {"noloss": "", "H3PO4": "*", "NH3": "-NH3", "H2O": "-H2O"}

    def _yield_group(rows):
        matches = []
        intensities = []
        
        for row in rows:
            loss = loss_to_mq.get(row.Fragment_Loss_Type, "")
            charge_str = f"({row.Fragment_Charge}+)" if int(row.Fragment_Charge) > 1 else ""
            matches.append(f"{row.Fragment_Type}{row.Fragment_Series_Number}{loss}{charge_str}")
            intensities.append(str(row.Relative_Intensity))
            
        first = rows[0]
        return {
            "pep":         first.Modified_Sequence,
            "z":           first.Precursor_Charge,
            "rt":          first.RT,
            "im":          first.IM,
            "matches":     ";".join(matches),
            "intensities": ";".join(intensities),
            "proteins":    first.Protein_Ids,
            "genes":       first.Genes,
        }

    def psm_iter():
        current_id = None
        current_rows = []
        rows_done = 0

        for batch in pf.iter_batches(batch_size=10000):
            chunk = batch.to_pandas()
            chunk.columns = chunk.columns.str.replace(".", "_", regex=False)
            for row in chunk.itertuples(index=False):
                pid = row.Precursor_Id
                
                if pid != current_id:
                    if current_rows:
                        yield _yield_group(current_rows)
                    current_id = pid
                    current_rows = [row]
                else:
                    current_rows.append(row)

            rows_done += len(chunk)
            print(f"  {rows_done}/{total_rows} rows ({100*rows_done//total_rows}%)", end="\r")

        if current_rows:
            yield _yield_group(current_rows)        

    _psm_stream_to_specset(psm_iter(), out_dir, max_no_files_per_directory, overwrite=overwrite)
    
    print("\nDone.")


def diannreport_to_specset(diannreport: str, out_dir: str, column: str = "", values: List[str] = [], max_no_files_per_directory: int = 5000, overwrite: bool = False, by_Run: bool = False) -> None:
    """
    Convert a DIA-NN report parquet file to a specset, or to one specset per run.

    The report is read once. Each row (one precursor in one run) becomes one PSM. With by_Run,
    the PSMs of each run go to their own specset, out_dir/<Run>.specset, where the characters
    / \\ : * ? " < > | and whitespace in the run name are replaced by "_" (and _2, _3, ... is
    added if two runs get the same name). Fragments are taken from the Fr.N.Id and
    Fr.N.Quantity columns, e.g. "y8-unknown^1/816.403259"; fragments with an empty id or a
    zero quantity are skipped. Decoy rows are skipped.

    Parameters
    ----------
    diannreport : str
        Path to the DIA-NN report parquet file.
    out_dir : str
        Output specset directory; with by_Run, the directory in which one <Run>.specset is created per run.
    column : str, optional
        Report column to select rows on, e.g. "Run". Default is "" (import all rows).
    values : list of str, optional
        With column, import only the rows whose value in column is one of these values
        (compared as text).
    max_no_files_per_directory : int, optional
        Maximum number of files per batch directory. Default is 5000.
    overwrite : bool, optional
        If True, delete an existing output specset before writing it. If False, raise
        FileExistsError when an output specset exists and is not empty. Default is False.
    by_Run : bool, optional
        If True, create one specset per run in out_dir. Default is False (one specset with all runs).

    Returns
    -------
    None
        Writes directly to disk.
    """

    import pyarrow.parquet as pq

    pf = pq.ParquetFile(diannreport)
    total_rows = pf.metadata.num_rows
    print(f"Total rows: {total_rows}")

    n_frags = sum(1 for name in pf.schema_arrow.names if re.fullmatch(r"Fr\.\d+\.Id", name))
    columns = ["Run", "Modified.Sequence", "Precursor.Charge", "RT", "IM", "Q.Value", "Decoy", "Ms2.Scan", "Protein.Ids", "Genes"]
    columns += [f"Fr.{i}.{field}" for i in range(n_frags) for field in ("Id", "Quantity")]
    if column:
        if column not in pf.schema_arrow.names:
            raise ValueError(f"column '{column}' not found in {diannreport}")
        if column not in columns:
            columns.append(column)
        selected = {str(v) for v in values}

    # fragment id: type, number, loss, charge, m/z, e.g. "y8-unknown^1/816.403259"
    frag_id_pattern = re.compile(r"^([a-z])(\d+)-([^^]+)\^(\d+)/")
    loss_to_mq = {"unknown": "", "noloss": "", "H3PO4": "*", "NH3": "-NH3", "H2O": "-H2O"}

    out_path = Path(out_dir)
    writers = {}  # run (None without by_Run) -> _specset_writer
    used_names = set()

    def specset_name(run):
        """Specset directory name for a run: special characters replaced by '_', made unique."""
        name = re.sub(r'[\\/:*?"<>|\s]', "_", str(run))
        unique, n = name, 1
        while unique in used_names:
            n += 1
            unique = f"{name}_{n}"
        used_names.add(unique)
        return f"{unique}.specset"

    if not by_Run:
        writers[None] = _specset_writer(out_path, max_no_files_per_directory, overwrite=overwrite)
    rows_done = 0

    try:
        for batch in pf.iter_batches(batch_size=10000, columns=columns):
            chunk = batch.to_pandas()
            for row_values in chunk.itertuples(index=False, name=None):
                row = dict(zip(columns, row_values))
                if row["Decoy"]:
                    continue
                if column and str(row[column]) not in selected:
                    continue

                matches = []
                intensities = []
                for i in range(n_frags):
                    frag_id = row[f"Fr.{i}.Id"]
                    quantity = row[f"Fr.{i}.Quantity"]
                    if not frag_id or not quantity or quantity != quantity:  # empty, zero or NaN
                        continue
                    m = frag_id_pattern.match(frag_id)
                    if m is None or m.group(3) not in loss_to_mq:
                        continue
                    frag_type, frag_number, frag_loss, frag_charge = m.groups()
                    charge_str = f"({frag_charge}+)" if int(frag_charge) > 1 else ""
                    matches.append(f"{frag_type}{frag_number}{loss_to_mq[frag_loss]}{charge_str}")
                    intensities.append(str(quantity))

                if not matches:
                    continue

                run = row["Run"]
                key = run if by_Run else None
                if key not in writers:
                    print()
                    writers[key] = _specset_writer(out_path / specset_name(run), max_no_files_per_directory, overwrite=overwrite)

                writers[key].add({
                    "pep":         row["Modified.Sequence"],
                    "z":           row["Precursor.Charge"],
                    "rt":          row["RT"],
                    "im":          row["IM"],
                    "matches":     ";".join(matches),
                    "intensities": ";".join(intensities),
                    "run":         run,
                    "scan":        row["Ms2.Scan"],
                    "qvalue":      row["Q.Value"],
                    "proteins":    row["Protein.Ids"],
                    "genes":       row["Genes"],
                })

            rows_done += len(chunk)
            print(f"  {rows_done}/{total_rows} rows ({100*rows_done//total_rows}%)", end="\r")
    finally:
        for writer in writers.values():
            writer.close()

    if by_Run:
        print(f"\nDone. {len(writers)} specsets written to {out_path}")
    else:
        print(f"\nDone. Specset written to {out_path}")


def psmlist_to_specset(psmlist_file: str, out_dir: str, max_files_in_directory: int = 5000, overwrite: bool = False) -> None:
    """
    Read a list of files containing stripped MaxQuant PSMs and write to specset format.

    Parameters
    ----------
    psmlist_file : str
        Path to a file containing a list of PSM file paths (one per line).
    out_dir : str
        Output specset directory path.
    max_files_in_directory : int, optional
        Maximum number of files per batch directory. Default is 5000.
    overwrite : bool, optional
        If True, delete an existing out_dir before writing. If False, raise
        FileExistsError when out_dir exists and is not empty. Default is False.

    Returns
    -------
    None
        Writes directly to disk.
    """

    def psm_iter():
        with open(psmlist_file, "r") as file:
            for line in file:
                fname = line.strip()
                print(fname)
                df = pd.read_csv(fname, sep="\t")
                df.columns = df.columns.str.replace(".", "_", regex=False).str.replace(" ", "_", regex=False)
                for psm in df.itertuples():
                    yield {
                        "pep":         psm.Modified_sequence[1:-1],  
                        "z":           psm.Charge,
                        "rt":          psm.Retention_time,
                        "im":          0,
                        "score":       psm.Score,
                        "id":          psm.Search_ID,
                        "matches":     psm.Matches,
                        "intensities": psm.Intensities,
                    }

    _psm_stream_to_specset(psm_iter(), out_dir, max_files_in_directory, overwrite=overwrite)
    print("Done.")


def _specset_stream_to_string(specset_dir: str, n: Optional[int] = None) -> io.StringIO:
    """
    Convert a specset directory into a single StringIO buffer stream.

    Parameters
    ----------
    specset_dir : str
        Path to the specset directory.
    n : int, optional
        Maximum number of entries to process. If None, processes all.

    Returns
    -------
    io.StringIO
        A string buffer containing the concatenated tabular data.
    """

    _cols = ["FragmentLossType", "FragmentNumber", "FragmentType", "FragmentCharge", "RelativeIntensity"]

    string_buffer = io.StringIO()
    first = True
    specset_path = Path(specset_dir)

    for location, file_name, _ in iter_specset_files(specset_dir, n):
        a = {}
        psm_path = specset_path / location / f"{file_name}.tsv"

        with open(psm_path) as psm_file:
            for psm_line in psm_file:
                k, v = psm_line.rstrip("\n").split(maxsplit=1)
                a[k] = v

        if "ms2" in a:
            frag_df = pd.read_csv(io.StringIO(a["ms2"].rstrip(";")), sep=",", lineterminator=";", header=None, names=_cols)
        else:
            frag_df = pd.DataFrame(columns=_cols)

        frag_df["ModifiedPeptide"] = a["pep"]
        frag_df["PrecursorCharge"] = int(a["z"])
        frag_df["RT"] = float(a["rt"]) if "rt" in a else 0
        frag_df["IM"] = float(a["im"]) if "im" in a else 0

        frag_df.to_csv(string_buffer, sep="\t", index=False, header=first)
        first = False

    return string_buffer

def specset_to_alphapeptdeep_string(specset_dir: str, n: Optional[int] = None) -> Iterator[io.StringIO]:
    """
    Yields string buffers mapped to the AlphaPeptDeep spec format.

    Parameters
    ----------
    specset_dir : str
        Path to the specset directory.
    n : int, optional
        Maximum number of entries to process. If None, processes all.

    Returns
    -------
    Iterator[io.StringIO]
        An iterator of string buffers representing AlphaPeptDeep fragments.
    """

    _cols = ["FragmentLossType", "FragmentNumber", "FragmentType", "FragmentCharge", "RelativeIntensity"]

    
    specset_path = Path(specset_dir)

    for location, file_name, _ in iter_specset_files(specset_dir, n):
        a = {}
        psm_path = specset_path / location / f"{file_name}.tsv"

        with open(psm_path) as psm_file:
            for psm_line in psm_file:
                parts = psm_line.rstrip("\n").split(maxsplit=1)
                if len(parts) == 2:
                    a[parts[0]] = parts[1]

        if "ms2" in a:
            frag_df = pd.read_csv(io.StringIO(a["ms2"].rstrip(";")), sep=",", lineterminator=";", header=None, names=_cols)
        else:
            frag_df = pd.DataFrame(columns=_cols)

        frag_df["ModifiedPeptide"] = a["pep"]
        frag_df["PrecursorCharge"] = int(a["z"])
        frag_df["RT"] = float(a["rt"]) if "rt" in a else 0
        frag_df["IonMobility"] = float(a["im"]) if "im" in a else 0

        string_buffer = io.StringIO()
        frag_df.to_csv(string_buffer, sep="\t", index=False, header=True)

        yield string_buffer
        
    
def _read_location(location_dir: str, extra: List[str] = []) -> Dict[str, List[pd.DataFrame]]:
    """
    Read all PSM files in a specific specset location directory.

    Parameters
    ----------
    location_dir : str
        Path to the specific location directory inside the specset.
    extra : list of str, optional
        Additional columns to extract from the PSM files.

    Returns
    -------
    Dict[str, List[pd.DataFrame]]
        A dictionary mapping precursor charge to a list of fragment DataFrames.
    """

    loc_path = Path(location_dir)
    result = {}

    with open(loc_path / "__index.tsv") as index2_file:
        for line2 in index2_file:
            line2 = line2.rstrip("\n")
            if not line2:
                continue
            parts = line2.split(maxsplit=1)
            file_name = parts[0]
            z = parts[1]
            if z not in result:
                result[z] = []

            a = {}
            with open(loc_path / f"{file_name}.tsv") as psm_file:
                for psm_line in psm_file:
                    parts = psm_line.rstrip("\n").split(maxsplit=1)
                    if len(parts) == 2:
                        a[parts[0]] = parts[1]
                    elif len(parts) == 1:
                        a[parts[0]] = ""

            _cols = ["FragmentLossType", "FragmentNumber", "FragmentType", "FragmentCharge", "RelativeIntensity"]

            if "ms2" in a:
                frag_df = pd.read_csv(io.StringIO(a["ms2"].rstrip(";")), sep=",", lineterminator=";", header=None, names=_cols)
            else:
                frag_df = pd.DataFrame(columns=_cols)

            frag_df["ModifiedPeptide"] = a["pep"]
            frag_df["PrecursorCharge"] = int(a["z"])
            frag_df["RT"] = float(a["rt"]) if "rt" in a else 0
            frag_df["IM"] = float(a["im"]) if "im" in a else 0

            for col in extra:
                frag_df[col] = a.get(col)

            result[z].append(frag_df)

    return result


def _consensus(frag_dfs: List[pd.DataFrame], consensus_method: Optional[str] = None) -> pd.DataFrame:
    """
    Select or aggregate fragment DataFrames using a consensus method.

    Parameters
    ----------
    frag_dfs : List[pd.DataFrame]
        A list of fragment DataFrames for a particular precursor.
    consensus_method : str, optional
        'mean': average RelativeIntensity across PSMs per fragment.
        'median': median RelativeIntensity across PSMs per fragment.
        any other string: select by max value of that column.
        None: select by most fragments.

    Returns
    -------
    pd.DataFrame
        A single, merged DataFrame representing the consensus.
    """

    if consensus_method == "mean" or consensus_method == "median":
        combined = pd.concat(frag_dfs, ignore_index=True)
        group_cols = ["FragmentLossType", "FragmentNumber", "FragmentType", "FragmentCharge"]
        agg = combined.groupby(group_cols, as_index=False)["RelativeIntensity"]
        agg = agg.mean() if consensus_method == "mean" else agg.median()
        meta = frag_dfs[0].drop(columns=["RelativeIntensity"]).iloc[[0]]
        return agg.merge(meta, on=group_cols, how="left")
    elif consensus_method is not None and consensus_method in frag_dfs[0].columns:
        return max(frag_dfs, key=lambda df: float(df[consensus_method].iloc[0]))
    else:
        return max(frag_dfs, key=lambda df: len(df))


def _to_mz_seq(pep: str) -> str:
    """
    Convert UniMod peptide notation to mz.py single-char encoding.

    Parameters
    ----------
    pep : str
        Modified sequence in UniMod notation.

    Returns
    -------
    str
        Simplified single-character sequence encoding.
    """
    seq = pep.replace("(UniMod:1)", "*")
    seq = re.sub(r'M\(UniMod:35\)', '1', seq)
    seq = re.sub(r'S\(UniMod:21\)', '2', seq)
    seq = re.sub(r'T\(UniMod:21\)', '3', seq)
    seq = re.sub(r'Y\(UniMod:21\)', '4', seq)
    seq = re.sub(r'C\(UniMod:4\)', 'C', seq)
    return seq


def specset_to_speclib(specset_dir: str, lib_name: str, additional_columns: Dict[str, str] = {}, consensus_method: Optional[str] = None, n: Optional[int] = None) -> None:
    """
    Convert a specset directory to a consolidated spectral library TSV.

    Parameters
    ----------
    specset_dir : str
        Path to the specset directory.
    lib_name : str
        Output file path for the spectral library TSV.
    additional_columns : dict, optional
        Keys are specset field names, values are column names to use in the library.
    consensus_method : str, optional
        Aggregation method ('mean', 'median', or column name).
    n : int, optional
        Maximum number of entries to process.

    Returns
    -------
    None
        Writes the library directly to disk.
    """

    # Calculate the amino acid mass values upon import of the library
    # as they will be used frequently
    aa_mass = generate_aa_mass()
    specset_path = Path(specset_dir)

    #dia-nn header
    #cols = [
    #    "Precursor.Id",
    #    "Stripped.Sequence",
    #    "Precursor.Charge",
    #    "Precursor.Mz",
    #    "Protein.Name",
    #    "Precursor.NormalisedRetentionTime",
    #    "Fragment.Ion",
    #    "Fragment.Charge",
    #    "Fragment.Mz",
    #    "Fragment.Intensity",
    #    "Q.Value",
    #    "Global.Q.Value"        
    #]

    #alphabase header
    cols = ["ModifiedPeptide",
            "PrecursorCharge",
            "RT",
            "IM",
            "FragmentType",
            "RelativeIntensity",
            "FragmentCharge",
            "FragmentNumber",
            "FragmentLossType",
            "PrecursorMz",
            "FragmentMz"] + list(additional_columns.values())

    if consensus_method is not None and consensus_method not in ("mean", "median") and consensus_method not in cols:
        cols.append(consensus_method)

    index_path = specset_path / "_index.tsv"
    with open(index_path) as f:
        total = sum(1 for line in f if line.strip())
    print(f"Number of peptides = {total}")

    n_peps = 0

    with open(lib_name, "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        
        for pep, location in iter_specset_locations(specset_dir, n):
            location_path = os.path.join(specset_dir, location)
            for z, frag_dfs in _read_location(location_path, extra=list(additional_columns.keys())).items():

                df = _consensus(frag_dfs, consensus_method)
                pep_seq = _to_mz_seq(df.iloc[0]["ModifiedPeptide"])
                
                if not all(c in aa_mass for c in pep_seq):
                    continue
                    
                precursor_mz = get_precursor_mz(pep_seq, int(z), aa_mass)

                for row in df.itertuples(index=False):
                    frag_mz = get_ion_mz(pep_seq, row.FragmentType, row.FragmentNumber, int(row.FragmentCharge), row.FragmentLossType, aa_mass)
                    r = ["_" + row.ModifiedPeptide.replace("(", "[").replace(")", "]") + "_",
                         str(row.PrecursorCharge),
                         str(row.RT),
                         str(row.IM),
                         row.FragmentType,
                         str(row.RelativeIntensity),
                         str(row.FragmentCharge),
                         str(row.FragmentNumber),
                         row.FragmentLossType,
                         str(precursor_mz),
                         str(frag_mz)] + [str(getattr(row, col, "")) for col in additional_columns.keys()]

                    f.write("\t".join(r) + "\n")

            n_peps += 1
            print(f"  {n_peps}/{total} ({100*n_peps//total}%)", end="\r")
    
    print("\nDone.")


def _specset_to_string(specset_dir: str, scratch: Optional[int] = None) -> io.StringIO:
    """
    Convert a specset directory to a string buffer containing all PSMs.

    Parameters
    ----------
    specset_dir : str
        Path to the specset directory.
    scratch : int, optional
        Total length hint for progress printing. Re-calculated if None.

    Returns
    -------
    io.StringIO
        A string buffer containing the formatted tabular specset output.
    """
    specset_path = Path(specset_dir)

    if scratch is None:
        with open(specset_path / "_index.tsv") as f:
            scratch = sum(1 for line in f if line.strip())

    display = max(scratch // 20, 1)
    string_buffer = io.StringIO()
    first = True
    
    for i, (pep, location) in enumerate(iter_specset_locations(specset_dir, scratch)):
        if i % display == 0:
            print(f"{int((i+1) * 100 / scratch)}%")        

        location_path = specset_path / location
        index2_path = location_path / "__index.tsv"

        if not index2_path.exists(): 
            continue
            
        with open(index2_path) as index2:
            for line2 in index2:
                line2 = line2.strip()
                if not line2: 
                    continue
                
                file_name = line2.split(maxsplit=1)[0]
                psm_path = location_path / f"{file_name}.tsv"
                
                a = pd.read_csv(psm_path, 
                                header=None, sep="\t", index_col=0)
                
                frag_df = pd.read_json(io.StringIO(a.loc["ms2",1]), orient='values')
                frag_df.columns = ["FragmentLossType", "FragmentNumber", "FragmentType", "FragmentCharge", "RelativeIntensity"]
                frag_df["ModifiedPeptide"] = a.loc["pep", 1]
                frag_df["PrecursorCharge"] = int(a.loc["z", 1])
                frag_df["RT"] = float(a.loc["ret", 1])
                frag_df["CCS"] = float(a.loc["ccs", 1])

                frag_df.to_csv(string_buffer, sep="\t", index=False, header=first)
                first = False

    return string_buffer


def merge_specsets(specset_a: str, specset_b: str, specset_out: str, max_files_in_directory: int = 5000, mode: str = "ignore", flush_threshold: int = 1000, overwrite: bool = False) -> None:
    """
    Merge specset_a and specset_b into a new specset at specset_out.

    specset_a and specset_b are left unchanged. PSMs are renumbered from 0 in
    specset_out, PSMs from specset_a first, then those from specset_b.

    Parameters
    ----------
    specset_a : str
        Path to the primary specset directory.
    specset_b : str
        Path to the secondary specset directory.
    specset_out : str
        Path to the output specset directory.
    max_files_in_directory : int, optional
        Maximum number of files per batch directory. Default is 5000.
    mode : str, optional
        'both'   : for overlapping ones, keep both PSMs in specset_a and specset_b.
        'ignore' : for overlapping ones, keep PSMs in specset_a, skip the ones from specset_b. Default is 'ignore'.
    flush_threshold : int, optional
        Number of PSMs to buffer before flushing to disk. Default is 1000.
    overwrite : bool, optional
        If True, delete an existing specset_out before writing. If False, raise
        FileExistsError when specset_out exists and is not empty. Default is False.

    Returns
    -------
    None
        Writes directly to disk in specset_out.
    """
    if mode not in ("both", "ignore"):
        raise ValueError(f"mode must be 'both' or 'ignore', got {mode!r}")

    path_a = Path(specset_a)
    path_b = Path(specset_b)
    out_path = Path(specset_out)

    for path in (path_a, path_b):
        if path.resolve() == out_path.resolve() or out_path.resolve() in path.resolve().parents:
            raise ValueError(f"specset_out must not be or contain an input specset: {path}; use merge_specsets_inplace for in-place merging")
    if out_path.exists() and any(out_path.iterdir()):
        if not overwrite:
            raise FileExistsError(f"specset_out is not empty: {out_path}")
        print(f"Removing {out_path}...")
        shutil.rmtree(out_path)
    out_path.mkdir(parents=True, exist_ok=True)

    locations_a = list(iter_specset_locations(specset_a))
    locations_b = list(iter_specset_locations(specset_b))
    if mode == "ignore":
        peps_a = {pep for pep, _ in locations_a}
        locations_b = [(pep, loc) for pep, loc in locations_b if pep not in peps_a]
    sources = [(path_a, pep, loc) for pep, loc in locations_a] + [(path_b, pep, loc) for pep, loc in locations_b]

    display = max(len(sources) // 20, 1)

    pep_dict = {}
    index_buffer = defaultdict(list)
    buffer_size = 0
    i = 0

    def flush_index_buffer():
        for dir_path, lines in index_buffer.items():
            with open(dir_path / "__index.tsv", "a") as f:
                f.write("".join(lines))
        index_buffer.clear()
        return 0

    for k, (src_path, pep, location) in enumerate(sources):
        if k % display == 0:
            print(f"{int((k+1) * 100 / len(sources))}%", end='\r')

        with open(src_path / location / "__index.tsv") as index2:
            for line2 in index2:
                parts = line2.split()
                if not parts:
                    continue
                file_name, z = parts[0], parts[1]

                if pep not in pep_dict:
                    start = (i // max_files_in_directory) * max_files_in_directory
                    pep_dict[pep] = f"{start}-{start + max_files_in_directory - 1}/{i}"
                    (out_path / pep_dict[pep]).mkdir(parents=True)

                dir_path = out_path / pep_dict[pep]
                shutil.copyfile(src_path / location / f"{file_name}.tsv", dir_path / f"{i}.tsv")
                index_buffer[dir_path].append(f"{i}\t{z}\n")
                buffer_size += 1
                if buffer_size >= flush_threshold:
                    buffer_size = flush_index_buffer()
                i += 1

    flush_index_buffer()

    with open(out_path / "_index.tsv", "w", newline="") as f:
        for key, value in pep_dict.items():
            f.write(f"{key}\t{value}\n")

    print("\nDone.")

def merge_specsets_inplace(specset_a: str, specset_b: str, max_files_in_directory: int = 5000, mode: str = "ignore", flush_threshold: int = 1000) -> None:
    """
    Merge specset_b into specset_a in place.

    Peptides new to specset_a are added under a new location.

    Parameters
    ----------
    specset_a : str
        Path to the primary specset directory to merge into.
    specset_b : str
        Path to the secondary specset directory to merge from.
    max_files_in_directory : int, optional
        Maximum number of files per batch directory. Default is 5000.
    flush_threshold : int, optional
        Number of PSMs to buffer before flushing to disk. Default is 1000.
    mode : str, optional
        'both'   : for overlapping ones, keep both PSMs in specset_a and specset_b.
        'ignore' : for overlapping ones, keep PSMs in specset_a, skip the ones from specset_b. Default is 'ignore'.

    Returns
    -------
    None
        Writes directly to disk in specset_a.
    """

    path_a = Path(specset_a)
    path_b = Path(specset_b)
    index_file = path_a / "_index.tsv"

    if not index_file.exists():
        raise FileNotFoundError(f"specset_a has no _index.tsv: {index_file}")

    pep_dict = {}
    with open(index_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split(maxsplit=1)
            if len(parts) == 2:
                pep_dict[parts[0]] = parts[1]

    existing_peps = set(pep_dict)

    max_i = -1
    for location in pep_dict.values():
        index2_path = path_a / location / "__index.tsv"
        if not index2_path.exists():
            raise FileNotFoundError(f"specset_a location '{location}' has no __index.tsv: {index2_path}")
        with open(index2_path) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                max_i = max(max_i, int(line.split(maxsplit=1)[0]))

    i = max_i + 1
    start = (i // max_files_in_directory) * max_files_in_directory
    folder_name = f"{start}-{start + max_files_in_directory - 1}"
    (path_a / folder_name).mkdir(parents=True, exist_ok=True)

    index_buffer = defaultdict(list)
    buffer_size = 0

    def flush_index_buffer():
        for dir_path, lines in index_buffer.items():
            with open(dir_path / "__index.tsv", "a") as f:
                f.write("".join(lines))
        index_buffer.clear()
        return 0

    index_b = path_b / "_index.tsv"
    if not index_b.exists():
        raise FileNotFoundError(f"specset_b has no _index.tsv: {index_b}")

    with open(index_b) as f:
        scratch = sum(1 for line in f if line.strip())
    display = max(scratch // 20, 1)
    k = 0

    with open(index_b) as index2_file:
        for line in index2_file:
            line = line.rstrip("\n")
            if not line:
                continue
            if k % display == 0:
                print(f"{int((k+1) * 100 / scratch)}%")
            k += 1
            parts = line.split(maxsplit=1)
            if len(parts) != 2:
                continue
            location = parts[1]

            content_b = path_b / location / "__index.tsv"
            if not content_b.exists():
                continue

            with open(content_b) as content_file:
                for line2 in content_file:
                    line2 = line2.rstrip("\n")
                    if not line2:
                        continue
                    file_name = line2.split(maxsplit=1)[0]

                    psm = {}
                    psm_file_path = path_b / location / f"{file_name}.tsv"
                    with open(psm_file_path) as psm_file:
                        for psm_line in psm_file:
                            kv = psm_line.rstrip("\n").split("\t", maxsplit=1)
                            psm[kv[0]] = kv[1] if len(kv) == 2 else ""

                    if "pep" not in psm:
                        raise KeyError(f"PSM file missing required 'pep' field: {psm_file_path}")
                    if "z" not in psm:
                        raise KeyError(f"PSM file missing required 'z' field: {psm_file_path}")
                    pep = psm["pep"]

                    if mode == "ignore" and pep in existing_peps:
                        continue

                    if pep not in pep_dict:
                        pep_dict[pep] = f"{folder_name}/{i}"
                        (path_a / pep_dict[pep]).mkdir(parents=True, exist_ok=True)

                    dir_path = path_a / pep_dict[pep]

                    with open(dir_path / f"{i}.tsv", "w", newline="") as f:
                        for key, value in psm.items():
                            f.write(f"{key}\t{value}\n")

                    index_buffer[dir_path].append(f"{i}\t{psm['z']}\n")
                    buffer_size += 1

                    if buffer_size >= flush_threshold:
                        buffer_size = flush_index_buffer()

                    i += 1

                    if i % max_files_in_directory == 0:
                        folder_name = f"{i}-{i + max_files_in_directory - 1}"
                        (path_a / folder_name).mkdir(parents=True, exist_ok=True)

    if index_buffer:
        flush_index_buffer()

    with open(index_file, "w", newline="") as f:
        for key, value in pep_dict.items():
            f.write(f"{key}\t{value}\n")

def summarize(specset_a: str, specset_b: Optional[str] = None) -> None:
    """
    Print summary statistics for one or two specset directories.

    Parameters
    ----------
    specset_a : str
        Path to the first specset directory.
    specset_b : str, optional
        Path to a second specset directory. If given, also reports the
        peptide overlap and differences between the two specsets.

    Returns
    -------
    None
        Prints statistics to stdout.
    """
    peps_a = {pep for pep, _ in iter_specset_locations(specset_a)}
    print(f"\n{specset_a}: {len(peps_a)} modified peptides")

    if specset_b is None:
        return

    peps_b = {pep for pep, _ in iter_specset_locations(specset_b)}
    print(f"\n{specset_b}: {len(peps_b)} modified peptides")

    overlap = peps_a & peps_b
    only_a = peps_a - peps_b
    only_b = peps_b - peps_a

    print(f"\nOverlap: {len(overlap)} modified peptides")
    print(f"Only in {specset_a}: {len(only_a)} modified peptides")
    print(f"Only in {specset_b}: {len(only_b)} modified peptides")