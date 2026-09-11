import re
import pandas as pd

def hello_world():
    print("Hello, world")


_FRAG_BAD = {"", "acK", "acK*", "H13 C8 N O5"}
_FRAG_SKIP = {"pY", "a2"}
_FRAG_PATTERN = re.compile(r'^([a-z])([0-9]+)(\*|\-H2O|\-NH3|\-CO|\-H4COS)?(\([0-9]+\+\))?$')
_FRAG_LOSS_MAP = {
    None: "noloss", "": "noloss",
    "*": "H3PO4", "-NH3": "NH3", "-H2O": "H2O", "-CO": "CO", "-H4COS": "H4COS",
}


def matches_intensities_to_psm(matches, intensities, sep=";", max_intensity=-1):
    """Parse MaxQuant fragment matches and intensities into a DataFrame.

    example:
    matches = "b3;y3-NH3;b4;y3;y4;b5;y5*;b6;b7;y5;b8;y17*(2+);b9;y18*(2+);y9*;y10*;y9;y10;y12*-NH3;y12*;y11;y13*-NH3;y13*;y12;y14*;y13;y15*-NH3;y15*;y14;b16;y15;y17*-NH3;y16;y18*-NH3;y17"
    intensities = "64.1662826538086;7.858124732971191;53.30729675292969;10.5747709274292;66.02859497070312;117.34032440185547;5.844237804412842;133.63040161132812;58.8831672668457;8.286745071411133;54.73167419433594;18.573741912841797;31.438446044921875;11.144254684448242;10.002944946289062;22.58058738708496;6.572253704071045;4.429126739501953;7.7240190505981445;34.01408767700195;6.000754356384277;31.58258056640625;25.00667953491211;19.00553321838379;18.720115661621094;27.294437408447266;7.71705436706543;15.690030097961426;9.860227584838867;3.571875810623169;28.576702117919922;6.858004570007324;14.86058235168457;4.143375396728516;6.000750541687012"

    losses such as *-NH3 are ignored.

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


def _psm_stream_to_specset(psm_iter, out_dir, max_files_in_directory=5000, max_intensity=-1):
    """Write a stream of PSMs to the specset directory format.

    Parameters
    ----------
    psm_iter : iterable of dict
        Each dict must have keys: pep, z, rt, im, matches, intensities        
    out_dir : str
        Output specset directory.
    max_files_in_directory : int
        Number of PSM files per batch folder.
    """
    import os

    core_keys = {"pep", "z", "rt", "im", "matches", "intensities"}

    os.makedirs(out_dir, exist_ok=True)
    index_file = os.path.join(out_dir, "_index.tsv")

    pep_dict = {}
    i = 0
    start = 0
    folder_name = f"0-{max_files_in_directory - 1}"
    os.makedirs(os.path.join(out_dir, folder_name), exist_ok=True)

    for psm in psm_iter:
        
        d = matches_intensities_to_psm(psm["matches"], psm["intensities"], max_intensity = max_intensity)

        if d is not None:
            pep = psm["pep"]        
            specset = {
                "pep": pep,
                "z":   psm["z"],
                "ms2": d.to_csv(sep=",", lineterminator=";", index=False, header=False),
                "rt": psm["rt"],
                "im": psm["im"],
            }
            for key, value in psm.items():
                if key not in core_keys:
                    specset[key] = value

            if pep not in pep_dict:
                pep_dict[pep] = f"{folder_name}/{i}"
                dir_name = os.path.join(out_dir, pep_dict[pep])
                os.makedirs(dir_name, exist_ok=True)
                content_name = os.path.join(dir_name, "__index.tsv")
                if os.path.exists(content_name):
                    os.remove(content_name)

            dir_name = os.path.join(out_dir, pep_dict[pep])

            with open(os.path.join(dir_name, f"{i}.tsv"), "w", newline="") as f:
                for key, value in specset.items():
                    f.write(f"{key}\t{value}\n")

            with open(os.path.join(dir_name, "__index.tsv"), "a") as f:
                f.write(f"{i}\t{psm['z']}\n")

            i += 1
            if i % max_files_in_directory == 0:
                start = i
                folder_name = f"{start}-{start + max_files_in_directory - 1}"
                os.makedirs(os.path.join(out_dir, folder_name), exist_ok=True)

    with open(index_file, "w", newline="") as f:
        for key, value in pep_dict.items():
            f.write(f"{key}\t{value}\n")


def psm_utils_to_specset(msms_file, out_dir, evidence_file=None, max_files_in_directory=5000):
    """Read a MaxQuant msms.txt file and write it to the specset directory format."""
    import pandas as pd
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
                "rt":         psm.retention_time,
                "im":         im,
                "score":       psm.score,
                "id":          psm.run,
                "matches":     psm.metadata["Matches"],
                "intensities": psm.metadata["Intensities"],
            }

    _psm_stream_to_specset(psm_iter(), out_dir, max_files_in_directory)


def maxquant_to_specset(txt_dir, out_dir, additional_columns={}, max_files_in_directory=5000, n=-1):
    """Read a MaxQuant msms.txt and evidence.txt file and write it to the specset directory format.

    additional_columns : dict
        Keys are column names in msms.txt, values are the names to use in the specset dict.
        Example: {"Gene Names": "genes"}
    """
    import pandas as pd

    print("Reading msms.txt...")
    #msms_df = pd.read_csv(txt_dir.rstrip("/") + "/msms.txt", sep="\t", usecols=["Modified sequence", "Charge", "Retention time", "Matches", "Intensities", "Evidence ID"] + list(additional_columns.keys()), dtype={"Matches": str, "Intensities": str}, keep_default_na=False)
    msms_df = pd.read_csv(txt_dir.rstrip("/") + "/msms.txt", sep="\t", usecols=["Modified sequence", "Charge", "Retention time", "Matches", "Intensities", "Evidence ID"] + list(additional_columns.keys()), dtype=str, keep_default_na=False)
    msms_df = msms_df.rename(columns=lambda c: c.replace(" ", "_"))

    print("Reading evidence.txt...")
    evidence_df = pd.read_csv(txt_dir.rstrip("/") + "/evidence.txt", sep="\t", usecols=["id", "1/K0"], index_col="id", dtype={"id": str})

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

    _psm_stream_to_specset(psm_iter(), out_dir, max_files_in_directory)


def easypqp_to_specset(easypqp_file, out_dir, max_no_files_per_directory=5000):
    import pandas as pd

    df = pd.read_csv(easypqp_file, sep="\t", usecols=["ModifiedPeptideSequence", "PrecursorCharge", "FragmentLossType", "FragmentSeriesNumber", "FragmentType", "FragmentCharge", "LibraryIntensity", "NormalizedRetentionTime", "PrecursorIonMobility", "ProteinId", "GeneName"])

    loss_to_mq = {"noloss": "", "H3PO4": "*", "NH3": "-NH3", "H2O": "-H2O"}

    def psm_iter():
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

    _psm_stream_to_specset(psm_iter(), out_dir, max_no_files_per_directory)


def diannlib_to_specset(diannlib_file, out_dir, max_no_files_per_directory=5000):

    import pyarrow.parquet as pq
    import pandas as pd

    pf = pq.ParquetFile(diannlib_file)
    total_rows = pf.metadata.num_rows
    print(f"Total rows: {total_rows}")

    loss_to_mq = {"noloss": "", "H3PO4": "*", "NH3": "-NH3", "H2O": "-H2O"}

    def _yield_group(rows):
        group = pd.DataFrame(rows)
        matches = []
        intensities = []
        for row in group.itertuples(index=False):
            loss = loss_to_mq.get(row.Fragment_Loss_Type, "")
            charge_str = f"({row.Fragment_Charge}+)" if int(row.Fragment_Charge) > 1 else ""
            matches.append(f"{row.Fragment_Type}{row.Fragment_Series_Number}{loss}{charge_str}")
            intensities.append(str(row.Relative_Intensity))
        return {
            "pep":         group["Modified_Sequence"].iloc[0],
            "z":           group["Precursor_Charge"].iloc[0],
            "rt":          group["RT"].iloc[0],
            "im":          group["IM"].iloc[0],
            "matches":     ";".join(matches),
            "intensities": ";".join(intensities),
            "proteins": group["Protein_Ids"].iloc[0],
            "genes":       group["Genes"].iloc[0],
        }

    def psm_iter():
        current_id = None
        current_rows = []
        rows_done = 0

        for batch in pf.iter_batches(batch_size=100_000):
            chunk = batch.to_pandas()
            chunk.columns = chunk.columns.str.replace(".", "_", regex=False)
            for _, row in chunk.iterrows():
                pid = row["Precursor_Id"]
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

    _psm_stream_to_specset(psm_iter(), out_dir, max_no_files_per_directory)
    
    print("\nDone.")



def psmlist_to_specset(psmlist_file, out_dir, max_files_in_directory=5000):
    """Read a list of files containing stripped MaxQuant PSMs and write to the specset directory format."""
    '''
    run
    . venv_specset/bin/activate
    python
    #import pandas as pd
    
    from specset.utils import psmlist_to_specset
    psmlist_file = "/d0/tp/sitenet/data2/_files-phospho-260413.tsv"
    out_dir = "/d0/tp/sitenet/data2/phosphopeptides.specset"
    psmlist_to_specset(psmlist_file, out_dir)
    print("Done")
    result = 20313563
    Done.
    '''
    import pandas as pd

    def psm_iter():
        with open(psmlist_file, "r") as file:
            for line in file:
                fname = line.strip()
                print(fname)
                df = pd.read_csv(fname, sep="\t")
                df.columns = df.columns.str.replace(".", "_", regex=False).str.replace(" ", "_", regex=False)
                for psm in df.itertuples():
                    yield {
                        "pep":         psm.Modified_sequence[1:-1],  # remove leading/trailing _
                        "z":           psm.Charge,
                        "rt":          psm.Retention_time,
                        "im":          0,
                        "score":       psm.Score,
                        "id":          psm.Search_ID,
                        "matches":     psm.Matches,
                        "intensities": psm.Intensities,
                    }

    _psm_stream_to_specset(psm_iter(), out_dir, max_files_in_directory)
    print("Done.")


def _specset_stream_to_string(specset_dir, n=None):
    '''
    first n entries
    '''
    import io
    import pandas as pd

    _cols = ["FragmentLossType", "FragmentNumber", "FragmentType", "FragmentCharge", "RelativeIntensity"]

    string_buffer = io.StringIO()
    first = True
    i = 0

    with open(f"{specset_dir}/_index.tsv") as index_file:
        for line in index_file:
            if n is not None and i >= n:
                break

            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split(maxsplit=1)
            if len(parts) != 2:
                continue

            location = parts[1]

            with open(f"{specset_dir}/{location}/__index.tsv") as index2_file:
                for line2 in index2_file:
                    line2 = line2.rstrip("\n")
                    if not line2:
                        continue
                    file_name = line2.split(maxsplit=1)[0]

                    a = {}
                    with open(f"{specset_dir}/{location}/{file_name}.tsv") as psm_file:
                        for psm_line in psm_file:
                            k, v = psm_line.rstrip("\n").split(maxsplit=1)
                            a[k] = v

                    _cols = ["FragmentLossType", "FragmentNumber", "FragmentType", "FragmentCharge", "RelativeIntensity"]

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

            i += 1

    return string_buffer

def specset_to_alphapeptdeep_string(specset_dir, n=None):
    '''
    first n entries
    '''
    import io
    import pandas as pd

    _cols = ["FragmentLossType", "FragmentNumber", "FragmentType", "FragmentCharge", "RelativeIntensity"]

    
    first = True
    i = 0

    with open(f"{specset_dir}/_index.tsv") as index_file:
        for line in index_file:
            if n is not None and i >= n:
                break

            line = line.rstrip("\n")
            if not line:
                continue            

            parts = line.split(maxsplit=1)
            if len(parts) != 2:
                continue

            location = parts[1]

            with open(f"{specset_dir}/{location}/__index.tsv") as index2_file:
                for line2 in index2_file:
                    line2 = line2.rstrip("\n")
                    if not line2:
                        continue
                    file_name = line2.split(maxsplit=1)[0]

                    a = {}
                    with open(f"{specset_dir}/{location}/{file_name}.tsv") as psm_file:
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
                    frag_df.to_csv(string_buffer, sep="\t", index=False, header=first)
                    #first = False
                    yield string_buffer

            i += 1

    #return string_buffer

def specset_to_speclib_alphapep(specset_dir, lib_name, scratch=None):
    """Convert a specset directory to a spectral library TSV (equivalent of psm_to_speclib in psm.py)."""
    """
    from specset.utils import specset_to_speclib
    specset_to_speclib("/home/tpham/data/ai_storage/data/tmp/out.specset", "/home/tpham/data/ai_storage/data/tmp/out.speclib.tsv")
    """

    from alphabase.spectral_library.reader import LibraryReaderBase
    from alphabase.peptide.fragment import create_fragment_mz_dataframe
    from alphabase.spectral_library.base import SpecLibBase
    from peptdeep.spec_lib.translate import translate_to_tsv, mod_to_unimod_dict

    string_buffer = _specset_to_string(specset_dir, scratch)

    lib = LibraryReaderBase()
    psm_df = lib.import_file(string_buffer)

    frag_types = lib.fragment_intensity_df.columns.to_list()
    fragment_mz_df = create_fragment_mz_dataframe(psm_df, frag_types)

    speclib = SpecLibBase()
    speclib.precursor_df = psm_df
    speclib._fragment_intensity_df = lib.fragment_intensity_df
    speclib._fragment_mz_df = fragment_mz_df

    if "decoy" in speclib.precursor_df.columns:
        speclib._precursor_df = speclib.precursor_df[speclib.precursor_df.decoy == 0]

    translate_to_tsv(
        speclib,
        lib_name,
        keep_k_highest_fragments=100000,
        min_frag_nAA=0,
        min_frag_mz=200,
        min_frag_intensity=0.001,
        batch_size=100000,
        translate_mod_dict=mod_to_unimod_dict,
        multiprocessing=False,
    )


def _read_location(location_dir, extra=[]):
    """Read all PSM files in a specset location directory.
    Returns a dict of list of fragment_df, one per charge in __index.tsv.
    """
    import io
    import pandas as pd

    result = {}
    with open(f"{location_dir}/__index.tsv") as index2_file:
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
            with open(f"{location_dir}/{file_name}.tsv") as psm_file:
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


def _consensus(frag_dfs, consensus_method=None):
    """Select or aggregate fragment DataFrames.

    consensus_method : None or str
        None or unrecognised string: select by most fragments.
        'mean': average RelativeIntensity across PSMs per fragment.
        'median': median RelativeIntensity across PSMs per fragment.
        any other string: select by max value of that column.
    """
    import pandas as pd

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


def _to_mz_seq(pep):
    """Convert UniMod peptide notation to mz.py single-char encoding."""
    import re
    seq = pep.replace("(UniMod:1)", "*")
    seq = re.sub(r'M\(UniMod:35\)', '1', seq)
    seq = re.sub(r'S\(UniMod:21\)', '2', seq)
    seq = re.sub(r'T\(UniMod:21\)', '3', seq)
    seq = re.sub(r'Y\(UniMod:21\)', '4', seq)
    seq = re.sub(r'C\(UniMod:4\)', 'C', seq)
    return seq


def specset_to_speclib(specset_dir, lib_name, additional_columns={}, consensus_method=None, n=None):
    """
    additional_columns : dict
        Keys are specset field names, values are the column names to use in the output library.
        Example: {"proteins": "Protein.Ids"}
    """
    """Convert a specset directory to a spectral library TSV """
    """
    from specset.utils import specset_to_speclib
    specset_to_speclib("/home/tpham/data/ai_storage/data/tmp/out.specset", "/home/tpham/data/ai_storage/data/tmp/out.speclib.tsv")
    """

    '''
    first n entries
    '''
    import io
    import pandas as pd
    from .mz import get_ion_mz, get_precursor_mz, generate_aa_mass

    # Calculate the amino acid mass values upon import of the library
    # as they will be used frequently
    aa_mass = generate_aa_mass()

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

    index_path = f"{specset_dir}/_index.tsv"
    with open(index_path) as f:
        total = sum(1 for line in f if line.strip())
    print(f"Total peptides: {total}")

    n_peps = 0

    with open(lib_name, "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")

        with open(index_path) as index_file:
            for line in index_file:
                if n is not None and n_peps >= n:
                    break

                line = line.rstrip("\n")
                if not line:
                    continue

                parts = line.split(maxsplit=1)

                if len(parts) != 2:
                    continue

                location = parts[1]

                for z, frag_dfs in _read_location(f"{specset_dir}/{location}", extra=list(additional_columns.keys())).items():

                    df = _consensus(frag_dfs, consensus_method)
                    pep = _to_mz_seq(df.iloc[0]["ModifiedPeptide"])
                    if not all(c in aa_mass for c in pep):
                        continue
                    precursor_mz = get_precursor_mz(pep, int(z), aa_mass)

                    for row in df.itertuples(index=False):
                        frag_mz = get_ion_mz(pep, row.FragmentType, row.FragmentNumber, int(row.FragmentCharge), row.FragmentLossType, aa_mass)
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


def _specset_to_string(specset_dir, scratch=None):
    import io
    import pandas as pd
    #from tqdm import tqdm

    index_df = pd.read_csv(specset_dir + "/_index.tsv", header=None, names=['peptide', 'location'], sep="\t")

    if scratch is None:
        scratch = len(index_df)

    display = max(scratch // 20, 1)
    string_buffer = io.StringIO()

    #for i in tqdm(range(scratch), desc="Reading specset"):
    for i in range(scratch):
        if i % display == 0:
            print(f"{int((i+1) * 100 / scratch)}%")        

        index2 = pd.read_csv(f"{specset_dir}/{index_df.iloc[i]['location']}/__index.tsv", sep = "\t", header = None, names=['file', 'charge'])
        
        for j in range(len(index2)):
            
            a = pd.read_csv(f"{specset_dir}/{index_df.iloc[i]['location']}/{index2.iloc[j]['file']}.tsv", 
                                  header=None, sep="\t", index_col=0)
            
            frag_df = pd.read_json(io.StringIO(a.loc["ms2",1]), orient='values')
            frag_df.columns = ["FragmentLossType", "FragmentNumber", "FragmentType", "FragmentCharge", "RelativeIntensity"]
            frag_df["ModifiedPeptide"] = a.loc["pep", 1]
            frag_df["PrecursorCharge"] = int(a.loc["z", 1])
            frag_df["RT"] = float(a.loc["ret", 1])
            #frag_df["IM"] = float(a.loc["ccs", 1])
            frag_df["CCS"] = float(a.loc["ccs", 1])

            frag_df.to_csv(string_buffer, sep="\t", index=False, header=(i == 0 and j == 0))

    return string_buffer



def specset_to_speclib_alphabase(specset_dir, lib_name, scratch=None):
    """Convert a specset directory to a spectral library TSV (equivalent of psm_to_speclib in psm.py)."""
    from alphabase.spectral_library.reader import LibraryReaderBase
    from alphabase.peptide.fragment import create_fragment_mz_dataframe
    from alphabase.spectral_library.base import SpecLibBase
    from peptdeep.spec_lib.translate import translate_to_tsv, mod_to_unimod_dict

    string_buffer = _specset_to_string(specset_dir, scratch)

    lib = LibraryReaderBase()
    psm_df = lib.import_file(string_buffer)

    #print(psm_df.columns)
    #if "RT" not in psm_df.columns:
    #    psm_df["RT"] = 0.0
    #if "mobility" not in psm_df.columns:
    #    psm_df["mobility"] = psm_df["IM"].astype(float) if "IM" in psm_df.columns else 0

    frag_types = lib.fragment_intensity_df.columns.to_list()
    fragment_mz_df = create_fragment_mz_dataframe(psm_df, frag_types)

    speclib = SpecLibBase()
    speclib.precursor_df = psm_df
    speclib._fragment_intensity_df = lib.fragment_intensity_df
    speclib._fragment_mz_df = fragment_mz_df

    if "decoy" in speclib.precursor_df.columns:
        speclib._precursor_df = speclib.precursor_df[speclib.precursor_df.decoy == 0]

    translate_to_tsv(
        speclib,
        lib_name,
        keep_k_highest_fragments=100000,
        min_frag_nAA=0,
        min_frag_mz=200,
        min_frag_intensity=0.001,
        batch_size=100000,
        translate_mod_dict=mod_to_unimod_dict,
        multiprocessing=False,
    )
