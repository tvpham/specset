import argparse
import sys


def hello_world():
    print("Hello, world")


def _key_value(s):
    """Parse a KEY=VALUE command line argument into a (key, value) tuple."""
    if "=" not in s:
        raise argparse.ArgumentTypeError(f"expected KEY=VALUE, got '{s}'")
    key, value = s.split("=", 1)
    return key, value


def _raw_stdout():
    """Return stdout with newline translation disabled, so it gets the same bytes as an output file."""
    sys.stdout.reconfigure(newline="")
    return sys.stdout


def _write_text(text, output=None):
    """Write text to a file, or to stdout if no file is given."""
    if output is None:
        _raw_stdout().write(text)
    else:
        with open(output, "w", encoding="utf-8", newline="") as f:
            f.write(text)


def _cmd_hello(args):
    hello_world()


def _cmd_iter_specset_locations(args):
    from .utils import iter_specset_locations

    lines = [f"{pep}\t{location}\n" for pep, location in iter_specset_locations(args.specset_dir, n=args.n)]
    _write_text("".join(lines), args.output)


def _cmd_iter_specset_files(args):
    from .utils import iter_specset_files

    lines = [f"{location}\t{file_name}\t{'' if charge is None else charge}\n"
             for location, file_name, charge in iter_specset_files(args.specset_dir, n=args.n)]
    _write_text("".join(lines), args.output)


def _cmd_matches_intensities_to_psm(args):
    from .utils import matches_intensities_to_psm

    df = matches_intensities_to_psm(args.matches, args.intensities, sep=args.sep, max_intensity=args.max_intensity)
    if df is None:
        sys.exit("No valid fragments.")
    _write_text(df.to_csv(sep="\t", index=False, lineterminator="\n"), args.output)


def _cmd_maxquant_to_specset(args):
    from .utils import maxquant_to_specset

    maxquant_to_specset(args.txt_dir, args.out_dir, additional_columns=dict(args.additional_column),
                        max_files_in_directory=args.max_files_in_directory, n=args.n, overwrite=args.overwrite)


def _cmd_easypqp_to_specset(args):
    from .utils import easypqp_to_specset

    easypqp_to_specset(args.easypqp_file, args.out_dir, max_no_files_per_directory=args.max_files_in_directory,
                       overwrite=args.overwrite)


def _cmd_diannlib_to_specset(args):
    from .utils import diannlib_to_specset

    diannlib_to_specset(args.diannlib_file, args.out_dir, max_no_files_per_directory=args.max_files_in_directory,
                        overwrite=args.overwrite)


def _cmd_psmlist_to_specset(args):
    from .utils import psmlist_to_specset

    psmlist_to_specset(args.psmlist_file, args.out_dir, max_files_in_directory=args.max_files_in_directory,
                       overwrite=args.overwrite)


def _cmd_specset_stream_to_string(args):
    from .utils import _specset_stream_to_string

    _write_text(_specset_stream_to_string(args.specset_dir, n=args.n).getvalue(), args.output)


def _cmd_specset_to_alphapeptdeep_string(args):
    from .utils import specset_to_alphapeptdeep_string

    f = _raw_stdout() if args.output is None else open(args.output, "w", encoding="utf-8", newline="")
    try:
        first = True
        for string_buffer in specset_to_alphapeptdeep_string(args.specset_dir, n=args.n):
            text = string_buffer.getvalue()
            # every buffer carries its own header; keep only the first one
            f.write(text if first else text.split("\n", 1)[1])
            first = False
    finally:
        if f is not sys.stdout:
            f.close()


def _cmd_read_location(args):
    import pandas as pd
    from .utils import _read_location

    frag_dfs = [df for dfs in _read_location(args.location_dir, extra=args.extra).values() for df in dfs]
    if not frag_dfs:
        sys.exit("No PSMs found.")
    _write_text(pd.concat(frag_dfs, ignore_index=True).to_csv(sep="\t", index=False, lineterminator="\n"), args.output)


def _cmd_to_mz_seq(args):
    from .utils import _to_mz_seq

    print(_to_mz_seq(args.pep))


def _cmd_specset_to_speclib(args):
    from .utils import specset_to_speclib

    specset_to_speclib(args.specset_dir, args.lib_name, additional_columns=dict(args.additional_column),
                       consensus_method=args.consensus_method, n=args.n)


def _cmd_specset_to_string(args):
    import contextlib
    from .utils import _specset_to_string

    # progress is printed to stdout; keep it out of the TSV
    with contextlib.redirect_stdout(sys.stderr):
        string_buffer = _specset_to_string(args.specset_dir, scratch=args.scratch)
    _write_text(string_buffer.getvalue(), args.output)


def _cmd_merge_specsets(args):
    from .utils import merge_specsets

    merge_specsets(args.specset_a, args.specset_b, args.specset_out, max_files_in_directory=args.max_files_in_directory,
                   mode=args.mode, flush_threshold=args.flush_threshold, overwrite=args.overwrite)


def _cmd_merge_specsets_inplace(args):
    from .utils import merge_specsets_inplace

    merge_specsets_inplace(args.specset_a, args.specset_b, max_files_in_directory=args.max_files_in_directory,
                           mode=args.mode, flush_threshold=args.flush_threshold)


def _cmd_summarize(args):
    from .utils import summarize

    summarize(args.specset_a, args.specset_b)


def _add_max_files(p):
    p.add_argument("--max-files-in-directory", type=int, default=5000, help="Number of PSM files per batch folder (default: 5000)")


def _add_output(p):
    p.add_argument("-o", "--output", default=None, help="Output file (default: stdout)")


def _add_overwrite(p):
    p.add_argument("--overwrite", action="store_true", help="Delete the output directory first if it exists and is not empty (default: fail)")


def _add_merge_options(p):
    _add_max_files(p)
    p.add_argument("--mode", choices=["both", "ignore"], default="both",
                   help="For peptides in both: 'both' keeps the PSMs from both, 'ignore' skips the ones from specset_b (default: both)")
    p.add_argument("--flush-threshold", type=int, default=1000, help="Number of PSMs to buffer before flushing the index to disk (default: 1000)")


def main():
    parser = argparse.ArgumentParser(prog="specset", description="Toolkit for manipulation of a set of mass spectra")
    subparsers = parser.add_subparsers(dest="command", required=True)

    p = subparsers.add_parser("hello", help="Print Hello, world")
    p.set_defaults(func=_cmd_hello)

    p = subparsers.add_parser("iter-specset-locations", help="List the (peptide, location) pairs of a specset")
    p.add_argument("specset_dir", help="Specset directory")
    p.add_argument("-n", type=int, default=None, help="Only list the first N peptides (default: all)")
    _add_output(p)
    p.set_defaults(func=_cmd_iter_specset_locations)

    p = subparsers.add_parser("iter-specset-files", help="List the (location, file name, charge) of every PSM file of a specset")
    p.add_argument("specset_dir", help="Specset directory")
    p.add_argument("-n", type=int, default=None, help="Only use the first N peptides (default: all)")
    _add_output(p)
    p.set_defaults(func=_cmd_iter_specset_files)

    p = subparsers.add_parser("matches-intensities-to-psm",help="Parse MaxQuant fragment matches and intensities into a fragment table")
    p.add_argument("matches", help="MaxQuant 'Matches' string, e.g. 'b3;y3-NH3;b4'")
    p.add_argument("intensities", help="MaxQuant 'Intensities' string, e.g. '64.1;7.8;53.3'")
    p.add_argument("--sep", default=";", help="Separator used in matches and intensities (default: ';')")
    p.add_argument("--max-intensity", type=float, default=-1, help="Rescale intensities so that the maximum is this value (default: no rescaling)")
    _add_output(p)
    p.set_defaults(func=_cmd_matches_intensities_to_psm)

    p = subparsers.add_parser("maxquant-to-specset", help="Convert a MaxQuant txt directory (msms.txt and evidence.txt) to a specset")
    p.add_argument("txt_dir", help="MaxQuant txt directory containing msms.txt and evidence.txt")
    p.add_argument("out_dir", help="Output specset directory")
    p.add_argument("--additional-column", action="append", default=[], type=_key_value, metavar="COLUMN=NAME",
                   help="msms.txt column to keep and the name to use in the specset, e.g. 'Gene Names=genes' (repeatable)")
    _add_max_files(p)
    p.add_argument("-n", type=int, default=-1, help="Only convert the first N PSMs (default: all)")
    _add_overwrite(p)
    p.set_defaults(func=_cmd_maxquant_to_specset)

    p = subparsers.add_parser("easypqp-to-specset", help="Convert an EasyPQP library TSV to a specset")
    p.add_argument("easypqp_file", help="EasyPQP library TSV file")
    p.add_argument("out_dir", help="Output specset directory")
    _add_max_files(p)
    _add_overwrite(p)
    p.set_defaults(func=_cmd_easypqp_to_specset)

    p = subparsers.add_parser("diannlib-to-specset", help="Convert a DIA-NN parquet library to a specset")
    p.add_argument("diannlib_file", help="DIA-NN library parquet file")
    p.add_argument("out_dir", help="Output specset directory")
    _add_max_files(p)
    _add_overwrite(p)
    p.set_defaults(func=_cmd_diannlib_to_specset)

    p = subparsers.add_parser("psmlist-to-specset", help="Convert a list of files containing stripped MaxQuant PSMs to a specset")
    p.add_argument("psmlist_file", help="Text file with one PSM file path per line")
    p.add_argument("out_dir", help="Output specset directory")
    _add_max_files(p)
    _add_overwrite(p)
    p.set_defaults(func=_cmd_psmlist_to_specset)

    p = subparsers.add_parser("specset-stream-to-string", help="Dump a specset as a fragment TSV")
    p.add_argument("specset_dir", help="Specset directory")
    p.add_argument("-n", type=int, default=None, help="Only use the first N peptides (default: all)")
    _add_output(p)
    p.set_defaults(func=_cmd_specset_stream_to_string)

    p = subparsers.add_parser("specset-to-alphapeptdeep-string", help="Dump a specset as an AlphaPeptDeep fragment TSV")
    p.add_argument("specset_dir", help="Specset directory")
    p.add_argument("-n", type=int, default=None, help="Only use the first N peptides (default: all)")
    _add_output(p)
    p.set_defaults(func=_cmd_specset_to_alphapeptdeep_string)

    p = subparsers.add_parser("read-location", help="Dump all PSMs in a specset location directory as a fragment TSV")
    p.add_argument("location_dir", help="Specset location directory (containing __index.tsv)")
    p.add_argument("--extra", action="append", default=[], metavar="FIELD", help="Additional specset field to include as a column (repeatable)")
    _add_output(p)
    p.set_defaults(func=_cmd_read_location)

    p = subparsers.add_parser("to-mz-seq", help="Convert UniMod peptide notation to mz.py single-char encoding")
    p.add_argument("pep", help="Peptide in UniMod notation, e.g. 'AC(UniMod:4)M(UniMod:35)K'")
    p.set_defaults(func=_cmd_to_mz_seq)

    p = subparsers.add_parser("specset-to-speclib", help="Convert a specset to a spectral library TSV")
    p.add_argument("specset_dir", help="Specset directory")
    p.add_argument("lib_name", help="Output spectral library TSV file")
    p.add_argument("--additional-column", action="append", default=[], type=_key_value, metavar="FIELD=COLUMN",
                   help="Specset field to keep and the column name to use in the library, e.g. 'proteins=Protein.Ids' (repeatable)")
    p.add_argument("--consensus-method", default=None,
                   help="'mean', 'median', or a specset field to select the PSM with the max value (default: PSM with most fragments)")
    p.add_argument("-n", type=int, default=None, help="Only use the first N peptides (default: all)")
    p.set_defaults(func=_cmd_specset_to_speclib)

    p = subparsers.add_parser("specset-to-string", help="Dump a legacy (JSON ms2) specset as a fragment TSV")
    p.add_argument("specset_dir", help="Specset directory")
    p.add_argument("--scratch", type=int, default=None, help="Only use the first N peptides (default: all)")
    _add_output(p)
    p.set_defaults(func=_cmd_specset_to_string)

    p = subparsers.add_parser("merge-specsets", help="Merge specset_a and specset_b into a new specset")
    p.add_argument("specset_a", help="Primary specset directory (left unchanged)")
    p.add_argument("specset_b", help="Secondary specset directory (left unchanged)")
    p.add_argument("specset_out", help="Output specset directory")
    _add_merge_options(p)
    _add_overwrite(p)
    p.set_defaults(func=_cmd_merge_specsets)

    p = subparsers.add_parser("merge-specsets-inplace", help="Merge specset_b into specset_a in place")
    p.add_argument("specset_a", help="Specset directory to merge into (modified in place)")
    p.add_argument("specset_b", help="Specset directory to merge from")
    _add_merge_options(p)
    p.set_defaults(func=_cmd_merge_specsets_inplace)

    p = subparsers.add_parser("summarize", help="Print summary statistics for one or two specsets")
    p.add_argument("specset_a", help="Specset directory")
    p.add_argument("specset_b", nargs="?", default=None, help="Second specset directory; also reports peptide overlap and differences")
    p.set_defaults(func=_cmd_summarize)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
