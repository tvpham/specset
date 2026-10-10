"""Command line interfaces"""

import argparse


def _key_value(s):
    """Parse a KEY=VALUE command line argument into a (key, value) tuple."""
    if "=" not in s:
        raise argparse.ArgumentTypeError(f"expected KEY=VALUE, got '{s}'")
    key, value = s.split("=", 1)
    return key, value


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


def _cmd_diannreport_to_specset(args):
    from .utils import diannreport_to_specset

    if args.column and not args.values:
        raise SystemExit("diannreport_to_specset: --values is required with --column")
    diannreport_to_specset(args.diannreport, args.out_dir, column=args.column or "", values=args.values,
                           max_no_files_per_directory=args.max_files_in_directory,
                           overwrite=args.overwrite, by_Run=args.by_Run)


def _cmd_psmlist_to_specset(args):
    from .utils import psmlist_to_specset

    psmlist_to_specset(args.psmlist_file, args.out_dir, max_files_in_directory=args.max_files_in_directory,
                       overwrite=args.overwrite)


def _cmd_specset_to_speclib(args):
    from .utils import specset_to_speclib

    specset_to_speclib(args.specset_dir, args.lib_name, additional_columns=dict(args.additional_column),
                       consensus_method=args.consensus_method, n=args.n)


def _cmd_merge_specsets(args):
    from .utils import merge_specsets

    merge_specsets(args.specset_a, args.specset_b, args.specset_out, max_files_in_directory=args.max_files_in_directory,
                   mode=args.mode, flush_threshold=args.flush_threshold, overwrite=args.overwrite)


def _cmd_summarize(args):
    from .utils import summarize

    summarize(args.specset_a, args.specset_b)


def _add_max_files(p):
    p.add_argument("--max-files-in-directory", type=int, default=5000, help="Number of PSM files per batch folder (default: 5000)")


def _add_overwrite(p):
    p.add_argument("--overwrite", action="store_true", help="Delete the output directory first if it exists and is not empty (default: False)")


def _add_merge_options(p):
    _add_max_files(p)
    p.add_argument("--mode", choices=["both", "ignore"], default="ignore",
                   help="For peptides in both: 'both' keeps the PSMs from both, 'ignore' skips the ones from specset_b (default: ignore)")
    p.add_argument("--flush-threshold", type=int, default=1000, help="Number of PSMs to buffer before flushing the index to disk (default: 1000)")


def main():

    parser = argparse.ArgumentParser(prog="specset", description="Toolkit for manipulation of a set of mass spectra")
    subparsers = parser.add_subparsers(dest="command", required=True)

    p = subparsers.add_parser("maxquant_to_specset", help="Convert a MaxQuant txt directory (msms.txt and evidence.txt) to a specset")
    p.add_argument("txt_dir", help="MaxQuant txt directory containing msms.txt and evidence.txt")
    p.add_argument("out_dir", help="Output specset directory")
    p.add_argument("--additional-column", action="append", default=[], type=_key_value, metavar="COLUMN=NAME",
                   help="msms.txt column to keep and the name to use in the specset, e.g. 'Gene Names=genes' (repeatable)")
    _add_max_files(p)
    p.add_argument("-n", type=int, default=-1, help="Only convert the first N PSMs (default: all)")
    _add_overwrite(p)
    p.set_defaults(func=_cmd_maxquant_to_specset)

    p = subparsers.add_parser("easypqp_to_specset", help="Convert an EasyPQP library TSV to a specset")
    p.add_argument("easypqp_file", help="EasyPQP library TSV file")
    p.add_argument("out_dir", help="Output specset directory")
    _add_max_files(p)
    _add_overwrite(p)
    p.set_defaults(func=_cmd_easypqp_to_specset)

    p = subparsers.add_parser("diannlib_to_specset", help="Convert a DIA-NN parquet library to a specset")
    p.add_argument("diannlib_file", help="DIA-NN library parquet file")
    p.add_argument("out_dir", help="Output specset directory")
    _add_max_files(p)
    _add_overwrite(p)
    p.set_defaults(func=_cmd_diannlib_to_specset)

    p = subparsers.add_parser("diannreport_to_specset", help="Convert a DIA-NN parquet report to a specset, or one specset per run")
    p.add_argument("diannreport", help="DIA-NN report parquet file")
    p.add_argument("out_dir", help="Output specset directory; with --by-Run, the directory in which one <Run>.specset is created per run")
    _add_max_files(p)
    p.add_argument("--by-Run", dest="by_Run", action="store_true", help="Create one specset per run (default: one specset with all runs)")
    p.add_argument("--column", default="", help="Report column to select rows on, e.g. 'Run' (default: import all rows)")
    p.add_argument("--values", nargs="+", default=[], metavar="VALUE",
                   help="With --column, import only the rows whose value in that column is one of these, e.g. --values plasma1 plasma2")
    p.add_argument("--overwrite", action="store_true",
                   help="Delete an output specset first if it exists and is not empty (default: fail)")
    p.set_defaults(func=_cmd_diannreport_to_specset)

    p = subparsers.add_parser("psmlist_to_specset", help="Convert a list of files containing stripped MaxQuant PSMs to a specset")
    p.add_argument("psmlist_file", help="Text file with one PSM file path per line")
    p.add_argument("out_dir", help="Output specset directory")
    _add_max_files(p)
    _add_overwrite(p)
    p.set_defaults(func=_cmd_psmlist_to_specset)

    p = subparsers.add_parser("specset_to_speclib", help="Convert a specset to a spectral library TSV")
    p.add_argument("specset_dir", help="Specset directory")
    p.add_argument("lib_name", help="Output spectral library TSV file")
    p.add_argument("--additional-column", action="append", default=[], type=_key_value, metavar="FIELD=COLUMN",
                   help="Specset field to keep and the column name to use in the library, e.g. 'proteins=Protein.Ids' (repeatable)")
    p.add_argument("--consensus-method", default=None,
                   help="'mean', 'median', or a specset field to select the PSM with the max value (default: PSM with most fragments)")
    p.add_argument("-n", type=int, default=None, help="Only use the first N peptides (default: all)")
    p.set_defaults(func=_cmd_specset_to_speclib)

    p = subparsers.add_parser("merge_specsets", help="Merge specset_a and specset_b into a new specset")
    p.add_argument("specset_a", help="Primary specset directory (left unchanged)")
    p.add_argument("specset_b", help="Secondary specset directory (left unchanged)")
    p.add_argument("specset_out", help="Output specset directory")
    _add_merge_options(p)
    _add_overwrite(p)
    p.set_defaults(func=_cmd_merge_specsets)

    p = subparsers.add_parser("summarize", help="Print summary statistics for one or two specsets")
    p.add_argument("specset_a", help="Specset directory")
    p.add_argument("specset_b", nargs="?", default=None, help="Second specset directory; also reports peptide overlap and differences")
    p.set_defaults(func=_cmd_summarize)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":

    main()
