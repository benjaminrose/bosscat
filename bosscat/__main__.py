"""command-line interface."""

import argparse

from .pipeline import run_pipeline


def parse_args():
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments with attributes ``input_files``, ``snid``,
        and ``output``.
    """
    parser = argparse.ArgumentParser(
        description="BOSS Classification Analysis Toolkit",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "input_files",
        nargs="+",
        help="BOSS coadd FITS file path(s)",
    )
    parser.add_argument(
        "--snid",
        default=True,
        action=argparse.BooleanOptionalAction,
        help="Run the SNID-SAGE classifier",
    )
    parser.add_argument(
        "-o", "--output",
        default="results.csv",
        help="Output CSV file path",
    )
    return parser.parse_args()


def main():
    """Entry point for the ``bosscat`` command-line script."""
    args = parse_args()
    run_pipeline(args.input_files, args.output, run_snid=args.snid)


if __name__ == "__main__":
    main()
