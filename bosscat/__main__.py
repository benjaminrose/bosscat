"""command-line interface."""

import argparse
import sys

from .classifiers.snid import make_snid
from .config import SnidConfig
from .pipeline import run_pipeline

CLASSIFIER_FACTORIES = {"snid": make_snid}


def _rebin_factor(text):
    """Parse ``--rebin-factor``: an integer >= 1."""
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"invalid integer: {text!r}") from None
    if value < 1:
        raise argparse.ArgumentTypeError(f"must be >= 1, got {value}")
    return value


def parse_args():
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments with attributes ``input_files``, ``classifier``,
        ``no_classifier``, ``zmax``, ``rebin_factor``, ``templates_dir`` and
        ``output``. Unspecified SNID options are None.
    """
    parser = argparse.ArgumentParser(
        description="BOSS Classification Analysis Toolkit",
    )
    parser.add_argument(
        "input_files",
        nargs="+",
        help="BOSS coadd FITS file path(s)",
    )
    parser.add_argument(
        "--classifier",
        action="append",
        choices=sorted(CLASSIFIER_FACTORIES),
        help="Classifier to run; repeat for several (default: snid)",
    )
    parser.add_argument(
        "--no-classifier",
        action="store_true",
        help="Run no classifier (metadata only)",
    )
    parser.add_argument(
        "--zmax", type=float,
        help="SNID maximum redshift (default: SnidConfig default)",
    )
    parser.add_argument(
        "--rebin-factor", type=_rebin_factor,
        help="SNID rebinning factor (default: SnidConfig default)",
    )
    parser.add_argument(
        "--templates-dir",
        help="SNID templates directory (default: $BOSSCAT_SNID_TEMPLATES)",
    )
    parser.add_argument(
        "-o", "--output",
        default="results.csv",
        help="Output CSV file path (default: %(default)s)",
    )
    args = parser.parse_args()
    if args.classifier is None:
        args.classifier = ["snid"]
    return args


def main():
    """Entry point for the ``bosscat`` command-line script."""
    args = parse_args()
    overrides = {
        key: getattr(args, key)
        for key in ("zmax", "rebin_factor", "templates_dir")
        if getattr(args, key) is not None
    }
    cfg = SnidConfig(**overrides)
    names = [] if args.no_classifier else list(dict.fromkeys(args.classifier))
    classifiers = [CLASSIFIER_FACTORIES[n](cfg) for n in names]
    try:
        run_pipeline(args.input_files, args.output, classifiers)
    except RuntimeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
