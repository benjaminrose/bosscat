"""command-line interface."""

import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="BOSS Classification Analysis Toolkit",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        # suggest_on_error=True,
    )
    # Accept one or more file paths as positional arguments
    parser.add_argument(
        "input_files",
        nargs="+",
        help="BOSS spectra file path",
    )
    parser.add_argument(
        "--snid",
        default=True,
        action=argparse.BooleanOptionalAction,
        help="",
    )
    # parser.add_argument(
    #     "--fits",
    #     default=True,
    #     action=argparse.BooleanOptionalAction,
    #     help="Save data to fits",
    # )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
