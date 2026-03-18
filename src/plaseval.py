#!/usr/bin/env python

"""PlasEval (Plasmids Evaluator), a method for evaluating the predictions from plasmid binning tools.

Two modes of PlasEval:
- eval: evaluates plasmid bins against a set of ground truth bins to provide precision-recall statistics
- compare: compares two sets of plasmid bins to quantify the dissimilarity between the two given sets
"""  # noqa: E501

import argparse
from pathlib import Path

import evaluate_bins as eb
import plasmid_comparison_main as pcm


def main() -> None:
    """Compute the main CLI entry point."""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="mode to be used", dest="mode")
    subparsers.required = True
    # Evaluate mode
    eval_parser = subparsers.add_parser("eval", help="evaluate precision and recall")
    eval_parser.add_argument("--pred", help="Path to predictions")
    eval_parser.add_argument("--gt", help="Path to contig to plasmid mapping file")
    eval_parser.add_argument(
        "--min_len",
        type=int,
        default=0,
        help="Minimum length of contigs",
    )
    eval_parser.add_argument("--out_file", help="Path to output file")
    eval_parser.add_argument("--log_file", help="Path to log file")
    # Compare mode
    comp_parser = subparsers.add_parser("comp", help="compare two sets of plasmid bins")
    comp_parser.add_argument("--l", help="Path to file with 1st set of plasmids")
    comp_parser.add_argument("--r", help="Path to file with 2nd set of plasmids")
    comp_parser.add_argument("--p", type=float, default=0.5, help="Weight exponent")
    comp_parser.add_argument(
        "--min_len",
        type=int,
        default=0,
        help="Minimum length of contigs",
    )
    comp_parser.add_argument(
        "--max_calls",
        type=int,
        default=10000000,
        help="Maximum number of recursive function calls",
    )
    comp_parser.add_argument("--out_file", help="Path to output file")
    comp_parser.add_argument("--log_file", help="Path to log file")
    args = parser.parse_args()

    if args.mode == "eval":
        eb.eval_mode(
            args.pred,
            args.gt,
            args.min_len,
            Path(args.out_file),
            Path(args.log_file),
        )
    if args.mode == "comp":
        pcm.comp_mode(
            args.l,
            args.r,
            args.p,
            args.min_len,
            args.max_calls,
            Path(args.out_file),
            Path(args.log_file),
        )


if __name__ == "__main__":
    main()
