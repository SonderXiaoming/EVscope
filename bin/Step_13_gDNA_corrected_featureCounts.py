#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: MIT
"""Apply a simple strand-opposite gDNA subtraction to featureCounts output.

The EVscope shell step now generates two featureCounts tables from the same BAM:
  - s=1 (forward/F1R2-compatible counts)
  - s=2 (reverse/F2R1-compatible counts)
For a forward-stranded library, sense RNA signal is represented by s=1 and the
opposite-strand signal is used as a gDNA proxy.  For a reverse-stranded library,
sense RNA signal is represented by s=2 and s=1 is the proxy.
"""

import argparse
from pathlib import Path

import pandas as pd


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Correct featureCounts read counts for strand-opposite gDNA signal."
    )
    parser.add_argument(
        "--strand",
        required=True,
        choices=["forward", "reverse"],
        help="Strand specificity of the RNA-seq library.",
    )
    parser.add_argument(
        "--forward_featureCounts_table",
        required=True,
        help="featureCounts output generated with -s 1.",
    )
    parser.add_argument(
        "--reverse_featureCounts_table",
        required=True,
        help="featureCounts output generated with -s 2.",
    )
    parser.add_argument("--output", required=True, help="Output corrected featureCounts table.")
    return parser.parse_args()


def read_featurecounts(file_path: str) -> pd.DataFrame:
    """Read a featureCounts table while ignoring comment lines."""
    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"featureCounts table not found: {file_path}")
    df = pd.read_csv(path, sep="\t", comment="#", header=0)
    if df.shape[1] < 7:
        raise ValueError(f"featureCounts table has too few columns: {file_path}")
    return df


def validate_compatible_tables(forward_df: pd.DataFrame, reverse_df: pd.DataFrame) -> None:
    """Fail early if the two featureCounts tables do not describe the same genes."""
    key_cols = list(forward_df.columns[:6])
    if list(reverse_df.columns[:6]) != key_cols:
        raise ValueError("featureCounts annotation columns differ between -s 1 and -s 2 tables")
    if len(forward_df) != len(reverse_df):
        raise ValueError("featureCounts tables have different numbers of rows")
    mismatch = (forward_df.iloc[:, :6].astype(str).values != reverse_df.iloc[:, :6].astype(str).values).any()
    if mismatch:
        raise ValueError("featureCounts tables are not row-aligned; refusing to subtract counts")


def correct_gdna_contamination(forward_df: pd.DataFrame, reverse_df: pd.DataFrame, strand: str) -> pd.DataFrame:
    """Calculate gDNA-corrected counts based on library strandedness."""
    validate_compatible_tables(forward_df, reverse_df)
    count_column = forward_df.columns[-1]
    forward_counts = pd.to_numeric(forward_df.iloc[:, -1], errors="raise")
    reverse_counts = pd.to_numeric(reverse_df.iloc[:, -1], errors="raise")

    corrected_df = forward_df.copy()
    if strand == "forward":
        corrected_counts = forward_counts - reverse_counts
    else:  # strand == "reverse"
        corrected_counts = reverse_counts - forward_counts

    corrected_df[count_column] = corrected_counts.clip(lower=0)
    corrected_df.rename(columns={count_column: "gDNA_corrected_counts"}, inplace=True)
    return corrected_df


def main() -> None:
    args = parse_arguments()
    forward_counts_df = read_featurecounts(args.forward_featureCounts_table)
    reverse_counts_df = read_featurecounts(args.reverse_featureCounts_table)
    corrected_counts_df = correct_gdna_contamination(forward_counts_df, reverse_counts_df, args.strand)
    corrected_counts_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
