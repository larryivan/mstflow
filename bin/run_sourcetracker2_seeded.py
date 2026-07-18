#!/usr/bin/env python3
"""Run SourceTracker2 2.0.1 with deterministic single-process randomness."""

import argparse
import glob
import os
import random

import numpy as np


def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--seed", type=int, required=True)
    args, sourcetracker_args = parser.parse_known_args()
    if not 0 <= args.seed <= 2_147_483_647:
        parser.error("--seed must be between 0 and 2147483647")
    if not sourcetracker_args or sourcetracker_args[0] != "gibbs":
        parser.error("the wrapped SourceTracker2 command must begin with 'gibbs'")

    expected_hash_seed = str(args.seed)
    if os.environ.get("PYTHONHASHSEED") != expected_hash_seed:
        parser.error(
            "PYTHONHASHSEED must be exported to the same value as --seed before Python starts"
        )

    np.random.seed(args.seed)
    random.seed(args.seed)

    # SourceTracker2 collates per-sample files from glob order. The standard
    # library does not promise that order, so sort it for stable output rows.
    original_glob = glob.glob
    glob.glob = lambda *glob_args, **glob_kwargs: sorted(
        original_glob(*glob_args, **glob_kwargs)
    )

    from sourcetracker._cli import cli

    cli.main(
        args=sourcetracker_args,
        prog_name="sourcetracker2",
        standalone_mode=True,
    )


if __name__ == "__main__":
    main()
