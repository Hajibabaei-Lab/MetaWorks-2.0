#!/usr/bin/env python3
import argparse
import gzip
import sys


def _open_input(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Prepare per-sample reads for pooled denoising: prepend the sample "
            "name to every FASTA header and write gzipped output to stdout."
        )
    )
    parser.add_argument("input", help="Input FASTA file (.fasta or .fasta.gz)")
    parser.add_argument(
        "--sample-name",
        required=True,
        help="Sample name prepended to each header as '{sample}_{original}'.",
    )
    args = parser.parse_args(argv)

    out = sys.stdout.buffer
    try:
        inf = _open_input(args.input)
    except OSError as e:
        sys.exit(f"Cannot open infile: {e}")
    with inf, gzip.GzipFile(fileobj=out, mode="wb", mtime=0) as gz:
        for line in inf:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                line = f">{args.sample_name}_{line[1:]}".replace("-", "_")
            gz.write((line + "\n").encode())


if __name__ == "__main__":
    main()
