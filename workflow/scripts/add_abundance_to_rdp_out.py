#!/usr/bin/env python3
import argparse
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Add read abundance from ESV table to RDP taxonomy output."
    )
    parser.add_argument("table_file", help="ESV abundance table")
    parser.add_argument("rdp_file", help="RDP taxonomy output file")
    parser.add_argument("header", nargs="?", default=None,
                        help="Optional header to print as first line")
    parser.add_argument("min_abund", nargs="?", type=int, default=3,
                        help="Minimum abundance threshold (default: 3)")
    args = parser.parse_args(argv)

    assignment = {}

    with open(args.rdp_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if not parts:
                continue
            idline = parts.pop(0)
            id_parts = idline.split(";")
            global_otu = id_parts[0]
            assignment_str = ",".join(parts)
            assignment[global_otu] = assignment_str

    table_data = {}
    headers = []

    with open(args.table_file, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        line = line.rstrip("\n")
        if i == 0:
            headers = line.split()
            if len(headers) >= 2:
                headers = headers[2:]
        else:
            parts = line.split()
            if not parts:
                continue
            global_otu = parts.pop(0).strip()
            if len(parts) != len(headers):
                print("Check FAILED", file=sys.stderr)
            sample_abundances = {}
            for j, abund_str in enumerate(parts):
                sample = headers[j]
                try:
                    abund = int(abund_str)
                except ValueError:
                    try:
                        abund = float(abund_str)
                    except ValueError:
                        abund = 0
                sample_abundances[sample] = abund
            table_data[global_otu] = sample_abundances

    if args.header:
        print(args.header)

    for global_otu, tax_assignment in assignment.items():
        if global_otu in table_data:
            sample_abundances = table_data[global_otu]
            for sample, abund in sample_abundances.items():
                if abund >= args.min_abund:
                    print(f"{global_otu},{sample},{abund},{tax_assignment}")
        else:
            print(f"Cannot find global_otu {global_otu} in table", file=sys.stderr)


if __name__ == "__main__":
    main()
