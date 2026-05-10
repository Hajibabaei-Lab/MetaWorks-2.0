#!/usr/bin/env python3
import argparse
import csv
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Map global IDs (ESV or OTU) to trial sequences in results CSV using VSEARCH .uc clustering output."
    )
    parser.add_argument("uc_file", help="VSEARCH .uc clustering file")
    parser.add_argument("results_csv", help="Results CSV file")
    parser.add_argument(
        "--column-name",
        default="GlobalESV",
        help="Column name to search for and append (default: GlobalESV)",
    )
    args = parser.parse_args(argv)

    id_map = {}

    try:
        with open(args.uc_file, "r") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) < 10:
                    continue
                record_type = fields[0]
                query_id = fields[8]
                if record_type == "H":
                    id_map[query_id] = fields[9]
                elif record_type == "N":
                    id_map[query_id] = "NoGlobalMatch"
    except FileNotFoundError:
        print(f"UC file not found: {args.uc_file}", file=sys.stderr)
        sys.exit(1)

    try:
        with open(args.results_csv, "r", newline="") as f:
            reader = csv.reader(f)
            header = next(reader, None)
            if header is None:
                return

            col_name = args.column_name
            for i, col in enumerate(header):
                if col == col_name:
                    header[i] = "TrialESV"
                    break
            header.append(col_name)

            writer = csv.writer(sys.stdout)
            writer.writerow(header)

            for row in reader:
                trial_id = row[0] if row else ""
                global_id = id_map.get(trial_id, "NoGlobalMatch")
                row.append(global_id)
                writer.writerow(row)
    except FileNotFoundError:
        print(f"Results CSV not found: {args.results_csv}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
