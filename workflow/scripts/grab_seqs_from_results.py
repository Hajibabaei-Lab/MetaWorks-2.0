import argparse
import csv
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Extract FASTA sequences from results.csv"
    )
    parser.add_argument(
        "results_csv",
        help="Path to results.csv file",
    )
    args = parser.parse_args(argv)

    with open(args.results_csv, newline="") as fh:
        reader = csv.reader(fh)
        next(reader, None)
        for row in reader:
            if len(row) < 4:
                print(
                    f"Warning: skipping row with only {len(row)} columns: {row}",
                    file=sys.stderr,
                )
                continue
            print(f">{row[0]}")
            print(row[3])


if __name__ == "__main__":
    main()
