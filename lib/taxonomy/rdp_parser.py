"""
RDP classifier output parser for MetaWorks pipeline
"""

import sys
from pathlib import Path
from typing import Dict, List, Optional, Union


class RDPParser:
    """Parse RDP classifier output files"""

    def __init__(self, filepath: Union[str, Path]):
        """
        Initialize RDP parser.

        Args:
            filepath: Path to RDP output file
        """
        self.filepath = Path(filepath)
        self.records: List[Dict] = []

    def parse(self) -> List[Dict]:
        """
        Parse RDP output file.

        Returns:
            List of taxonomic assignment records
        """
        if not self.filepath.exists():
            raise FileNotFoundError(f"RDP file not found: {self.filepath}")

        with open(self.filepath, "r") as f:
            for line in f:
                if line.strip():
                    record = self._parse_line(line.strip())
                    if record:
                        self.records.append(record)

        return self.records

    def _parse_line(self, line: str) -> Optional[Dict]:
        """
        Parse a single line of RDP output.

        RDP format: ID\tRoot\trank\tBP\tKingdom\trank\tBP\t...

        Args:
            line: Line from RDP output

        Returns:
            Parsed record as dictionary
        """
        parts = line.split("\t")

        if len(parts) < 3:
            return None

        record = {"id": parts[0], "taxonomy": {}}

        # Parse taxonomic ranks (triplets of: name, rank, bootstrap)
        i = 1
        while i + 2 < len(parts):
            taxon_name = parts[i]
            rank = parts[i + 1]
            bootstrap = parts[i + 2]

            record["taxonomy"][rank.lower()] = {
                "name": taxon_name,
                "rank": rank,
                "bootstrap": float(bootstrap) if bootstrap.replace(".", "").isdigit() else 0.0,
            }

            i += 3

        return record

    def filter_by_bootstrap(self, min_bootstrap: float = 0.5) -> List[Dict]:
        """
        Filter records by minimum bootstrap value.

        Args:
            min_bootstrap: Minimum bootstrap threshold

        Returns:
            Filtered records
        """
        filtered = []

        for record in self.records:
            filtered_taxonomy = {}

            for rank, data in record["taxonomy"].items():
                if data["bootstrap"] >= min_bootstrap:
                    filtered_taxonomy[rank] = data

            if filtered_taxonomy:
                filtered.append({"id": record["id"], "taxonomy": filtered_taxonomy})

        return filtered

    def to_csv(self, ranks: Optional[List[str]] = None) -> str:
        """
        Convert to CSV format.

        Args:
            ranks: List of ranks to include (if None, includes all)

        Returns:
            CSV formatted string
        """
        if not self.records:
            return ""

        # Determine ranks if not provided
        if ranks is None:
            all_ranks = set()
            for record in self.records:
                all_ranks.update(record["taxonomy"].keys())
            ranks = sorted(all_ranks)

        # Build CSV
        lines = []

        # Header
        header = ["ID"]
        for rank in ranks:
            header.extend([f"{rank}_name", f"{rank}_bootstrap"])
        lines.append(",".join(header))

        # Data rows
        for record in self.records:
            row = [record["id"]]

            for rank in ranks:
                if rank in record["taxonomy"]:
                    tax_data = record["taxonomy"][rank]
                    row.append(tax_data["name"])
                    row.append(str(tax_data["bootstrap"]))
                else:
                    row.extend(["", ""])

            lines.append(",".join(row))

        return "\n".join(lines)


def parse_rdp_output(
    filepath: Union[str, Path], min_bootstrap: float = 0.0, output_format: str = "dict"
) -> Union[List[Dict], str]:
    """
    Parse RDP output file.

    Args:
        filepath: Path to RDP output file
        min_bootstrap: Minimum bootstrap threshold
        output_format: Output format ('dict' or 'csv')

    Returns:
        Parsed data in requested format
    """
    parser = RDPParser(filepath)
    records = parser.parse()

    if min_bootstrap > 0:
        records = parser.filter_by_bootstrap(min_bootstrap)

    if output_format == "csv":
        parser.records = records
        return parser.to_csv()

    return records


def main():
    """CLI interface for RDP parser"""
    if len(sys.argv) < 2:
        sys.exit("Usage: python -m lib.taxonomy.rdp_parser <rdp_file> [min_bootstrap] [format]")

    filepath = sys.argv[1]
    min_bootstrap = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    output_format = sys.argv[3] if len(sys.argv) > 3 else "csv"

    try:
        result = parse_rdp_output(filepath, min_bootstrap, output_format)

        if output_format == "csv":
            print(result)
        else:
            import json

            print(json.dumps(result, indent=2))
    except Exception as e:
        sys.exit(f"Error: {e}")


if __name__ == "__main__":
    main()
