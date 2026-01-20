#!/usr/bin/env python3
"""
Convert VSEARCH SINTAX tabbed output to MetaWorks' RDP-like `rdp.out.tmp` TSV.

MetaWorks downstream steps (pseudogene filtering, taxonomy.csv creation) expect
`rdp.out.tmp` to have an exact column count depending on the marker.

This converter:
- reads a SINTAX `--tabbedout` file
- extracts `tax=...` assignments when present (or best-effort from other columns)
- maps ranks into the marker-specific RDP-like column layout
- fills missing ranks with empty names and 0.0 confidence
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


@dataclass(frozen=True)
class RankSpec:
    name_col: str
    rank_col: str
    bp_col: str
    rank_label: str
    accepted_codes: Tuple[str, ...]


def marker_rank_specs(marker: str) -> List[RankSpec]:
    # These specs mirror the column layouts used by:
    # - python_scripts/rdp_tsv_to_csv.py
    # - python_scripts/filter_rdp_taxonomy.py
    tax3_abund6 = {
        "COI",
        "rbcL_eukaryota",
        "rbcL_landPlant",
        "12S_fish",
        "12S_vertebrate",
        "16S_vertebrate",
    }
    tax3_abund7 = {"16S", "28S_fungi"}
    tax3_abund9 = {"18S_eukaryota"}

    # Extra marker layout used by filter_rdp_taxonomy.py
    with_subkingdom = {"rbcL_diatom"}
    its_fungi = {"ITS_fungi"}

    marker = marker.strip()

    if marker in tax3_abund6:
        return [
            RankSpec("Root", "RootRank", "rBP", "root", ("root", "r")),
            RankSpec("SuperKingdom", "SuperKingdomRank", "skBP", "superkingdom", ("sk", "superkingdom", "d")),
            RankSpec("Kingdom", "KingdomRank", "kBP", "kingdom", ("k", "kingdom")),
            RankSpec("Phylum", "PhylumRank", "pBP", "phylum", ("p", "phylum")),
            RankSpec("Class", "ClassRank", "cBP", "class", ("c", "class")),
            RankSpec("Order", "OrderRank", "oBP", "order", ("o", "order")),
            RankSpec("Family", "FamilyRank", "fBP", "family", ("f", "family")),
            RankSpec("Genus", "GenusRank", "gBP", "genus", ("g", "genus")),
            RankSpec("Species", "SpeciesRank", "sBP", "species", ("s", "species")),
        ]

    if marker in tax3_abund7:
        return [
            RankSpec("Domain", "DomainRank", "dBP", "domain", ("d", "domain")),
            RankSpec("Phylum", "PhylumRank", "pBP", "phylum", ("p", "phylum")),
            RankSpec("Class", "ClassRank", "cBP", "class", ("c", "class")),
            RankSpec("Order", "OrderRank", "oBP", "order", ("o", "order")),
            RankSpec("Family", "FamilyRank", "fBP", "family", ("f", "family")),
            RankSpec("Genus", "GenusRank", "gBP", "genus", ("g", "genus")),
        ]

    if marker in tax3_abund9:
        return [
            RankSpec("Root", "RootRank", "rBP", "root", ("root", "r")),
            RankSpec("Domain", "DomainRank", "dBP", "domain", ("d", "domain")),
            RankSpec("Kingdom", "KingdomRank", "kBP", "kingdom", ("k", "kingdom")),
            RankSpec("Phylum", "PhylumRank", "pBP", "phylum", ("p", "phylum")),
            RankSpec("Class", "ClassRank", "cBP", "class", ("c", "class")),
            RankSpec("Order", "OrderRank", "oBP", "order", ("o", "order")),
            RankSpec("Family", "FamilyRank", "fBP", "family", ("f", "family")),
            RankSpec("Genus", "GenusRank", "gBP", "genus", ("g", "genus")),
        ]

    if marker in with_subkingdom:
        return [
            RankSpec("Root", "RootRank", "rBP", "root", ("root", "r")),
            RankSpec("Domain", "DomainRank", "dBP", "domain", ("d", "domain")),
            RankSpec("Kingdom", "KingdomRank", "kBP", "kingdom", ("k", "kingdom")),
            RankSpec("SubKingdom", "SubKingdomRank", "skBP", "subkingdom", ("sk", "subkingdom")),
            RankSpec("Phylum", "PhylumRank", "pBP", "phylum", ("p", "phylum")),
            RankSpec("Class", "ClassRank", "cBP", "class", ("c", "class")),
            RankSpec("Order", "OrderRank", "oBP", "order", ("o", "order")),
            RankSpec("Family", "FamilyRank", "fBP", "family", ("f", "family")),
            RankSpec("Genus", "GenusRank", "gBP", "genus", ("g", "genus")),
            RankSpec("Species", "SpeciesRank", "sBP", "species", ("s", "species")),
        ]

    if marker in its_fungi:
        return [
            RankSpec("Root", "RootRank", "rBP", "root", ("root", "r")),
            RankSpec("Kingdom", "KingdomRank", "kBP", "kingdom", ("k", "kingdom")),
            RankSpec("Phylum", "PhylumRank", "pBP", "phylum", ("p", "phylum")),
            RankSpec("Class", "ClassRank", "cBP", "class", ("c", "class")),
            RankSpec("Order", "OrderRank", "oBP", "order", ("o", "order")),
            RankSpec("Family", "FamilyRank", "fBP", "family", ("f", "family")),
            RankSpec("Genus", "GenusRank", "gBP", "genus", ("g", "genus")),
            RankSpec("Species", "SpeciesRank", "sBP", "species", ("s", "species")),
        ]

    raise SystemExit(
        f"Unsupported marker {marker!r} for SINTAX conversion. "
        "Add it to python_scripts/sintax_to_rdp_out.py marker_rank_specs()."
    )


_TOKEN_RE = re.compile(
    r"""
    ^
    (?P<code>[A-Za-z_]+)
    (?:
        :
        |
        __
    )
    (?P<name>.*?)
    (?:
        \((?P<conf>[0-9]*\.?[0-9]+)\)
    )?
    $
    """,
    re.VERBOSE,
)


def _normalize_code(code: str) -> str:
    code = code.strip().lower()
    if code in {"domain", "d"}:
        return "d"
    if code in {"superkingdom", "sk"}:
        return "sk"
    if code in {"subkingdom", "sub_k", "subk"}:
        return "subkingdom"
    if code in {"kingdom", "k"}:
        return "k"
    if code in {"phylum", "p"}:
        return "p"
    if code in {"class", "c"}:
        return "c"
    if code in {"order", "o"}:
        return "o"
    if code in {"family", "f"}:
        return "f"
    if code in {"genus", "g"}:
        return "g"
    if code in {"species", "s"}:
        return "s"
    if code in {"root", "r"}:
        return "root"
    return code


def parse_sintax_tax(tax: str) -> Dict[str, Tuple[str, float]]:
    """
    Parse a SINTAX `tax=...` string into {rank_code: (name, confidence)}.

    Supported token forms (best-effort):
    - k:Bacteria(1.0)
    - k__Bacteria(1.0)
    - p:Proteobacteria
    """
    tax = tax.strip()
    tax = tax.removeprefix("tax=").strip()
    tax = tax.strip(";")

    ranks: Dict[str, Tuple[str, float]] = {}
    if not tax:
        return ranks

    # Split on commas and semicolons; tolerate whitespace.
    tokens = [t.strip() for t in re.split(r"[;,]", tax) if t.strip()]
    for token in tokens:
        m = _TOKEN_RE.match(token)
        if not m:
            continue
        code = _normalize_code(m.group("code") or "")
        name = (m.group("name") or "").strip()
        if not name:
            continue
        conf_str = m.group("conf")
        conf = float(conf_str) if conf_str is not None else 0.0
        ranks[code] = (name, conf)
    return ranks


def extract_tax_field(columns: List[str]) -> str:
    """
    Best-effort extraction of the taxonomy string from a SINTAX tabbed output line.
    """
    for col in columns[1:]:
        if "tax=" in col:
            return col[col.index("tax=") :]
    # Fallback: sometimes the taxonomy is the 2nd column
    if len(columns) > 1:
        return columns[1]
    return ""


def extract_strand(columns: List[str]) -> Optional[str]:
    """
    Try to find strand (+ or -) in a SINTAX output line.
    """
    for col in columns[1:3]:
        if col in {"+", "-"}:
            return col
    return None


def format_rdp_like_row(
    query_id: str, strand: str, marker: str, ranks: Dict[str, Tuple[str, float]]
) -> str:
    specs = marker_rank_specs(marker)

    parts: List[str] = [query_id, strand]

    for spec in specs:
        chosen: Optional[Tuple[str, float]] = None
        for code in spec.accepted_codes:
            normalized = _normalize_code(code)
            if normalized in ranks:
                chosen = ranks[normalized]
                break

        if spec.rank_label == "root":
            # Root is always present in RDP output; default it if missing.
            if chosen is None:
                chosen = ("Root", 1.0)

        if chosen is None:
            name, conf = ("", 0.0)
        else:
            name, conf = chosen

        parts.extend([name, spec.rank_label, f"{conf:.4f}"])

    return "\t".join(parts)


def iter_lines(path: Path) -> Iterable[List[str]]:
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            yield line.split("\t")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="VSEARCH SINTAX --tabbedout file")
    ap.add_argument("--marker", required=True, help="Marker (controls output columns)")
    ap.add_argument("--output", required=True, help="Output path for rdp.out.tmp")
    ap.add_argument(
        "--default-strand",
        default="+",
        choices=["+", "-"],
        help="Strand to use if not present in input",
    )
    args = ap.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as out:
        for cols in iter_lines(input_path):
            query_id = cols[0]
            strand = extract_strand(cols) or args.default_strand
            tax_field = extract_tax_field(cols)
            ranks = parse_sintax_tax(tax_field)
            out.write(format_rdp_like_row(query_id, strand, args.marker, ranks) + "\n")


if __name__ == "__main__":
    main()
