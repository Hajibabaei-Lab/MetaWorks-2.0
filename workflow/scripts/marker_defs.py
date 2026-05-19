"""Canonical marker definitions for the MetaWorks 2.0 pipeline.

SINGLE SOURCE OF TRUTH for all marker-to-column-layout mappings.
Import from this module instead of defining marker lists locally.

Usage::

    from marker_defs import get_header, get_condition, get_all_markers

    header = get_header("COI")
    header_pseudo = get_header("COI", pseudogene=True)
"""

from __future__ import annotations

from typing import Dict, List, Tuple


# ---------------------------------------------------------------------------
# 1. Condition headers – full column strings including the 5-field prefix
# ---------------------------------------------------------------------------

HEADERS: Dict[str, str] = {
    "condition6": (
        "GlobalESV,SampleName,ESVsize,ESVseq,Strand,"
        "Root,RootRank,rBP,"
        "SuperKingdom,SuperKingdomRank,skBP,"
        "Kingdom,KingdomRank,kBP,"
        "Phylum,PhylumRank,pBP,"
        "Class,ClassRank,cBP,"
        "Order,OrderRank,oBP,"
        "Family,FamilyRank,fBP,"
        "Genus,GenusRank,gBP,"
        "Species,SpeciesRank,sBP"
    ),
    "condition7": (
        "GlobalESV,SampleName,ESVsize,ESVseq,Strand,"
        "Domain,DomainRank,dBP,"
        "Phylum,PhylumRank,pBP,"
        "Class,ClassRank,cBP,"
        "Order,OrderRank,oBP,"
        "Family,FamilyRank,fBP,"
        "Genus,GenusRank,gBP"
    ),
    "condition8": (
        "GlobalESV,SampleName,ESVsize,ESVseq,Strand,"
        "Root,RootRank,rBP,"
        "Kingdom,KingdomRank,kBP,"
        "Phylum,PhylumRank,pBP,"
        "Class,ClassRank,cBP,"
        "Order,OrderRank,oBP,"
        "Family,FamilyRank,fBP,"
        "Genus,GenusRank,gBP,"
        "Species,SpeciesRank,sBP"
    ),
    "condition9": (
        "GlobalESV,SampleName,ESVsize,ESVseq,Strand,"
        "Root,RootRank,rBP,"
        "Domain,DomainRank,dBP,"
        "Kingdom,KingdomRank,kBP,"
        "Phylum,PhylumRank,pBP,"
        "Class,ClassRank,cBP,"
        "Order,OrderRank,oBP,"
        "Family,FamilyRank,fBP,"
        "Genus,GenusRank,gBP"
    ),
    "condition10": (
        "GlobalESV,SampleName,ESVsize,ESVseq,Strand,"
        "Root,RootRank,rBP,"
        "Domain,DomainRank,dBP,"
        "Kingdom,KingdomRank,kBP,"
        "SubKingdom,SubKingdomRank,skBP,"
        "Phylum,PhylumRank,pBP,"
        "Class,ClassRank,cBP,"
        "Order,OrderRank,oBP,"
        "Family,FamilyRank,fBP,"
        "Genus,GenusRank,gBP,"
        "Species,SpeciesRank,sBP"
    ),
    "condition11": (
        "GlobalESV,SampleName,ESVsize,ORFseq,Strand,"
        "Root,RootRank,rBP,"
        "Domain,DomainRank,dBP,"
        "Kingdom,KingdomRank,kBP,"
        "SubKingdom,SubKingdomRank,skBP,"
        "Phylum,PhylumRank,pBP,"
        "Class,ClassRank,cBP,"
        "Order,OrderRank,oBP,"
        "Family,FamilyRank,fBP,"
        "Genus,GenusRank,gBP,"
        "Species,SpeciesRank,sBP"
    ),
    "condition12": (
        "GlobalESV,SampleName,ESVsize,ORFseq,Strand,"
        "Root,RootRank,rBP,"
        "SuperKingdom,SuperKingdomRank,skBP,"
        "Kingdom,KingdomRank,kBP,"
        "Phylum,PhylumRank,pBP,"
        "Class,ClassRank,cBP,"
        "Order,OrderRank,oBP,"
        "Family,FamilyRank,fBP,"
        "Genus,GenusRank,gBP,"
        "Species,SpeciesRank,sBP"
    ),
}

# Number of prefix columns added by downstream scripts (GlobalESV..Strand)
_PREFIX_LEN = 5


# ---------------------------------------------------------------------------
# 2. Marker-to-condition mappings
# ---------------------------------------------------------------------------

MARKER_TO_CONDITION: Dict[str, str] = {
    "COI": "condition6",
    "rbcL_eukaryota": "condition6",
    "rbcL_landPlant": "condition6",
    "12S_fish": "condition6",
    "12S_vertebrate": "condition6",
    "16S_vertebrate": "condition6",
    "16S": "condition7",
    "28S_fungi": "condition7",
    "ITS_fungi": "condition8",
    "ITS_plants": "condition8",
    "18S_eukaryota": "condition9",
    "18S_diatom": "condition10",
    "rbcL_diatom": "condition10",
}

MARKER_TO_CONDITION_PSEUDOGENE: Dict[str, str] = {
    "COI": "condition12",
    "rbcL_eukaryota": "condition12",
    "rbcL_landPlant": "condition12",
    "rbcL_diatom": "condition11",
}


# ---------------------------------------------------------------------------
# 3. Backward-compatibility marker lists
# ---------------------------------------------------------------------------

TAX3_ABUND6_MARKERS: List[str] = [
    "COI",
    "rbcL_eukaryota",
    "rbcL_landPlant",
    "12S_fish",
    "12S_vertebrate",
    "16S_vertebrate",
]

TAX3_ABUND7_MARKERS: List[str] = ["16S", "28S_fungi"]

TAX3_ABUND8_MARKERS: List[str] = ["ITS_fungi", "ITS_plants"]

TAX3_ABUND9_MARKERS: List[str] = ["18S_eukaryota"]

TAX3_ABUND10_MARKERS: List[str] = ["18S_diatom", "rbcL_diatom"]

ORF3_TAX4_ABUND11_MARKERS: List[str] = ["rbcL_diatom"]

ORF3_TAX4_ABUND12_MARKERS: List[str] = ["COI", "rbcL_landPlant", "rbcL_eukaryota"]


# ---------------------------------------------------------------------------
# 4. Rank specifications for SINTAX conversion
# ---------------------------------------------------------------------------

# (rank_name, code_letter) per condition – used by sintax_to_rdp_out.py
_RANK_SPECS: Dict[str, List[Tuple[str, str]]] = {
    "condition6": [
        ("Root", "ro"),
        ("SuperKingdom", "sk"),
        ("Kingdom", "k"),
        ("Phylum", "p"),
        ("Class", "c"),
        ("Order", "o"),
        ("Family", "f"),
        ("Genus", "g"),
        ("Species", "s"),
    ],
    "condition7": [
        ("Domain", "d"),
        ("Phylum", "p"),
        ("Class", "c"),
        ("Order", "o"),
        ("Family", "f"),
        ("Genus", "g"),
    ],
    "condition8": [
        ("Root", "ro"),
        ("Kingdom", "k"),
        ("Phylum", "p"),
        ("Class", "c"),
        ("Order", "o"),
        ("Family", "f"),
        ("Genus", "g"),
        ("Species", "s"),
    ],
    "condition9": [
        ("Root", "ro"),
        ("Domain", "d"),
        ("Kingdom", "k"),
        ("Phylum", "p"),
        ("Class", "c"),
        ("Order", "o"),
        ("Family", "f"),
        ("Genus", "g"),
    ],
    "condition10": [
        ("Root", "ro"),
        ("Domain", "d"),
        ("Kingdom", "k"),
        ("SubKingdom", "sk"),
        ("Phylum", "p"),
        ("Class", "c"),
        ("Order", "o"),
        ("Family", "f"),
        ("Genus", "g"),
        ("Species", "s"),
    ],
    "condition11": [
        ("Root", "ro"),
        ("Domain", "d"),
        ("Kingdom", "k"),
        ("SubKingdom", "sk"),
        ("Phylum", "p"),
        ("Class", "c"),
        ("Order", "o"),
        ("Family", "f"),
        ("Genus", "g"),
        ("Species", "s"),
    ],
    "condition12": [
        ("Root", "ro"),
        ("SuperKingdom", "sk"),
        ("Kingdom", "k"),
        ("Phylum", "p"),
        ("Class", "c"),
        ("Order", "o"),
        ("Family", "f"),
        ("Genus", "g"),
        ("Species", "s"),
    ],
}


# ---------------------------------------------------------------------------
# 5. Public convenience functions
# ---------------------------------------------------------------------------


def get_condition(marker: str, pseudogene: bool = False) -> str:
    """Return the condition key for *marker*.

    Parameters
    ----------
    marker:
        Marker name (e.g. ``"COI"``).
    pseudogene:
        When *True*, return the pseudogene-specific condition for markers
        that have a different layout (ORFseq instead of ESVseq).
    """
    marker = marker.strip()
    if pseudogene and marker in MARKER_TO_CONDITION_PSEUDOGENE:
        return MARKER_TO_CONDITION_PSEUDOGENE[marker]
    if marker in MARKER_TO_CONDITION:
        return MARKER_TO_CONDITION[marker]
    raise ValueError(f"Unknown marker: {marker!r}")


def get_header(marker: str, pseudogene: bool = False) -> str:
    """Return the full CSV header string for *marker*."""
    return HEADERS[get_condition(marker, pseudogene)]


def get_rdp_csv_header(marker: str, pseudogene: bool = False) -> str:
    """Return the taxonomy-only header (without the 5-column prefix).

    Suitable for scripts like *rdp_tsv_to_csv.py* that prepend their own
    ``GlobalESV,SampleName,ESVsize,ESVseq/ORFseq,Strand`` columns.
    """
    full = get_header(marker, pseudogene)
    parts = full.split(",")
    return ",".join(parts[_PREFIX_LEN:])


def get_all_markers() -> List[str]:
    """Return all supported marker names (sorted)."""
    return sorted(MARKER_TO_CONDITION.keys())


def get_rank_specs(marker: str, pseudogene: bool = False) -> List[Tuple[str, str]]:
    """Return rank specifications for *marker* as ``(rank_name, code)`` tuples.

    Each tuple maps a taxonomy rank name to its short code used in SINTAX
    taxonomy strings.  For example::

        >>> get_rank_specs("COI")
        [("Root", "ro"), ("SuperKingdom", "sk"), ("Kingdom", "k"), ...]
    """
    condition = get_condition(marker, pseudogene)
    return list(_RANK_SPECS[condition])
