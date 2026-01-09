# rules/utils.smk

import os
import glob
import pandas as pd
import re

# Extracts sample name from FASTQ file name
def extract_sample_name_from_file(filename):
    base = os.path.basename(filename)
    return re.sub(r"_L001_R[12]_001\.fastq\.gz$", "", base)

# Determines the path to the FASTQ file for a given sample and read
def get_fastq_path(sample, read):
    if config["input"]["sample_source"] == "csv":
        return SAMPLE_PATHS[sample].format(read=read)
    elif config["input"]["sample_source"] == "folder":
        return os.path.join(config["input"]["fastq_dir"], f"{sample}_L001_R{read}_001.fastq.gz")
    else:
        raise ValueError("Invalid sample_source config value.")

# Load sample information based on input method
if config["input"]["sample_source"] == "csv":
    samples_df = pd.read_csv(config["input"]["samples_csv"], sep=",")
    SAMPLES_UNIQUE = samples_df["sample"].tolist()
    SAMPLE_PATHS = dict(zip(samples_df["sample"], samples_df["path"]))

elif config["input"]["sample_source"] == "folder":
    FASTQ_FOLDER = config["input"]["fastq_dir"]
    samples = glob.glob(os.path.join(FASTQ_FOLDER, "*_R1_001.fastq.gz"))
    SAMPLES_UNIQUE = [extract_sample_name_from_file(f) for f in samples]
else:
    raise ValueError("Unknown sample_source: must be 'csv' or 'folder'")

# Output directory for convenience
OUTPUT_DIR = config["pipeline"]["output_dir"]

# RDP classifier options function
def rdp_options(config):
    """Generate RDP classifier options based on new config structure."""
    if config["classification"]["use_custom_classifier"]:
        return f"-t {config['classification']['classifier_path']}"
    else:
        # Using built-in classifier
        if config["classification"]["builtin_classifier"] == "fungallsu":
            return "-g fungallsu"
        elif config["classification"]["builtin_classifier"] == "fungalits_unite":
            return "-g fungalits_unite"
        elif config["classification"]["builtin_classifier"] == "fungalits_warcup":
            return "-g fungalits_warcup"
        else:
            print("ERROR: Unknown builtin_classifier specified")
            return ""


# Marker categories for downstream logic
tax3_abund6 = ['COI', 'rbcL_eukaryota', 'rbcL_landPlant', '12S_fish', '12S_vertebrate', '16S_vertebrate']
tax3_abund7 = ['16S', '28S_fungi']
tax3_abund9 = ['18S_eukaryota']
tax3_abund10 = ['18S_diatom', 'rbcL_diatom']
orf3_tax4_abund12 = ['COI', 'rbcL_landPlant', 'rbcL_eukaryota']
orf3_tax4_abund11 = ['rbcL_diatom']

HEADERS = {
    "condition6": "GlobalESV,SampleName,ESVsize,ESVseq,Strand,Root,RootRank,rBP,SuperKingdom,SuperKingdomRank,skBP,Kingdom,KingdomRank,kBP,Phylum,PhylumRank,pBP,Class,ClassRank,cBP,Order,OrderRank,oBP,Family,FamilyRank,fBP,Genus,GenusRank,gBP,Species,SpeciesRank,sBP",
    "condition7": "GlobalESV,SampleName,ESVsize,ESVseq,Strand,Domain,DomainRank,dBP,Phylum,PhylumRank,pBP,Class,ClassRank,cBP,Order,OrderRank,oBP,Family,FamilyRank,fBP,Genus,GenusRank,gBP",
    "condition8": "GlobalESV,SampleName,ESVsize,ESVseq,Strand,Root,RootRank,rBP,Kingdom,KingdomRank,kBP,Phylum,PhylumRank,pBP,Class,ClassRank,cBP,Order,OrderRank,oBP,Family,FamilyRank,fBP,Genus,GenusRank,gBP,Species,SpeciesRank,sBP",
    "condition9": "GlobalESV,SampleName,ESVsize,ESVseq,Strand,Root,RootRank,rBP,Domain,DomainRank,dBP,Kingdom,KingdomRank,kBP,Phylum,PhylumRank,pBP,Class,ClassRank,cBP,Order,OrderRank,oBP,Family,FamilyRank,fBP,Genus,GenusRank,gBP",
    "condition10": "GlobalESV,SampleName,ESVsize,ESVseq,Strand,Root,RootRank,rBP,Domain,DomainRank,dBP,Kingdom,KingdomRank,kBP,SubKingdom,SubKingdomRank,skBP,Phylum,PhylumRank,pBP,Class,ClassRank,cBP,Order,OrderRank,oBP,Family,FamilyRank,fBP,Genus,GenusRank,gBP,Species,SpeciesRank,sBP",
    "condition11": "GlobalESV,SampleName,ESVsize,ORFseq,Strand,Root,RootRank,rBP,Domain,DomainRank,dBP,Kingdom,KingdomRank,kBP,SubKingdom,SubKingdomRank,skBP,Phylum,PhylumRank,pBP,Class,ClassRank,cBP,Order,OrderRank,oBP,Family,FamilyRank,fBP,Genus,GenusRank,gBP,Species,SpeciesRank,sBP",
    "condition12": "GlobalESV,SampleName,ESVsize,ORFseq,Strand,Root,RootRank,rBP,SuperKingdom,SuperKingdomRank,skBP,Kingdom,KingdomRank,kBP,Phylum,PhylumRank,pBP,Class,ClassRank,cBP,Order,OrderRank,oBP,Family,FamilyRank,fBP,Genus,GenusRank,gBP,Species,SpeciesRank,sBP",
    "default": "",
}

def condition_key(config):
    """Determine condition key based on marker type."""
    marker = config["classification"]["marker"]
    if marker in tax3_abund6:
        return "condition6"
    elif marker in tax3_abund7:
        return "condition7"
    elif marker in tax3_abund9:
        return "condition9"
    elif marker in tax3_abund10:
        return "condition10"
    elif marker in orf3_tax4_abund11:
        return "condition11"
    elif marker in orf3_tax4_abund12:
        return "condition12"
    else:
        print("Unknown marker for results configuration.")
        return "default"


def header_value(config):
    """Resolve the header string for results based on marker."""
    key = condition_key(config)
    cfg_headers = config.get("HEADER", {}) if isinstance(config, dict) else {}
    return cfg_headers.get(key, HEADERS.get(key, ""))
