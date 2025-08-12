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
    if config["sample_source"] == "csv":
        return SAMPLE_PATHS[sample].format(read=read)
    elif config["sample_source"] == "folder":
        return os.path.join(config["fastq_dir"], f"{sample}_L001_R{read}_001.fastq.gz")
    else:
        raise ValueError("Invalid sample_source config value.")

# Load sample information based on input method
if config["sample_source"] == "csv":
    samples_df = pd.read_csv(config["samples_csv"], sep=",")
    SAMPLES_UNIQUE = samples_df["sample"].tolist()
    SAMPLE_PATHS = dict(zip(samples_df["sample"], samples_df["path"]))

elif config["sample_source"] == "folder":
    FASTQ_FOLDER = config["fastq_dir"]
    samples = glob.glob(os.path.join(FASTQ_FOLDER, "*_R1_001.fastq.gz"))
    SAMPLES_UNIQUE = [extract_sample_name_from_file(f) for f in samples]
else:
    raise ValueError("Unknown sample_source: must be 'csv' or 'folder'")

# ITSx output configuration
if config["ITSpart"] == "ITS1":
    itsx_out = config["dir"] + "/ITSx_out.ITS1.fasta"
    itsx_out2 = config["dir"] + "/ITSx_out.ITS1.fasta.2"
elif config["ITSpart"] == "ITS2":
    itsx_out = config["dir"] + "/ITSx_out.ITS2.fasta"
    itsx_out2 = config["dir"] + "/ITSx_out.ITS2.fasta.2"

def rdp_options(config):
    if config["RDP"]["custom"] == "yes":
        return f"-t {config['RDP']['t']}"
    elif config["RDP"]["custom"] == "no":
        if config["marker"] == "16S":
            return f"-c {config['RDP']['c']} -f {config['RDP']['f']}"
        elif config["marker"] == "28S_fungi":
            return f"-g {config['RDP']['g']}"
        else:
            print("ERROR: Indicate whether you are analyzing a prokaryote 16S or fungal ITS or 28S (LSU) marker in the generic_config.yaml file")
            return ""
    else:
        print("ERROR: Indicate whether you are working with a custom reference sequence database in the generic_config.yaml file")
        return ""


# Marker categories for downstream logic
tax3_abund6 = ['COI', 'rbcL_eukaryota', 'rbcL_landPlant', '12S_fish', '12S_vertebrate', '16S_vertebrate']
tax3_abund7 = ['16S', '28S_fungi']
tax3_abund9 = ['18S_eukaryota']
tax3_abund10 = ['18S_diatom', 'rbcL_diatom']
orf3_tax4_abund12 = ['COI', 'rbcL_landPlant', 'rbcL_eukaryota']
orf3_tax4_abund11 = ['rbcL_diatom']

def condition_key(config):
    if config["marker"] in tax3_abund6:
        return "condition6"
    elif config["marker"] in tax3_abund7:
        return "condition7"
    elif config["marker"] in tax3_abund9:
        return "condition9"
    elif config["marker"] in tax3_abund10:
        return "condition10"
    elif config["marker"] in orf3_tax4_abund11:
        return "condition11"
    elif config["marker"] in orf3_tax4_abund12:
        return "condition12"
    else:
        print("Unknown marker for results configuration.")
        return "default"
