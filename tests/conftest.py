import gzip

import pytest


@pytest.fixture
def sample_fasta_gz(tmp_path):
    p = tmp_path / "test.fasta.gz"
    with gzip.open(p, "wt") as f:
        f.write(">seq1\nACGTACGT\n>seq2\nTGCATGCA\n")
    return p


@pytest.fixture
def sample_fastq_gz(tmp_path):
    p = tmp_path / "test.fastq.gz"
    with gzip.open(p, "wt") as f:
        f.write(
            "@read1\nACGTACGT\n+\nIIIIIIII\n"
            "@read2\nTGCATGCA\n+\nIIIIIIII\n"
        )
    return p


@pytest.fixture
def sample_esv_table(tmp_path):
    p = tmp_path / "esv.table"
    p.write_text(
        "#OTU ID\tSampleA\tSampleB\n"
        "ESV1\t10\t5\n"
        "ESV2\t3\t7\n"
        "ESV3\t0\t12\n"
    )
    return p


@pytest.fixture
def sample_rdp_out(tmp_path):
    p = tmp_path / "rdp.out.tmp"
    p.write_text(
        "ESV1\tBacteria\t1.0\tProteobacteria\t0.9\n"
        "ESV2\tBacteria\t0.8\tFirmicutes\t0.7\n"
    )
    return p


@pytest.fixture
def sample_fasta_multiline(tmp_path):
    p = tmp_path / "multiline.fasta"
    p.write_text(">ESV1\nACGT\nACGT\n>ESV2\nTGCATGCA\n>ESV3\nAAA\n")
    return p


@pytest.fixture
def sample_orfs_fasta(tmp_path):
    p = tmp_path / "orfs.fasta"
    p.write_text(
        ">ESV1\nACGTACGT\n>ESV2\nTGCATGCA\n"
    )
    return p


@pytest.fixture
def sample_taxon_file(tmp_path):
    p = tmp_path / "taxon.zotus"
    p.write_text("ESV1\nESV3\n")
    return p
