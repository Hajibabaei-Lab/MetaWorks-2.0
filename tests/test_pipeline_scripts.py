import subprocess
import sys

import pytest


class TestRenameFastaGzip:
    def test_has_main_function(self):
        from workflow.scripts.rename_fasta_gzip import main
        assert callable(main)

    def test_renames_fasta_headers(self, sample_fasta, capsys):
        from workflow.scripts.rename_fasta_gzip import main
        main([str(sample_fasta)])
        captured = capsys.readouterr()
        assert ">test_seq1\n" in captured.out
        assert ">test_seq2\n" in captured.out

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/rename_fasta_gzip.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0
        assert "usage" in result.stderr.lower() or "error" in result.stderr.lower()


class TestFastqGzStats:
    def test_has_main_function(self):
        from workflow.scripts.fastq_gz_stats import main
        assert callable(main)

    def test_outputs_stats(self, sample_fastq_gz, capsys):
        from workflow.scripts.fastq_gz_stats import main
        main([str(sample_fastq_gz)])
        captured = capsys.readouterr()
        parts = captured.out.strip().split("\t")
        assert len(parts) == 7
        assert parts[0] == str(sample_fastq_gz)
        assert parts[1] == "2"

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/fastq_gz_stats.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0


class TestFastaGzStats:
    def test_has_main_function(self):
        from workflow.scripts.fasta_gz_stats import main
        assert callable(main)

    def test_outputs_stats(self, sample_fasta_gz, capsys):
        from workflow.scripts.fasta_gz_stats import main
        main([str(sample_fasta_gz)])
        captured = capsys.readouterr()
        parts = captured.out.strip().split("\t")
        assert len(parts) == 7
        assert parts[1] == "2"

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/fasta_gz_stats.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0


class TestAddSeqsToTax3:
    def test_has_main_function(self):
        from workflow.scripts.add_seqs_to_tax3 import main
        assert callable(main)

    def test_combines_fasta_with_rdp(self, sample_fasta_multiline, sample_rdp_out, capsys):
        from workflow.scripts.add_seqs_to_tax3 import main
        main([str(sample_fasta_multiline), str(sample_rdp_out)])
        captured = capsys.readouterr()
        lines = [line for line in captured.out.strip().split("\n") if line]
        assert len(lines) == 2
        assert lines[0].startswith("ESV1\t")
        assert "ACGTACGT" in lines[0]

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/add_seqs_to_tax3.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0


class TestAddSeqsToTax4:
    def test_has_main_function(self):
        from workflow.scripts.add_seqs_to_tax4 import main
        assert callable(main)

    def test_filters_rdp_by_fasta(self, sample_fasta_multiline, sample_rdp_out, capsys):
        from workflow.scripts.add_seqs_to_tax4 import main
        main([str(sample_fasta_multiline), str(sample_rdp_out)])
        captured = capsys.readouterr()
        lines = [line for line in captured.out.strip().split("\n") if line]
        assert len(lines) == 2
        assert lines[0].startswith("ESV1\t")
        assert "ACGTACGT" in lines[0]

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/add_seqs_to_tax4.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0


class TestGetTaxonOnly:
    def test_has_main_function(self):
        from workflow.scripts.get_taxon_only import main
        assert callable(main)

    def test_extracts_matching_sequences(self, sample_taxon_file, sample_fasta_multiline, capsys):
        from workflow.scripts.get_taxon_only import main
        main([str(sample_taxon_file), str(sample_fasta_multiline)])
        captured = capsys.readouterr()
        lines = captured.out.strip().split("\n")
        assert any(">ESV1" in line for line in lines)
        assert any(">ESV3" in line for line in lines)
        assert ">ESV2" not in captured.out

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/get_taxon_only.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0


class TestParseOrfs3:
    def test_has_main_function(self):
        from workflow.scripts.parse_orfs3 import main
        assert callable(main)

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/parse_orfs3.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0


class TestRc:
    def test_has_main_function(self):
        pytest.importorskip("Bio")
        from workflow.scripts.rc import main
        assert callable(main)

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/rc.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0


class TestFilterRdpTaxonomy:
    def test_has_main_function(self):
        pytest.importorskip("pandas")
        from workflow.scripts.filter_rdp_taxonomy import main
        assert callable(main)

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/filter_rdp_taxonomy.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0


class TestFilterEsvTable:
    def test_has_main_function(self):
        pytest.importorskip("pandas")
        from workflow.scripts.filter_ESV_table import main
        assert callable(main)

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/filter_ESV_table.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0


class TestFilterRdp:
    def test_has_main_function(self):
        from workflow.scripts.filter_rdp import main
        assert callable(main)

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/filter_rdp.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0


class TestRdpTsvToCsv:
    def test_has_main_function(self):
        pytest.importorskip("pandas")
        from workflow.scripts.rdp_tsv_to_csv import main
        assert callable(main)

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/rdp_tsv_to_csv.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0


class TestAddAbundanceToRdpOut:
    def test_has_main_function(self):
        from workflow.scripts.add_abundance_to_rdp_out import main
        assert callable(main)

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/add_abundance_to_rdp_out.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0


class TestMergeEsvTables:
    def test_has_main_function(self):
        pytest.importorskip("pandas")
        from workflow.scripts.merge_esv_tables import main
        assert callable(main)

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/merge_esv_tables.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0


class TestParseOrfs4:
    def test_has_main_function(self):
        from workflow.scripts.parse_orfs4 import main
        assert callable(main)

    def test_accepts_four_positional_args(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/parse_orfs4.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0
        assert "usage" in result.stderr.lower() or "nt_file" in result.stderr.lower()

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/parse_orfs4.py", "a", "b"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0


class TestFormatAdaptersAnchored:
    def test_has_main_function(self):
        pytest.importorskip("Bio")
        pytest.importorskip("pandas")
        from workflow.scripts.formatAdapters_anchored import main
        assert callable(main)

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/formatAdapters_anchored.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0


class TestFormatAdaptersAnchoredFilename:
    def test_has_main_function(self):
        pytest.importorskip("Bio")
        from workflow.scripts.formatAdapters_anchored_filename import main
        assert callable(main)

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/formatAdapters_anchored_filename.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0


class TestGrabSeqsFromResults:
    def test_grab_seqs_has_main_function(self):
        from workflow.scripts.grab_seqs_from_results import main
        assert callable(main)

    def test_grab_seqs_basic(self, tmp_path, capsys):
        csv_path = tmp_path / "results.csv"
        csv_path.write_text(
            "GlobalESV,SampleName,ESVsize,ESVseq,Strand,kingdom\n"
            "Zotu_1,sample1,100,ACGTACGT,+,Bacteria\n"
            "Zotu_2,sample2,50,TGCATGCA,+,Fungi\n"
        )
        from workflow.scripts.grab_seqs_from_results import main
        main([str(csv_path)])
        captured = capsys.readouterr()
        assert ">Zotu_1\n" in captured.out
        assert "ACGTACGT\n" in captured.out
        assert ">Zotu_2\n" in captured.out
        assert "TGCATGCA\n" in captured.out

    def test_grab_seqs_skips_header(self, tmp_path, capsys):
        csv_path = tmp_path / "results.csv"
        csv_path.write_text(
            "GlobalESV,SampleName,ESVsize,ESVseq,Strand,kingdom\n"
            "Zotu_1,sample1,100,ACGTACGT,+,Bacteria\n"
        )
        from workflow.scripts.grab_seqs_from_results import main
        main([str(csv_path)])
        captured = capsys.readouterr()
        assert ">GlobalESV" not in captured.out
        assert ">Zotu_1\n" in captured.out

    def test_grab_seqs_short_row_warning(self, tmp_path, capsys):
        csv_path = tmp_path / "results.csv"
        csv_path.write_text(
            "GlobalESV,SampleName,ESVsize,ESVseq\n"
            "Zotu_1,sample1,100\n"
        )
        from workflow.scripts.grab_seqs_from_results import main
        main([str(csv_path)])
        captured = capsys.readouterr()
        assert captured.out == ""
        assert "Warning" in captured.err
        assert "skipping row" in captured.err

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/grab_seqs_from_results.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0


class TestMapGlobalEsvsToResults:
    def test_map_global_esvs_has_main_function(self):
        from workflow.scripts.map_global_esvs_to_results import main
        assert callable(main)

    def test_map_global_esvs_basic(self, tmp_path, capsys):
        uc_path = tmp_path / "global.uc"
        uc_path.write_text(
            "H\t0\t100\t100.0\t100.0\t0\t0\t0\tZotu_1\tsha1abc123\n"
            "H\t0\t100\t100.0\t100.0\t0\t0\t0\tZotu_2\tsha1def456\n"
            "N\t0\t100\t*\t*\t0\t0\t0\tZotu_3\t*\n"
        )
        csv_path = tmp_path / "results.csv"
        csv_path.write_text(
            "GlobalESV,SampleName,ESVsize,ESVseq,Strand,kingdom,phylum,class,order,family,genus,species\n"
            "Zotu_1,sample1,100,ACGT,+,Bacteria,Arthropoda,Insecta,Diptera,Mosquito,Aedes,albopictus\n"
            "Zotu_2,sample1,50,TGCA,+,Bacteria,Nematoda,Chromadorea,Rhabditida,Caenorhabditis,elegans,\n"
        )
        from workflow.scripts.map_global_esvs_to_results import main
        main([str(uc_path), str(csv_path)])
        captured = capsys.readouterr()
        lines = captured.out.strip().split("\n")
        assert len(lines) == 3
        header = lines[0]
        assert "TrialESV" in header
        assert "GlobalESV" in header
        row1 = lines[1]
        assert "sha1abc123" in row1
        row2 = lines[2]
        assert "sha1def456" in row2

    def test_map_global_esvs_no_match(self, tmp_path, capsys):
        uc_path = tmp_path / "global.uc"
        uc_path.write_text(
            "N\t0\t100\t*\t*\t0\t0\t0\tZotu_3\t*\n"
        )
        csv_path = tmp_path / "results.csv"
        csv_path.write_text(
            "GlobalESV,SampleName,ESVsize,ESVseq\n"
            "Zotu_3,sample2,30,GGCC\n"
        )
        from workflow.scripts.map_global_esvs_to_results import main
        main([str(uc_path), str(csv_path)])
        captured = capsys.readouterr()
        lines = captured.out.strip().split("\n")
        assert len(lines) == 2
        assert "NoGlobalMatch" in lines[1]

    def test_map_global_esvs_empty_uc(self, tmp_path, capsys):
        uc_path = tmp_path / "global.uc"
        uc_path.write_text("")
        csv_path = tmp_path / "results.csv"
        csv_path.write_text(
            "GlobalESV,SampleName,ESVsize,ESVseq\n"
            "Zotu_1,sample1,100,ACGT\n"
        )
        from workflow.scripts.map_global_esvs_to_results import main
        main([str(uc_path), str(csv_path)])
        captured = capsys.readouterr()
        lines = captured.out.strip().split("\n")
        assert len(lines) == 2
        assert "NoGlobalMatch" in lines[1]

    def test_map_global_esvs_empty_results(self, tmp_path, capsys):
        uc_path = tmp_path / "global.uc"
        uc_path.write_text(
            "H\t0\t100\t100.0\t100.0\t0\t0\t0\tZotu_1\tsha1abc123\n"
        )
        csv_path = tmp_path / "results.csv"
        csv_path.write_text(
            "GlobalESV,SampleName,ESVsize,ESVseq\n"
        )
        from workflow.scripts.map_global_esvs_to_results import main
        main([str(uc_path), str(csv_path)])
        captured = capsys.readouterr()
        lines = [line for line in captured.out.strip().split("\n") if line]
        assert len(lines) == 1
        assert "TrialESV" in lines[0]
        assert "GlobalESV" in lines[0]

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/map_global_esvs_to_results.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0


class TestMapGlobalOtusToResults:
    def test_map_global_otus_has_main_function(self):
        from workflow.scripts.map_global_otus_to_results import main
        assert callable(main)

    def test_map_global_otus_basic(self, tmp_path, capsys):
        uc_content = "H\t0\t100\t100.0\t100.0\t0\t0\t0\tZotu_1\tOtu_001\n"
        csv_content = "GlobalOTU,SampleName,ESVsize,ESVseq\nZotu_1,sample1,100,ACGT\n"
        uc_path = tmp_path / "global.uc"
        csv_path = tmp_path / "results.csv"
        uc_path.write_text(uc_content)
        csv_path.write_text(csv_content)
        from workflow.scripts.map_global_otus_to_results import main
        main([str(uc_path), str(csv_path)])
        captured = capsys.readouterr()
        lines = captured.out.strip().split("\n")
        assert "TrialESV" in lines[0]
        assert "GlobalOTU" in lines[0]
        assert "Otu_001" in lines[1]

    def test_map_global_otus_no_match(self, tmp_path, capsys):
        uc_content = "N\t0\t100\t*\t*\t0\t0\t0\tZotu_1\t*\n"
        csv_content = "GlobalOTU,SampleName,ESVsize,ESVseq\nZotu_1,sample1,100,ACGT\n"
        uc_path = tmp_path / "global.uc"
        csv_path = tmp_path / "results.csv"
        uc_path.write_text(uc_content)
        csv_path.write_text(csv_content)
        from workflow.scripts.map_global_otus_to_results import main
        main([str(uc_path), str(csv_path)])
        captured = capsys.readouterr()
        lines = captured.out.strip().split("\n")
        assert "NoGlobalMatch" in lines[1]

    def test_map_global_otus_empty_uc(self, tmp_path, capsys):
        csv_content = "GlobalOTU,SampleName,ESVsize,ESVseq\nZotu_1,sample1,100,ACGT\n"
        uc_path = tmp_path / "global.uc"
        csv_path = tmp_path / "results.csv"
        uc_path.write_text("")
        csv_path.write_text(csv_content)
        from workflow.scripts.map_global_otus_to_results import main
        main([str(uc_path), str(csv_path)])
        captured = capsys.readouterr()
        lines = captured.out.strip().split("\n")
        assert "NoGlobalMatch" in lines[1]

    def test_missing_args_exits(self):
        result = subprocess.run(
            [sys.executable, "workflow/scripts/map_global_otus_to_results.py"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0
