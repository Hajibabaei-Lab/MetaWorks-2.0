import shlex
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SMK_PATH = REPO_ROOT / "workflow" / "rules" / "pseudogene" / "pseudogene.smk"

RDP_LINES = [
    "Zotu1\t\tcellularOrganisms\tcellularOrganisms\t1.0\tEukaryota\tsuperkingdom\t1.0"
    "\tMetazoa\tkingdom\t1.0\tArthropoda\tphylum\t1.0\tInsecta\tclass\t1.0",
    "Zotu2\t\tcellularOrganisms\tcellularOrganisms\t1.0\tEukaryota\tsuperkingdom\t1.0"
    "\tMetazoa\tkingdom\t1.0\tChordata\tphylum\t1.0\tMammalia\tclass\t1.0",
    "Zotu3\t\tcellularOrganisms\tcellularOrganisms\t1.0\tEukaryota\tsuperkingdom\t1.0"
    "\tMetazoa\tkingdom\t1.0\tArthropoda\tphylum\t1.0\tCrustacea\tclass\t1.0",
    "Zotu4\t\tcellularOrganisms\tcellularOrganisms\t1.0\tEukaryota\tsuperkingdom\t1.0"
    "\tMetazoa\tkingdom\t1.0\tNematoda\tphylum\t1.0\tChromadorea\tclass\t1.0",
    "Zotu5\t\tcellularOrganisms\tcellularOrganisms\t1.0\tEukaryota\tsuperkingdom\t1.0"
    "\tMetazoa\tkingdom\t1.0\tunclassified\tphylum\t1.0"
    "\tenvironmental sample related to Arthropoda-like sequences\tclass\t1.0",
]
RDP_LINES_TAB_ONLY = RDP_LINES[:4]


@pytest.fixture
def rdp_out(tmp_path):
    p = tmp_path / "rdp.out.tmp"
    p.write_text("\n".join(RDP_LINES) + "\n")
    return p


@pytest.fixture
def rdp_out_tab_only(tmp_path):
    p = tmp_path / "rdp.out.tmp"
    p.write_text("\n".join(RDP_LINES_TAB_ONLY) + "\n")
    return p


def normalize_param(param_value):
    return " ".join(shlex.split(param_value)) or "''"


def build_fixed_command(param_value, input_path, output_path):
    normalized = normalize_param(param_value)
    return (
        "set -euo pipefail; set +o pipefail; "
        f"grep {normalized} {input_path} | awk '{{print $1}}' > \"{output_path}\" || true"
    )


def build_buggy_command(param_value, input_path, output_path):
    return (
        "set -euo pipefail; set +o pipefail; "
        f'grep "{param_value}" {input_path} | awk \'{{print $1}}\' > "{output_path}" || true'
    )


def run_command(cmd):
    return subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)


def read_ids(path):
    return [line for line in path.read_text().splitlines() if line]


class TestFixedRuleCommand:
    @pytest.mark.parametrize(
        ("param_value", "expected_ids"),
        [
            ("-e Arthropoda", ["Zotu1", "Zotu3", "Zotu5"]),
            ("Arthropoda", ["Zotu1", "Zotu3", "Zotu5"]),
            ("-v Chordata", ["Zotu1", "Zotu3", "Zotu4", "Zotu5"]),
            ("-e Chordata", ["Zotu2"]),
            ("-e Bacteria", []),
            ("  -e   Arthropoda  ", ["Zotu1", "Zotu3", "Zotu5"]),
        ],
    )
    def test_yields_expected_id_subset(self, rdp_out, tmp_path, param_value, expected_ids):
        output = tmp_path / "taxon.zotus"
        result = run_command(build_fixed_command(param_value, rdp_out, output))
        assert result.returncode == 0, result.stderr
        assert output.exists()
        assert read_ids(output) == expected_ids

    def test_no_match_yields_empty_file_not_crash(self, rdp_out, tmp_path):
        output = tmp_path / "taxon.zotus"
        result = run_command(build_fixed_command("-e Bacteria", rdp_out, output))
        assert result.returncode == 0, result.stderr
        assert output.exists()
        assert output.read_text() == ""


class TestQuotedFormRegression:
    def test_old_quoted_single_arg_matches_nothing(self, rdp_out_tab_only):
        old = run_command(f'grep "-e Arthropoda" {rdp_out_tab_only}')
        assert old.returncode == 1
        assert old.stdout == ""
        new = run_command(f"grep -e Arthropoda {rdp_out_tab_only}")
        assert new.returncode == 0
        assert len(new.stdout.strip().splitlines()) > 0

    def test_old_quoted_effective_pattern_requires_leading_space(self, rdp_out):
        old = run_command(f'grep "-e Arthropoda" {rdp_out}')
        matched = [line.split("\t")[0] for line in old.stdout.splitlines() if line]
        assert matched == ["Zotu5"]
        assert "Zotu1" not in matched
        assert "Zotu3" not in matched

    def test_old_full_rule_command_loses_all_real_arthropoda(self, rdp_out_tab_only, tmp_path):
        output = tmp_path / "taxon.zotus"
        result = run_command(build_buggy_command("-e Arthropoda", rdp_out_tab_only, output))
        assert result.returncode == 0, result.stderr
        assert read_ids(output) == []
        result = run_command(build_fixed_command("-e Arthropoda", rdp_out_tab_only, output))
        assert result.returncode == 0, result.stderr
        assert read_ids(output) == ["Zotu1", "Zotu3"]


class TestParamNormalization:
    @pytest.mark.parametrize(
        ("raw", "normalized"),
        [
            ("-e Arthropoda", "-e Arthropoda"),
            ("  -e   Arthropoda  ", "-e Arthropoda"),
            ("Arthropoda", "Arthropoda"),
            ("  -v   Chordata  ", "-v Chordata"),
        ],
    )
    def test_join_shlex_split_collapses_whitespace(self, raw, normalized):
        assert " ".join(shlex.split(raw)) == normalized

    def test_empty_param_falls_back_to_match_all_token(self):
        assert normalize_param("") == "''"
        assert normalize_param("   ") == "''"

    def test_empty_taxon_param_matches_all_ids(self, rdp_out, tmp_path):
        output = tmp_path / "taxon.zotus"
        result = run_command(build_fixed_command("", rdp_out, output))
        assert result.returncode == 0, result.stderr
        assert read_ids(output) == [line.split("\t")[0] for line in RDP_LINES]


class TestPseudogeneSmkSourceGuard:
    def _source(self):
        return SMK_PATH.read_text()

    def _rule_block(self, name):
        source = self._source()
        start = source.index(f"rule {name}:")
        next_rule = source.find("rule ", start + 1)
        return source[start:] if next_rule == -1 else source[start:next_rule]

    def test_taxon_params_normalized_with_shlex_split(self):
        source = self._source()
        assert 'shlex.split(PSEUDOGENE_CONFIG.get("taxon1"' in source
        assert 'shlex.split(PSEUDOGENE_CONFIG.get("taxon2"' in source

    def test_taxon_params_not_quoted_in_grep(self):
        source = self._source()
        assert 'grep "{params.taxon1}"' not in source
        assert 'grep "{params.taxon2}"' not in source

    def test_taxon_params_expanded_unquoted(self):
        source = self._source()
        assert "grep {params.taxon1} {input}" in source
        assert "grep {params.taxon1} {input} | grep {params.taxon2}" in source

    def test_hmmscan_guards_against_empty_orf_input(self):
        hmmscan_block = self._rule_block("hmmscan")
        assert '[ -s \\"{input.orf}\\" ]' in hmmscan_block


class TestFilterRdpEmptyHmm:
    SCRIPT = REPO_ROOT / "workflow" / "scripts" / "filter_rdp.py"

    def test_empty_hmm_and_orfs_exits_cleanly(self, tmp_path):
        hmm = tmp_path / "hmm.txt"
        orfs = tmp_path / "orfs.fasta.nt.filtered.hmm"
        rdp = tmp_path / "rdp.out.tmp"
        hmm.write_text("")
        orfs.write_text("")
        rdp.write_text("\n".join(RDP_LINES) + "\n")
        result = subprocess.run(
            [sys.executable, str(self.SCRIPT), str(hmm), str(orfs), str(rdp)],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr
        assert result.stdout == ""

    def test_populated_hmm_still_filters(self, tmp_path):
        hmm = tmp_path / "hmm.txt"
        orfs = tmp_path / "orfs.fasta.nt.filtered.hmm"
        rdp = tmp_path / "rdp.out.tmp"
        hmm.write_text(
            "# hmmscan tblout header\n"
            "bold\t-\tZotu1\t-\t1e-20\t85.0\n"
            "bold\t-\tZotu3\t-\t1e-30\t95.0\n"
        )
        orfs.write_text(">Zotu1\nACGTACGTAC\n>Zotu3\nACGTACGTAA\n")
        rdp.write_text("\n".join(RDP_LINES[:3]) + "\n")
        result = subprocess.run(
            [sys.executable, str(self.SCRIPT), str(hmm), str(orfs), str(rdp)],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr
        assert "Zotu1" in result.stdout
        assert "Zotu2" not in result.stdout
