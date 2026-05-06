import importlib
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from lib.config.loader import load_yaml
from lib.config.schema_builder import build_config_schema

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULTS_PATH = REPO_ROOT / "config" / "defaults.yaml"
PROFILES_DIR = REPO_ROOT / "config" / "presets"


def _load_defaults():
    return load_yaml(str(DEFAULTS_PATH))


def _load_profile(name):
    return load_yaml(str(PROFILES_DIR / f"{name}.yaml"))


def load_app(monkeypatch, **env):
    for key, value in env.items():
        if value is None:
            monkeypatch.delenv(key, raising=False)
        else:
            monkeypatch.setenv(key, value)

    import api.config as config_module
    import api.job_manager as job_manager_module
    import api.main as main_module

    importlib.reload(config_module)
    importlib.reload(job_manager_module)
    main_module = importlib.reload(main_module)
    return main_module.app


class TestSchemaTopLevelStructure:
    def test_schema_has_all_sections(self):
        defaults = _load_defaults()
        profile = _load_profile("coi")
        schema = build_config_schema(defaults, profile, "coi")

        expected_sections = [
            "pipeline",
            "input",
            "modules",
            "trimming",
            "denoising",
            "classification",
            "pseudogene_filtering",
            "stats",
            "output",
        ]
        for section in expected_sections:
            assert section in schema["sections"], f"Missing section: {section}"

    def test_schema_top_level_keys(self):
        defaults = _load_defaults()
        profile = _load_profile("coi")
        schema = build_config_schema(defaults, profile, "coi")

        assert schema["profile"] == "coi"
        assert schema["marker"] == "COI"
        assert schema["workflow"] == "esv"
        assert isinstance(schema["sections"], dict)


class TestSchemaPipelineSection:
    def test_schema_pipeline_section(self):
        defaults = _load_defaults()
        profile = _load_profile("coi")
        schema = build_config_schema(defaults, profile, "coi")

        pipeline = schema["sections"]["pipeline"]
        assert pipeline["label"] == "Pipeline Settings"
        fields = pipeline["fields"]

        assert "parallel_jobs" in fields
        assert "name" not in fields, "'name' should be excluded (system-managed)"
        assert "output_dir" not in fields, "'output_dir' should be excluded (system-managed)"

        pj = fields["parallel_jobs"]
        assert pj["type"] == "integer"
        assert pj["default"] == 4
        assert pj["constraints"]["ge"] == 1
        assert pj["constraints"]["le"] == 32


class TestSchemaModulesSection:
    def test_schema_modules_section(self):
        defaults = _load_defaults()
        profile = _load_profile("coi")
        schema = build_config_schema(defaults, profile, "coi")

        modules = schema["sections"]["modules"]
        assert modules["label"] == "Module Toggles"
        fields = modules["fields"]

        boolean_toggles = [
            "trimming",
            "denoising",
            "classification",
            "pseudogene_filtering",
            "stats",
        ]
        for toggle in boolean_toggles:
            assert toggle in fields, f"Missing module toggle: {toggle}"
            assert fields[toggle]["type"] == "boolean"
            assert isinstance(fields[toggle]["default"], bool)

        assert "classification_engine" in fields
        ce = fields["classification_engine"]
        assert ce["type"] == "select"
        assert set(ce["options"]) == {"rdp", "sintax"}


class TestSchemaTrimmingSection:
    def test_schema_trimming_section(self):
        defaults = _load_defaults()
        profile = _load_profile("coi")
        schema = build_config_schema(defaults, profile, "coi")

        trimming = schema["sections"]["trimming"]
        assert trimming["label"] == "Trimming Parameters (Cutadapt)"
        assert trimming["enabled_by"] == "modules.trimming"
        assert trimming["collapsed"] is True
        fields = trimming["fields"]

        expected_fields = [
            "adapters",
            "min_length",
            "quality_cutoff",
            "error_rate",
            "min_adapter_overlap",
            "max_n_bases",
            "enable_rc",
        ]
        for field_name in expected_fields:
            assert field_name in fields, f"Missing trimming field: {field_name}"

        assert fields["adapters"]["type"] == "file_ref"
        assert fields["min_length"]["type"] == "integer"
        assert fields["min_length"]["constraints"]["ge"] == 10
        assert fields["min_length"]["constraints"]["le"] == 10000
        assert fields["enable_rc"]["type"] == "boolean"
        assert fields["enable_rc"]["default"] is True


class TestSchemaClassificationSection:
    def test_schema_classification_has_groups(self):
        defaults = _load_defaults()
        profile = _load_profile("coi")
        schema = build_config_schema(defaults, profile, "coi")

        classification = schema["sections"]["classification"]
        assert classification["label"] == "Classification Parameters"
        assert classification["enabled_by"] == "modules.classification"
        assert classification["collapsed"] is True
        fields = classification["fields"]

        for group_name in ["rdp", "sintax"]:
            assert group_name in fields, f"Missing classification group: {group_name}"
            group = fields[group_name]
            assert group["type"] == "group"
            assert "visible_when" in group
            assert "fields" in group
            assert isinstance(group["fields"], dict)

        assert (
            fields["rdp"]["visible_when"]
            == "modules.classification_engine == 'rdp'"
        )
        assert (
            fields["sintax"]["visible_when"]
            == "modules.classification_engine == 'sintax'"
        )

    def test_schema_classification_rdp_fields(self):
        defaults = _load_defaults()
        profile = _load_profile("coi")
        schema = build_config_schema(defaults, profile, "coi")

        rdp = schema["sections"]["classification"]["fields"]["rdp"]
        rdp_fields = rdp["fields"]

        assert "memory_gb" in rdp_fields
        assert rdp_fields["memory_gb"]["type"] == "integer"
        assert rdp_fields["memory_gb"]["default"] == 20

        assert "use_custom_classifier" in rdp_fields
        assert rdp_fields["use_custom_classifier"]["type"] == "boolean"

        assert "classifier_path" in rdp_fields
        assert rdp_fields["classifier_path"]["type"] == "file_ref"

        assert "builtin_classifier" in rdp_fields

    def test_schema_classification_sintax_fields(self):
        defaults = _load_defaults()
        profile = _load_profile("coi")
        schema = build_config_schema(defaults, profile, "coi")

        sintax = schema["sections"]["classification"]["fields"]["sintax"]
        sintax_fields = sintax["fields"]

        assert "db_fasta" in sintax_fields
        assert sintax_fields["db_fasta"]["type"] == "file_ref"

        assert "cutoff" in sintax_fields

        assert "threads" in sintax_fields
        assert sintax_fields["threads"]["type"] == "integer"

    def test_schema_classification_marker_field(self):
        defaults = _load_defaults()
        profile = _load_profile("coi")
        schema = build_config_schema(defaults, profile, "coi")

        classification = schema["sections"]["classification"]
        assert "marker" in classification["fields"]
        marker = classification["fields"]["marker"]
        assert marker["default"] == "COI"

    def test_schema_classification_min_confidence(self):
        defaults = _load_defaults()
        profile = _load_profile("coi")
        schema = build_config_schema(defaults, profile, "coi")

        classification = schema["sections"]["classification"]
        assert "min_confidence" in classification["fields"]
        mc = classification["fields"]["min_confidence"]
        assert mc["type"] == "float"
        assert mc["default"] == 0.8


class TestSchemaPseudogeneSection:
    def test_schema_pseudogene_filtering(self):
        defaults = _load_defaults()
        profile = _load_profile("coi")
        schema = build_config_schema(defaults, profile, "coi")

        pg = schema["sections"]["pseudogene_filtering"]
        assert pg["label"] == "Pseudogene Filtering"
        assert pg["enabled_by"] == "modules.pseudogene_filtering"
        assert pg["collapsed"] is True

        fields = pg["fields"]
        assert "method" in fields
        method = fields["method"]
        assert method["type"] == "select"
        assert set(method["options"]) == {"hmm", "orf"}

    def test_schema_pseudogene_all_fields(self):
        defaults = _load_defaults()
        profile = _load_profile("coi")
        schema = build_config_schema(defaults, profile, "coi")

        pg_fields = schema["sections"]["pseudogene_filtering"]["fields"]

        expected = [
            "method",
            "grep_type",
            "taxon1",
            "taxon2",
            "hmm_profile",
            "genetic_code",
            "orf_start_codon",
            "min_orf_length",
            "ignore_nested_orfs",
            "strand",
        ]
        for field_name in expected:
            assert field_name in pg_fields, f"Missing pseudogene field: {field_name}"

        assert pg_fields["hmm_profile"]["type"] == "file_ref"
        assert pg_fields["genetic_code"]["type"] == "select"
        assert pg_fields["ignore_nested_orfs"]["type"] == "boolean"


class TestSchemaProfileMerge:
    def test_schema_profile_merge_16s(self):
        defaults = _load_defaults()
        profile_16s = _load_profile("16s")
        schema = build_config_schema(defaults, profile_16s, "16s")

        assert schema["marker"] == "16S"

        classification = schema["sections"]["classification"]
        marker_field = classification["fields"]["marker"]
        assert marker_field["default"] == "16S"

    def test_schema_profile_merge_16s_genetic_code(self):
        defaults = _load_defaults()
        profile_16s = _load_profile("16s")
        schema = build_config_schema(defaults, profile_16s, "16s")

        pg = schema["sections"]["pseudogene_filtering"]
        assert pg["fields"]["genetic_code"]["default"] == 11

    def test_schema_profile_merge_16s_builtin_classifier(self):
        defaults = _load_defaults()
        profile_16s = _load_profile("16s")
        schema = build_config_schema(defaults, profile_16s, "16s")

        rdp = schema["sections"]["classification"]["fields"]["rdp"]
        assert rdp["fields"]["builtin_classifier"]["default"] == "rdp_train_16s"

    def test_schema_profile_merge_16s_method(self):
        defaults = _load_defaults()
        profile_16s = _load_profile("16s")
        schema = build_config_schema(defaults, profile_16s, "16s")

        pg = schema["sections"]["pseudogene_filtering"]
        assert pg["fields"]["method"]["default"] == "orf"

    def test_schema_profile_merge_coi_classifier_path(self):
        defaults = _load_defaults()
        profile_coi = _load_profile("coi")
        schema = build_config_schema(defaults, profile_coi, "coi")

        rdp = schema["sections"]["classification"]["fields"]["rdp"]
        assert rdp["fields"]["classifier_path"]["default"] == "runtime/classifiers/COI.properties"

    def test_schema_profile_merge_coi_hmm_profile(self):
        defaults = _load_defaults()
        profile_coi = _load_profile("coi")
        schema = build_config_schema(defaults, profile_coi, "coi")

        pg = schema["sections"]["pseudogene_filtering"]
        assert pg["fields"]["hmm_profile"]["default"] == "config/hmm/bold.hmm"


class TestSchemaInputSection:
    def test_schema_input_section(self):
        defaults = _load_defaults()
        profile = _load_profile("coi")
        schema = build_config_schema(defaults, profile, "coi")

        input_section = schema["sections"]["input"]
        assert input_section["label"] == "Input Settings"
        fields = input_section["fields"]

        assert "sample_source" in fields
        assert fields["sample_source"]["type"] == "select"
        assert fields["sample_source"]["default"] == "folder"
        assert fields["sample_source"]["options"] == ["folder", "csv"]

        assert "fastq_dir" in fields
        assert fields["fastq_dir"]["type"] == "file_ref"

        assert "samples_csv" in fields
        assert fields["samples_csv"]["type"] == "file_ref"


class TestSchemaNullDefaults:
    def test_schema_null_defaults(self):
        defaults = _load_defaults()
        profile_16s = _load_profile("16s")
        schema = build_config_schema(defaults, profile_16s, "16s")

        classification = schema["sections"]["classification"]

        rdp = classification["fields"]["rdp"]
        assert rdp["fields"]["classifier_path"]["default"] is None

        sintax = classification["fields"]["sintax"]
        assert sintax["fields"]["db_fasta"]["default"] is None
        assert sintax["fields"]["cutoff"]["default"] is None

    def test_schema_null_pseudogene_hmm_default(self):
        defaults = _load_defaults()
        schema = build_config_schema(defaults, {}, "minimal")

        pg = schema["sections"]["pseudogene_filtering"]
        assert pg["fields"]["hmm_profile"]["default"] is None


class TestSchemaWorkflowParam:
    def test_schema_default_workflow(self):
        defaults = _load_defaults()
        profile = _load_profile("coi")
        schema = build_config_schema(defaults, profile, "coi")
        assert schema["workflow"] == "esv"

    def test_schema_otu_workflow(self):
        defaults = _load_defaults()
        profile = _load_profile("coi")
        schema = build_config_schema(defaults, profile, "coi", workflow="otu")
        assert schema["workflow"] == "otu"


class TestSchemaSectionEnabledBy:
    @pytest.mark.parametrize(
        "section_name,expected_enabled_by",
        [
            ("trimming", "modules.trimming"),
            ("denoising", "modules.denoising"),
            ("classification", "modules.classification"),
            ("itsx_extraction", "modules.itsx_extraction"),
            ("pseudogene_filtering", "modules.pseudogene_filtering"),
            ("stats", "modules.stats"),
        ],
    )
    def test_section_enabled_by(self, section_name, expected_enabled_by):
        defaults = _load_defaults()
        profile = _load_profile("coi")
        schema = build_config_schema(defaults, profile, "coi")

        section = schema["sections"][section_name]
        assert section["enabled_by"] == expected_enabled_by
        assert section["collapsed"] is True

    def test_pipeline_no_enabled_by(self):
        defaults = _load_defaults()
        profile = _load_profile("coi")
        schema = build_config_schema(defaults, profile, "coi")

        pipeline = schema["sections"]["pipeline"]
        assert "enabled_by" not in pipeline
        assert "collapsed" not in pipeline

    def test_input_no_enabled_by(self):
        defaults = _load_defaults()
        profile = _load_profile("coi")
        schema = build_config_schema(defaults, profile, "coi")

        input_section = schema["sections"]["input"]
        assert "enabled_by" not in input_section
        assert "collapsed" not in input_section


class TestSchemaEndpoint:
    def test_api_schema_endpoint_coi(self, monkeypatch):
        app = load_app(
            monkeypatch,
            METAWORKS_SERVE_LEGACY_UI="0",
            METAWORKS_CORS_ALLOWED_ORIGINS="http://localhost:5173",
        )
        client = TestClient(app)

        response = client.get("/api/config/schema", params={"profile": "coi"})
        assert response.status_code == 200

        data = response.json()
        assert data["profile"] == "coi"
        assert data["marker"] == "COI"
        assert data["workflow"] == "esv"
        assert "sections" in data
        assert isinstance(data["sections"], dict)
        assert len(data["sections"]) > 0

        for section_name, section_data in data["sections"].items():
            assert "label" in section_data, f"Section '{section_name}' missing 'label'"
            assert "fields" in section_data, f"Section '{section_name}' missing 'fields'"

    def test_api_schema_endpoint_16s(self, monkeypatch):
        app = load_app(
            monkeypatch,
            METAWORKS_SERVE_LEGACY_UI="0",
            METAWORKS_CORS_ALLOWED_ORIGINS="http://localhost:5173",
        )
        client = TestClient(app)

        response = client.get("/api/config/schema", params={"profile": "16s"})
        assert response.status_code == 200

        data = response.json()
        assert data["profile"] == "16s"
        assert data["marker"] == "16S"

    def test_api_schema_endpoint_its(self, monkeypatch):
        app = load_app(
            monkeypatch,
            METAWORKS_SERVE_LEGACY_UI="0",
            METAWORKS_CORS_ALLOWED_ORIGINS="http://localhost:5173",
        )
        client = TestClient(app)

        response = client.get("/api/config/schema", params={"profile": "its"})
        assert response.status_code == 200

        data = response.json()
        assert data["profile"] == "its"
        assert data["marker"] == "ITS_fungi"

    def test_api_schema_endpoint_not_found(self, monkeypatch):
        app = load_app(
            monkeypatch,
            METAWORKS_SERVE_LEGACY_UI="0",
            METAWORKS_CORS_ALLOWED_ORIGINS="http://localhost:5173",
        )
        client = TestClient(app)

        response = client.get(
            "/api/config/schema", params={"profile": "nonexistent"}
        )
        assert response.status_code == 404

    def test_api_schema_default_profile(self, monkeypatch):
        app = load_app(
            monkeypatch,
            METAWORKS_SERVE_LEGACY_UI="0",
            METAWORKS_CORS_ALLOWED_ORIGINS="http://localhost:5173",
        )
        client = TestClient(app)

        response = client.get("/api/config/schema")
        assert response.status_code == 200

        data = response.json()
        assert data["profile"] == "coi"

    def test_api_schema_workflow_param(self, monkeypatch):
        app = load_app(
            monkeypatch,
            METAWORKS_SERVE_LEGACY_UI="0",
            METAWORKS_CORS_ALLOWED_ORIGINS="http://localhost:5173",
        )
        client = TestClient(app)

        response = client.get(
            "/api/config/schema",
            params={"profile": "coi", "workflow": "otu"},
        )
        assert response.status_code == 200

        data = response.json()
        assert data["workflow"] == "otu"

    def test_api_schema_in_openapi(self, monkeypatch):
        app = load_app(
            monkeypatch,
            METAWORKS_SERVE_LEGACY_UI="0",
            METAWORKS_CORS_ALLOWED_ORIGINS="http://localhost:5173",
        )
        client = TestClient(app)

        schema = client.get("/api/openapi.json").json()
        assert "/api/config/schema" in schema["paths"]

        schema_endpoint = schema["paths"]["/api/config/schema"]
        assert "get" in schema_endpoint

    def test_api_schema_response_matches_model(self, monkeypatch):
        from api.schemas import ConfigSchemaResponse

        app = load_app(
            monkeypatch,
            METAWORKS_SERVE_LEGACY_UI="0",
            METAWORKS_CORS_ALLOWED_ORIGINS="http://localhost:5173",
        )
        client = TestClient(app)

        response = client.get("/api/config/schema", params={"profile": "coi"})
        assert response.status_code == 200

        parsed = ConfigSchemaResponse(**response.json())
        assert parsed.profile == "coi"
        assert parsed.marker == "COI"
        assert len(parsed.sections) >= 9

    def test_api_schema_sections_have_expected_structure(self, monkeypatch):
        from api.schemas import SectionSchema

        app = load_app(
            monkeypatch,
            METAWORKS_SERVE_LEGACY_UI="0",
            METAWORKS_CORS_ALLOWED_ORIGINS="http://localhost:5173",
        )
        client = TestClient(app)

        response = client.get("/api/config/schema", params={"profile": "coi"})
        data = response.json()

        for section_name, section_raw in data["sections"].items():
            section = SectionSchema(**section_raw)
            assert section.label
            assert isinstance(section.fields, dict)
            if section.enabled_by:
                assert section.collapsed is not None


class TestSchemaITSxSection:
    def test_itsx_section_present_in_schema(self):
        defaults = _load_defaults()
        schema = build_config_schema(defaults, {}, "default", "esv")
        assert "itsx_extraction" in schema["sections"]

    def test_itsx_section_has_correct_fields(self):
        defaults = _load_defaults()
        schema = build_config_schema(defaults, {}, "default", "esv")
        itsx = schema["sections"]["itsx_extraction"]
        assert "its_part" in itsx["fields"]
        assert "threads" in itsx["fields"]

    def test_itsx_section_enabled_by_module_toggle(self):
        defaults = _load_defaults()
        schema = build_config_schema(defaults, {}, "default", "esv")
        itsx = schema["sections"]["itsx_extraction"]
        assert itsx["enabled_by"] == "modules.itsx_extraction"

    def test_itsx_its_part_is_select(self):
        defaults = _load_defaults()
        schema = build_config_schema(defaults, {}, "default", "esv")
        field = schema["sections"]["itsx_extraction"]["fields"]["its_part"]
        assert field["type"] == "select"
        assert "ITS1" in field["options"]
        assert "ITS2" in field["options"]

    def test_itsx_threads_is_integer(self):
        defaults = _load_defaults()
        schema = build_config_schema(defaults, {}, "default", "esv")
        field = schema["sections"]["itsx_extraction"]["fields"]["threads"]
        assert field["type"] == "integer"
        assert field["default"] == 4

    def test_itsx_modules_toggle_in_schema(self):
        defaults = _load_defaults()
        schema = build_config_schema(defaults, {}, "default", "esv")
        modules_fields = schema["sections"]["modules"]["fields"]
        assert "itsx_extraction" in modules_fields
        assert modules_fields["itsx_extraction"]["type"] == "boolean"
        assert modules_fields["itsx_extraction"]["default"] is False
