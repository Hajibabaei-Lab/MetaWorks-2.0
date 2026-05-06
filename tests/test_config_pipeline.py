import importlib
from pathlib import Path

import pytest

from lib.config.config_manager import ConfigManager, ResolvedConfig
from lib.config.loader import ConfigError
from lib.config.module_registry import (
    clustering_enabled,
    is_module_enabled,
)

REPO_ROOT = Path(__file__).resolve().parent.parent


def _make_manager() -> ConfigManager:
    return ConfigManager(repo_root=str(REPO_ROOT))


class TestIsModuleEnabledCanonical:
    def test_default_true_for_trimming(self):
        config = {"modules": {}}
        assert is_module_enabled(config, "trimming") is True

    def test_default_true_for_denoising(self):
        config = {"modules": {}}
        assert is_module_enabled(config, "denoising") is True

    def test_default_false_for_clustering(self):
        config = {"modules": {}}
        assert is_module_enabled(config, "clustering") is False

    def test_default_false_for_pseudogene_filtering(self):
        config = {"modules": {}}
        assert is_module_enabled(config, "pseudogene_filtering") is False

    def test_default_false_for_stats_when_disabled(self):
        config = {"modules": {"stats": False}}
        assert is_module_enabled(config, "stats") is False

    def test_explicit_override(self):
        config = {"modules": {"clustering": True}}
        assert is_module_enabled(config, "clustering") is True

    def test_string_true(self):
        config = {"modules": {"clustering": "true"}}
        assert is_module_enabled(config, "clustering") is True

    def test_string_yes(self):
        config = {"modules": {"clustering": "yes"}}
        assert is_module_enabled(config, "clustering") is True

    def test_modules_not_dict(self):
        config = {"modules": "invalid"}
        assert is_module_enabled(config, "trimming") is True

    def test_no_modules_key(self):
        config = {}
        assert is_module_enabled(config, "trimming") is True


class TestClusteringEnabled:
    def test_clustering_off_by_default(self):
        assert clustering_enabled({}) is False
        assert clustering_enabled({"modules": {}}) is False

    def test_clustering_on_when_enabled(self):
        assert clustering_enabled({"modules": {"clustering": True}}) is True


class TestResolvedConfigFrozen:
    def test_frozen_raises_on_assignment(self):
        resolved = ResolvedConfig(
            merged={"pipeline": {"name": "esv"}},
            profile=None,
            workflow="esv",
            module_configs={},
        )
        with pytest.raises(AttributeError):
            resolved.workflow = "otu"

    def test_frozen_raises_on_merged_mutation(self):
        resolved = ResolvedConfig(
            merged={"pipeline": {"name": "esv"}},
            profile=None,
            workflow="esv",
            module_configs={},
        )
        with pytest.raises(AttributeError):
            resolved.merged = {}

    def test_export_returns_copy(self):
        data = {"pipeline": {"name": "esv"}}
        resolved = ResolvedConfig(
            merged=data, profile=None, workflow="esv", module_configs={}
        )
        exported = resolved.export_for_workflow()
        assert exported == data
        assert exported is not data

    def test_invalid_workflow_type_raises(self):
        resolved = ResolvedConfig(
            merged={}, profile=None, workflow="esv", module_configs={}
        )
        with pytest.raises(ConfigError, match="Invalid workflow type"):
            resolved.export_for_workflow("bogus")


class TestConfigManagerMergeReturnsResolved:
    def test_merge_returns_resolved_config(self):
        manager = _make_manager()
        manager.load_defaults_config()
        manager.load_module_configs()
        resolved = manager.merge(workflow="esv")
        assert isinstance(resolved, ResolvedConfig)
        assert resolved.workflow == "esv"
        assert resolved.merged["pipeline"]["name"] == "esv"

    def test_merge_without_defaults_raises(self):
        manager = _make_manager()
        with pytest.raises(ConfigError, match="Defaults config not loaded"):
            manager.merge()

    def test_merge_without_modules_raises(self):
        manager = _make_manager()
        manager.load_defaults_config()
        with pytest.raises(ConfigError, match="Module configs not loaded"):
            manager.merge()

    def test_merge_with_profile(self):
        manager = _make_manager()
        manager.load_defaults_config()
        manager.load_module_configs()
        manager.load_profile("coi")
        resolved = manager.merge(workflow="esv")
        assert resolved.profile == "coi"
        assert resolved.merged["classification"]["marker"] == "COI"

    def test_merge_with_output_dir(self):
        manager = _make_manager()
        manager.load_defaults_config()
        manager.load_module_configs()
        resolved = manager.merge(workflow="otu", output_dir="/tmp/test_output")
        assert resolved.merged["pipeline"]["output_dir"] == "/tmp/test_output"

    def test_merge_default_output_dir_esv(self):
        manager = _make_manager()
        manager.load_defaults_config()
        manager.load_module_configs()
        resolved = manager.merge(workflow="esv")
        assert resolved.merged["pipeline"]["output_dir"] == "ESV_results"

    def test_merge_default_output_dir_otu(self):
        manager = _make_manager()
        manager.load_defaults_config()
        manager.load_module_configs()
        resolved = manager.merge(workflow="otu")
        assert resolved.merged["pipeline"]["output_dir"] == "OTU_results"


class TestConfigManagerLoadReturnsResolved:
    def test_load_returns_resolved(self):
        resolved = ConfigManager.load(profile="coi", repo_root=str(REPO_ROOT))
        assert isinstance(resolved, ResolvedConfig)
        assert resolved.profile == "coi"
        assert resolved.workflow == "esv"

    def test_load_from_dict_returns_resolved(self):
        resolved = ConfigManager.load_from_dict(
            profile="coi", workflow="esv", repo_root=str(REPO_ROOT)
        )
        assert isinstance(resolved, ResolvedConfig)


class TestResolvedConfigValidate:
    def test_validate_valid_config(self):
        resolved = ConfigManager.load(profile="coi", repo_root=str(REPO_ROOT))
        errors = resolved.validate()
        assert errors == []

    def test_validate_catches_missing_dependency(self):
        manager = _make_manager()
        manager.load_defaults_config()
        manager.load_module_configs()
        manager.load_profile("coi")
        manager.user_config = {
            "modules": {
                "classification": True,
                "denoising": False,
            }
        }
        resolved = manager.merge(workflow="esv")
        errors = resolved.validate()
        assert any("requires modules.denoising" in e for e in errors)


class TestResolvedConfigGetModuleConfig:
    def test_get_module_config_trimming(self):
        resolved = ConfigManager.load(profile="coi", repo_root=str(REPO_ROOT))
        trimming = resolved.get_module_config("trimming")
        assert isinstance(trimming, dict)
        assert "min_length" in trimming

    def test_get_module_config_unknown_raises(self):
        resolved = ConfigManager.load(profile="coi", repo_root=str(REPO_ROOT))
        with pytest.raises(ConfigError, match="Module not found"):
            resolved.get_module_config("nonexistent_module")


class TestSystemConfigRemoved:
    def test_no_system_config_import(self):
        import lib.config.schemas as schemas

        assert not hasattr(schemas, "SystemConfig")

    def test_no_load_system_config_method(self):
        manager = _make_manager()
        assert not hasattr(manager, "load_system_config")
        assert not hasattr(manager, "system_config")
        assert not hasattr(manager, "system_config_path")

    def test_no_merged_config_attribute(self):
        manager = _make_manager()
        assert not hasattr(manager, "merged_config")


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


class TestITSxModule:
    def test_itsx_disabled_by_default(self):
        config = {"modules": {}}
        assert is_module_enabled(config, "itsx_extraction") is False

    def test_itsx_enabled_when_toggle_true(self):
        config = {"modules": {"itsx_extraction": True}}
        assert is_module_enabled(config, "itsx_extraction") is True

    def test_itsx_registry_entry_exists(self):
        from lib.config.module_registry import MODULE_REGISTRY

        assert "itsx_extraction" in MODULE_REGISTRY
        entry = MODULE_REGISTRY["itsx_extraction"]
        assert entry["config_section"] == "itsx_extraction"
        assert entry["stage"] == "denoising"
        assert "denoising" in entry["depends_on"]
        assert entry["module"]["enabled_by_default"] is False

    def test_itsx_included_when_enabled(self):
        from lib.config.module_registry import _resolve_included_module_names

        config = {
            "modules": {
                "trimming": True,
                "denoising": True,
                "classification": True,
                "itsx_extraction": True,
            }
        }
        names = _resolve_included_module_names(config)
        assert "itsx_extraction" in names

    def test_itsx_not_included_when_disabled(self):
        from lib.config.module_registry import _resolve_included_module_names

        config = {
            "modules": {
                "trimming": True,
                "denoising": True,
                "classification": True,
            }
        }
        names = _resolve_included_module_names(config)
        assert "itsx_extraction" not in names


class TestITSxConfig:
    def test_defaults_yaml_has_itsx_section(self):
        resolved = ConfigManager.load(repo_root=str(REPO_ROOT))
        itsx = resolved.merged.get("itsx_extraction", {})
        assert itsx.get("its_part") == "ITS2"
        assert itsx.get("threads") == 4

    def test_defaults_yaml_has_itsx_toggle(self):
        resolved = ConfigManager.load(repo_root=str(REPO_ROOT))
        assert resolved.merged["modules"]["itsx_extraction"] is False

    def test_itsx_config_model_validates_its_part(self):
        from lib.config.schemas import ITSxConfig

        config = ITSxConfig(its_part="ITS1", threads=4)
        assert config.its_part == "ITS1"

    def test_itsx_config_model_rejects_bad_its_part(self):
        from lib.config.schemas import ITSxConfig
        from pydantic import ValidationError

        with pytest.raises(ValidationError, match="its_part"):
            ITSxConfig(its_part="bad", threads=4)

    def test_module_selection_has_itsx_toggle(self):
        from lib.config.schemas import ModuleSelection

        ms = ModuleSelection()
        assert ms.itsx_extraction is False

    def test_module_selection_itsx_can_be_true(self):
        from lib.config.schemas import ModuleSelection

        ms = ModuleSelection(itsx_extraction=True)
        assert ms.itsx_extraction is True

    def test_merged_config_has_itsx_section(self):
        resolved = ConfigManager.load(profile="coi", repo_root=str(REPO_ROOT))
        assert "itsx_extraction" in resolved.merged


class TestITSxProfileMerge:
    def test_its_profile_enables_itsx(self):
        resolved = ConfigManager.load(profile="its", repo_root=str(REPO_ROOT))
        assert resolved.merged["modules"]["itsx_extraction"] is True

    def test_its_profile_sets_marker(self):
        resolved = ConfigManager.load(profile="its", repo_root=str(REPO_ROOT))
        assert resolved.merged["classification"]["marker"] == "ITS_fungi"

    def test_its_profile_has_itsx_config(self):
        resolved = ConfigManager.load(profile="its", repo_root=str(REPO_ROOT))
        itsx = resolved.merged["itsx_extraction"]
        assert itsx["its_part"] == "ITS2"
        assert itsx["threads"] == 4

    def test_its_plants_profile_enables_itsx(self):
        resolved = ConfigManager.load(profile="its_plants", repo_root=str(REPO_ROOT))
        assert resolved.merged["modules"]["itsx_extraction"] is True

    def test_its_plants_profile_sets_marker(self):
        resolved = ConfigManager.load(profile="its_plants", repo_root=str(REPO_ROOT))
        assert resolved.merged["classification"]["marker"] == "ITS_plants"

    def test_non_its_profile_keeps_itsx_disabled(self):
        resolved = ConfigManager.load(profile="coi", repo_root=str(REPO_ROOT))
        assert resolved.merged["modules"]["itsx_extraction"] is False


class TestITSxGetSequencesInput:
    def test_itsx_active_returns_stripped_path(self):
        config = {
            "pipeline": {"output_dir": "/tmp/test"},
            "modules": {"clustering": False, "itsx_extraction": True},
            "itsx_extraction": {"its_part": "ITS2"},
        }
        assert is_module_enabled(config, "itsx_extraction")
        assert not clustering_enabled(config)

    def test_itsx_inactive_returns_denoised(self):
        config = {
            "pipeline": {"output_dir": "/tmp/test"},
            "modules": {"clustering": False, "itsx_extraction": False},
        }
        assert not is_module_enabled(config, "itsx_extraction")
        assert not clustering_enabled(config)

    def test_clustering_takes_priority_over_itsx(self):
        config = {
            "pipeline": {"output_dir": "/tmp/test"},
            "modules": {"clustering": True, "itsx_extraction": True},
            "itsx_extraction": {"its_part": "ITS2"},
        }
        assert clustering_enabled(config)
        assert is_module_enabled(config, "itsx_extraction")


class TestApiRoutesWithResolvedConfig:
    def test_render_config_returns_yaml(self, monkeypatch):
        from fastapi.testclient import TestClient

        app = load_app(
            monkeypatch,
            METAWORKS_SERVE_LEGACY_UI="0",
            METAWORKS_CORS_ALLOWED_ORIGINS="http://localhost:5173",
        )
        client = TestClient(app)

        response = client.post(
            "/api/configs/render",
            json={"profile": "coi", "workflow": "esv"},
        )
        assert response.status_code == 200
        data = response.json()
        assert "merged" in data
        assert "pipeline" in data["merged"]

    def test_default_config_sections(self, monkeypatch):
        from fastapi.testclient import TestClient

        app = load_app(
            monkeypatch,
            METAWORKS_SERVE_LEGACY_UI="0",
            METAWORKS_CORS_ALLOWED_ORIGINS="http://localhost:5173",
        )
        client = TestClient(app)

        response = client.get("/api/configs/defaults/esv/sections")
        assert response.status_code == 200
        data = response.json()
        assert "sections" in data
        assert "modules" in data["sections"]
