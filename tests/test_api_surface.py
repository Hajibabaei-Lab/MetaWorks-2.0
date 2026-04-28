import importlib

from fastapi.testclient import TestClient


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


def test_api_prefix_and_legacy_aliases(monkeypatch):
    app = load_app(
        monkeypatch,
        METAWORKS_SERVE_LEGACY_UI="0",
        METAWORKS_CORS_ALLOWED_ORIGINS="http://localhost:5173",
    )
    client = TestClient(app)

    assert client.get("/api/health").status_code == 200
    assert client.get("/api/health").text == "ok"
    assert client.get("/health").status_code == 200
    assert client.get("/api/docs").status_code == 200

    schema = client.get("/api/openapi.json").json()
    assert "/api/health" in schema["paths"]
    assert "/api/runs" in schema["paths"]
    assert "/health" not in schema["paths"]
    assert "/runs" not in schema["paths"]


def test_backend_root_404s_when_legacy_ui_is_disabled(monkeypatch):
    app = load_app(monkeypatch, METAWORKS_SERVE_LEGACY_UI="0")
    client = TestClient(app)

    assert client.get("/").status_code == 404


def test_legacy_ui_mount_still_available(monkeypatch):
    app = load_app(monkeypatch, METAWORKS_SERVE_LEGACY_UI="1")
    client = TestClient(app)

    response = client.get("/")
    assert response.status_code == 200
    assert "MetaWorks" in response.text


def test_settings_paths_exposes_frontend_defaults(monkeypatch):
    app = load_app(monkeypatch, METAWORKS_SERVE_LEGACY_UI="0")
    client = TestClient(app)

    payload = client.get("/api/settings/paths").json()

    assert "default_runtime" in payload
    assert "container_image" in payload
    assert payload["container_image"]
