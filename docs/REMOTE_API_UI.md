# Control Node API and Web UI (single node)

This document captures the single-node control stack for MetaWorks: a FastAPI backend plus a Vue 3 SPA frontend that can launch Snakemake-based ESV/OTU runs using host conda, Docker, or experimental Apptainer configuration without changing the UI.

## Architecture

- **FastAPI service** in `api/`
  - Endpoints for run submission/status/log tail/artifact packaging/cancel.
  - Config templating (`/configs/defaults/{workflow}`, `/configs/render`).
  - Managed assets for classifiers/adapters (upload/list/delete).
  - Serves the API at `/api/*`; the Vue 3 SPA frontend is served separately by Caddy (split deployment) or during dev by the Vite dev server.
  - Tracks run metadata, PIDs/return codes, and log paths in `runtime/state/runs.json`.
- **Runner** (single node, pluggable): builds and launches Snakemake commands based on `runtime` (`conda`, `docker`, or `apptainer`). Uses local `subprocess.Popen` today; can be swapped for a scheduler submitter later without UI changes.
- **State + data layout** (host-visible so each runtime sees the same files):
  - `runtime/runs/` per-run working directories (rendered configs, logs).
  - `runtime/artifacts/` packaged tarballs.
  - `runtime/classifiers/` and `runtime/adapters/` managed uploads.
  - `runtime/cache/` default Apptainer/Singularity cache/prefix for container pulls.

## Runtime options (pluggable)

- `conda`: Runs Snakemake directly in the active environment:
  `snakemake -s <Snakefile> --configfile <rendered> --cores N --printshellcmds`
- `docker`: Runs Snakemake through the Docker runtime using the configured image and bind paths.
- `apptainer`: Experimental path using similar bind/cache concepts to Docker. The image can be a `.sif` on disk or a `docker://` URI, depending on site policy.

Run and data paths are automatically bound (`<run_dir>:<run_dir>`, `<input_dir>:<input_dir>`); additional binds, container image URI, and cache/prefix can be provided per submission.

## API surface (runtime-aware)

- `POST /runs`: submit a run. Payload includes `workflow`, `run_name`, `runtime`, optional `container_image`, `bind_paths`, `cache_dir`, plus config/input fields.
- `GET /runs/{run_id}`: return run metadata (status, PID/return code, timestamps, command, paths).
- `GET /runs/{run_id}/logs`: tail the Snakemake log from `runtime/runs/<run>/logs/snakemake.log`.
- `POST /runs/{run_id}/cancel`: terminate the recorded PID and mark cancelled.
- `GET /runs/{run_id}/artifacts`: package the run directory to `runtime/artifacts/<run>.tar.gz`.
- Config + asset endpoints remain the same (`/configs/*`, `/classifiers`, `/adapters`).
- `GET /health`: liveness.

## Web UI

- Vue 3 + TypeScript SPA in `frontend/`, built with Vite.
- In the split deployment, Caddy serves the frontend and proxies `/api/*` to the FastAPI backend.
- In development, the Vite dev server runs on port 5173 and proxies API requests to port 8000.
- Adds a runtime selector (`conda`, `docker`, or `apptainer`) plus optional fields for container image URI, bind paths, and cache/prefix; shown only for containerized runtimes.
- Status/log/artifact polling works against the runtime-aware status/log endpoints.
- Upload/manage classifier and adapter assets from the browser.

### Using uploaded adapters and classifiers

- Uploaded adapters land at `runtime/adapters/<filename>`; classifiers land at `runtime/classifiers/<filename>` on the host. Because the repo is auto-bound to `/MetaWorks` for container runs, these become `/MetaWorks/runtime/adapters/<filename>` and `/MetaWorks/runtime/classifiers/<filename>` inside the container.
- In the Config cards:
  - **Adapters**: set `trimming.adapters` to `/MetaWorks/runtime/adapters/<your_adapter.fasta>` or another bound path if you stored it elsewhere and added a bind.
  - **Classifiers**: set `classification.rdp.use_custom_classifier` to `yes` and `classification.rdp.classifier_path` to `/MetaWorks/runtime/classifiers/<your_classifier.properties>` or your bound path.
- If you upload outside the repo, add a bind in "Bind paths" (for example, `/data/classifiers:/data/classifiers`) and reference the inside-container path you chose.

## Running on a single node

1. Ensure the host has Snakemake plus the chosen runtime available, or use the Docker Compose release stack where the backend image supplies the runtime.
2. Start the API + UI:
   ```bash
   python -m uvicorn api.main:app --host 0.0.0.0 --port 8000
   ```
3. Browse to `http://<host>:8000`.
4. Configure defaults via `METAWORKS_` env vars (see `api/config.py`), for example:
   - `METAWORKS_RUN_ROOT=/data/metaworks/runs`
   - `METAWORKS_CONTAINER_IMAGE=docker://metaworks:latest`
   - `METAWORKS_SINGULARITY_CACHE=/data/metaworks/cache`
   - `METAWORKS_DEFAULT_RUNTIME=apptainer` to favor Apptainer without editing the UI.
   - `METAWORKS_ALLOWED_RUNTIMES=docker,apptainer` to constrain web-submitted runtime modes.
   - `METAWORKS_RETENTION_POLICY=until_download` (`until_download|immediate|manual`) to control run file cleanup behavior.

## Future-proofing notes

- Snakemake remains the orchestrator; runner swaps (local `Popen` vs scheduler submit) stay behind the same API surface.
- Apptainer readiness is configuration-driven and site-dependent: provide image URI/path, bind paths, and cache dir through the API/UI fields or defaults, then validate on the target host.
- Conda remains supported by keeping `runtime=conda` and omitting container-specific fields.

## Deployment Quickstarts

For end-user Docker Compose setup, use [deploy/README.md](../deploy/README.md). For local source, server, and experimental HPC patterns, use [DEPLOYMENT_GUIDE.md](DEPLOYMENT_GUIDE.md). This document keeps the API/runtime contract separate from those operational walkthroughs.
