# Control Node API and Web UI (single node)

This document captures the single-node control stack for MetaWorks: a FastAPI backend plus a static web UI that can launch Snakemake-based ESV/OTU runs using host conda, Docker, or Apptainer without changing the UI.

## Architecture

- **FastAPI service** in `api/`
  - Endpoints for run submission/status/log tail/artifact packaging/cancel.
  - Config templating (`/configs/defaults/{workflow}`, `/configs/render`).
  - Managed assets for classifiers/adapters (upload/list/delete).
  - Serves the static UI from `ui/`.
  - Tracks run metadata, PIDs/return codes, and log paths in `runtime/state/runs.json`.
- **Runner** (single node, pluggable): builds and launches Snakemake commands based on `runtime` (conda|docker|apptainer). Uses local `subprocess.Popen` today; can be swapped for a scheduler submitter later without UI changes.
- **State + data layout** (host-visible so each runtime sees the same files):
  - `runtime/runs/` per-run working directories (rendered configs, logs).
  - `runtime/artifacts/` packaged tarballs.
  - `runtime/classifiers/` and `runtime/adapters/` managed uploads.
  - `runtime/cache/` default Apptainer/Singularity cache/prefix (for container pulls).

## Runtime options (pluggable)

- `conda` – Runs Snakemake directly on the host:
  `snakemake -s <Snakefile> --configfile <rendered> --cores N --printshellcmds`
- `docker` – Runs Snakemake with `--use-singularity`, pointing to `docker://<image>`:
  `--singularity-args "--bind <runs>:<runs>,<input>:<input>,<extra>" --singularity-prefix <cache> --container-image docker://<image>`
- `apptainer` – Same flags as docker but the image can be a `.sif` on disk or a `docker://` URI. Uses the cache/prefix dir for pulls.

Run and data paths are automatically bound (`<run_dir>:<run_dir>`, `<input_dir>:<input_dir>`); additional binds, container image URI, and cache/prefix can be provided per submission.

## API surface (runtime-aware)

- `POST /runs` – submit a run. Payload includes `workflow`, `run_name`, `runtime`, optional `container_image`, `bind_paths`, `cache_dir`, plus config/input fields.
- `GET /runs/{run_id}` – return run metadata (status, PID/return code, timestamps, command, paths).
- `GET /runs/{run_id}/logs` – tail the Snakemake log from `runtime/runs/<run>/logs/snakemake.log`.
- `POST /runs/{run_id}/cancel` – SIGTERM the recorded PID and mark cancelled.
- `GET /runs/{run_id}/artifacts` – package the run dir to `runtime/artifacts/<run>.tar.gz`.
- Config + asset endpoints remain the same (`/configs/*`, `/classifiers`, `/adapters`).
- `GET /health` – liveness.

## Web UI

- Static HTML/JS/CSS in `ui/` served by FastAPI.
- Adds a runtime selector (conda|docker|apptainer) plus optional fields for container image URI, bind paths, and cache/prefix; shown only for containerized runtimes.
- Status/log/artifact polling is unchanged and works against the runtime-aware status/log endpoints.
- Upload/manage classifier and adapter assets from the browser.

### Using uploaded adapters and classifiers

- Uploaded adapters land at `runtime/adapters/<filename>`; classifiers land at `runtime/classifiers/<filename>` on the host. Because the repo is auto-bound to `/workspace` for container runs, these become `/workspace/runtime/adapters/<filename>` and `/workspace/runtime/classifiers/<filename>` inside the container.
- In the Config cards:
  - **Adapters**: set `CUTADAPT.fasta` to `/workspace/runtime/adapters/<your_adapter.fasta>` (or another bound path if you stored it elsewhere and added a bind).
  - **Classifiers**: set `RDP.custom` to `yes` and `RDP.t` to `/workspace/runtime/classifiers/<your_classifier.properties>` (or your bound path).
- If you upload outside the repo, add a bind in “Bind paths” (e.g., `/data/classifiers:/data/classifiers`) and reference the inside-container path you chose.

## Running on a single node

1. Ensure the host has Snakemake plus the chosen runtimes (conda env active; Docker daemon or Apptainer binary installed).
2. Start the API + UI:
   ```bash
   python -m uvicorn api.main:app --host 0.0.0.0 --port 8000
   ```
3. Browse to `http://<host>:8000`.
4. Configure defaults via `METAWORKS_` env vars (see `api/config.py`), e.g.:
   - `METAWORKS_RUN_ROOT=/data/metaworks/runs`
   - `METAWORKS_CONTAINER_IMAGE=docker://metaworks:latest`
   - `METAWORKS_SINGULARITY_CACHE=/data/metaworks/cache`
   - `METAWORKS_DEFAULT_RUNTIME=apptainer` to favor Apptainer without editing the UI.

## Future two-container deployment (plan)

- Split into **backend** (FastAPI + runner) and **frontend** (static UI via nginx/Caddy).
- Bind-mount into the backend container: `runtime/runs`, `runtime/artifacts`, `runtime/classifiers`, `runtime/adapters`, input data directories, and the Apptainer binary plus image (for `runtime=apptainer`).
- Expose backend on port `8000`; frontend serves static files and proxies to backend.
- Keep the same API/UI shape so local PIDs can later be replaced with scheduler IDs without changing the frontend.

## Future-proofing notes

- Snakemake remains the orchestrator; runner swaps (local `Popen` vs scheduler submit) stay behind the same API surface.
- Apptainer readiness is configuration-driven: provide image URI/path, bind paths, and cache dir through the API/UI fields or defaults.
- Conda remains supported by keeping `runtime=conda` and omitting container-specific fields.

## Local Docker quickstart (container drives the workflow)

Use this path when you already have the MetaWorks image locally and want Snakemake to run inside the container (no host Snakemake needed):

1. Prereqs
   - Docker installed and running.
   - MetaWorks repo checked out locally.
   - MetaWorks image available (e.g., `metaworks:latest`).

2. Start the API/UI from the repo root:
   ```bash
   python -m uvicorn api.main:app --host 0.0.0.0 --port 8000
   ```

3. Prepare paths (absolute):
   - Repo root: `/abs/path/to/MetaWorks-2.0` (mounted to `/workspace` inside the container).
   - Input data: e.g., `/abs/path/to/MetaWorks-2.0/testing/testing_data` (test FASTQs).
   - Optional adapters/classifiers: under `runtime/adapters` / `runtime/classifiers` or elsewhere (bind if outside the repo).
   - Cache/prefix (optional): `/abs/path/to/MetaWorks-2.0/runtime/cache`.

4. Submit via the UI (http://localhost:8000 → “Submit Run”):
   - Workflow: `ESV` (or OTU).
   - Runtime: `Docker`.
   - Container image URI: `metaworks:latest` (or your tag).
   - Bind paths (one per line, absolute host paths). Example on Windows:
     ```
     D:/MetaWorks-2.0:/workspace
     ```
     (The repo is auto-mounted to `/workspace`; the bundled test data lives at `/workspace/testing/testing_data` and works without extra binds. Add more binds only if your data/classifiers live outside the repo.)
   - Cache/prefix: `/abs/path/to/MetaWorks-2.0/runtime/cache` (or leave blank to use default).
   - Input directory: `/workspace/testing/testing_data` (bundled test data).
   - Samples CSV: leave blank for the bundled test.
   - Overrides JSON: optional; for the test adapters you can set
     ```json
     {"CUTADAPT": {"fasta": "/workspace/testing/adapters_anchored.fasta"}}
     ```
   - Leave cores/mem default; dry-run off.

5. Monitor
   - UI “Tail logs” shows `runtime/runs/<run_id>/logs/snakemake.log`.
   - Status updates include the exact `docker run ... snakemake ...` command.
   - Package artifacts when done → `runtime/artifacts/<run_id>.tar.gz`.

Notes
- The backend auto-binds the repo to `/workspace` inside the container so `Snakefile`/configs resolve. Any paths outside the repo must be bound explicitly and referenced by their inside-container path.
- Apptainer runtime follows the same pattern but uses `apptainer exec ... snakemake ...` with `-B` binds and the cache/prefix for pulls.
- Conda runtime runs Snakemake on the host if you prefer that mode.
