# Self-Contained Split Deployment

This deployment is optimized for the lowest-friction user experience:

- `frontend` serves the standalone Vue UI
- `backend` runs FastAPI and executes both MetaWorks pipelines from its own Conda environment
- runtime state persists in a Docker named volume
- users only need Docker and Docker Compose

## What Users Need To Install

- Docker Desktop or Docker Engine with Compose support

They do **not** need to install Python, Conda, Snakemake, FastAPI, or the MetaWorks pipeline
dependencies on the host.

## Input Data

- Quick test dataset: already bundled in the backend image at `/MetaWorks/tests/testing_data`
- User data: place FASTQ files under `deploy/inputs/` and use `/data/...` paths in the UI

The `deploy/inputs` directory is bind-mounted read-only into the backend container.

## Start

```bash
cd deploy
cp .env.example .env
docker compose up --build
```

Open `http://localhost:8080`.

If the standalone UI repo is not checked out as `../../metaworks-ui`, update
`METAWORKS_UI_CONTEXT` in `.env`.

## SSH-Only Server Access

For a lab server, the safest default is to keep the web UI bound to server-localhost only:

```env
METAWORKS_FRONTEND_BIND=127.0.0.1:8080:8080
```

That means the app is not directly reachable from the network. Users access it through SSH port
forwarding:

```bash
ssh -L 8080:127.0.0.1:8080 <username>@<lab-server>
```

Then they open `http://127.0.0.1:8080` in their local browser.

This works well when:

- users already have SSH accounts on the server
- you want SSH itself to be the access gate
- you do not want to expose the UI publicly

If you later want the UI reachable without an SSH tunnel, change `METAWORKS_FRONTEND_BIND` to
`8080:8080` or put Caddy behind your institution's existing reverse proxy/auth layer.

## Development Override

Expose the backend directly on port `8000` while still running the split stack:

```bash
docker compose -f docker-compose.yml -f docker-compose.dev.yml up --build
```

## Notes

- The default runtime in this deployment is `conda`, so workflow execution stays inside the backend
  container and avoids host Docker socket/path issues.
- Runtime state, uploads, logs, and artifacts are persisted in the `metaworks-runtime` named
  volume.
- The backend remains independently usable through `uvicorn api.main:app` or direct `/api/*`
  access outside Compose if you need advanced runtimes later.
