# MetaWorks 2.0 Deployment Guide

MetaWorks runs as two Docker containers:

- **backend** — Conda environment with the Snakemake pipeline, Python scripts, and FastAPI API.
  Defaults to an interactive bash shell. When run via Compose, it starts the API server instead.
- **frontend** — Caddy web server hosting the Vue UI. Proxies `/api/*` requests to the backend.

Runtime state persists in a Docker named volume. Users only need Docker and Docker Compose.

## Prerequisites

- Docker Desktop or Docker Engine with Compose support

You do **not** need Python, Conda, Snakemake, FastAPI, or any pipeline dependencies on the host.

## Choosing a Compose File

| File | Purpose | Requires |
|------|---------|-----------|
| `docker-compose.ghcr.yml` | **End users** — pulls pre-built images from GitHub | Docker only |
| `docker-compose.yml` | **Developers** — builds images from local source code | Cloned repo |

Both files produce the same running stack. The ghcr file is faster to start (no build) but
requires a GitHub Release to have been published. The build file lets you test local changes.

### Option A: Pre-built Images (recommended for users)

```bash
cd deploy
cp .env.example .env
docker compose -f docker-compose.ghcr.yml pull
docker compose -f docker-compose.ghcr.yml up
```

This pulls `ghcr.io/hajibabaei-lab/metaworks-backend` and
`ghcr.io/hajibabaei-lab/metaworks-frontend` from GitHub Container Registry. No source code needed.

To pin a specific release version, edit the `image:` tags in the compose file (e.g.
`metaworks-backend:v2.0.0` instead of `metaworks-backend:latest`).

### Option B: Build from Source (for developers)

```bash
git clone https://github.com/Hajibabaei-Lab/MetaWorks-2.0.git
cd MetaWorks-2.0/deploy
cp .env.example .env
docker compose up --build
```

This builds both images locally from the `Dockerfile` and `frontend/Dockerfile`. Use this when
you are modifying the pipeline, API, or frontend.

### Development Override

Expose the backend API directly on port 8000 (useful for debugging with `curl` or Swagger UI):

```bash
docker compose -f docker-compose.yml -f docker-compose.dev.yml up --build
```

The API docs are then available at `http://localhost:8000/api/docs`.

## Input Data

- **Quick test**: a small COI dataset is bundled in the backend image at
  `/MetaWorks/tests/testing_data`
- **Your data**: place FASTQ files under `deploy/inputs/` and reference them as `/data/...` paths

The `deploy/inputs` directory is bind-mounted read-only into the backend container at `/data`.

## CLI-Only Pipeline Use (No Web UI)

The backend image defaults to an interactive bash shell. You can run Snakemake directly without
the frontend or API:

```bash
docker run -it --rm \
  -v $(pwd)/my-data:/data \
  ghcr.io/hajibabaei-lab/metaworks-backend:latest

# Inside the container:
snakemake --snakefile workflow/Snakefile \
  --configfile tests/testing_data/user_config.yaml \
  --cores 4
```

The image includes `nano` and `vim` for editing config files.

## Network Access

### Local machine

Open `http://localhost:8080` after starting the stack.

### SSH-Only Server (lab server)

The default binds the frontend to localhost only, so it is not directly reachable from the
network. Users access it through SSH port forwarding:

```bash
ssh -L 8080:127.0.0.1:8080 <username>@<lab-server>
```

Then open `http://127.0.0.1:8080` in the local browser.

This works well when users already have SSH accounts and you want SSH to be the access gate.

To expose the UI without SSH, change `METAWORKS_FRONTEND_BIND` in `.env`:

```env
METAWORKS_FRONTEND_BIND=8080:8080
```

Or put Caddy behind your institution's existing reverse proxy/auth layer.

## Configuration

Copy and edit the environment file:

```bash
cp .env.example .env
```

Available variables:

| Variable | Default | Description |
|----------|---------|-------------|
| `METAWORKS_FRONTEND_BIND` | `127.0.0.1:8080:8080` | Host:port binding for the frontend |
| `METAWORKS_DEFAULT_RUNTIME` | `conda` | Default Snakemake runtime. In release Compose deployments this means the backend container's bundled Conda environment, not host Conda. |
| `METAWORKS_ALLOWED_RUNTIMES` | `conda` | Comma-separated list of allowed Snakemake runtimes |
| `METAWORKS_CORS_ALLOWED_ORIGINS` | `http://localhost:5173,...` | Allowed CORS origins for the API |

## Architecture

```
Browser → :8080 (frontend/Caddy)
              ├── /api/* → backend:8000 (FastAPI)
              └── /*     → Vue SPA (/srv)

Backend container:
  /MetaWorks/           — pipeline code, scripts, config
  /MetaWorks/runtime/   — named volume (classifiers, runs, state, artifacts)
  /data/                — bind mount from deploy/inputs/ (read-only)
```

## Notes

- The default runtime is `conda`, so pipeline execution stays inside the backend container's bundled MetaWorks Conda environment. Users do not need Conda or Snakemake on the host.
- Runtime state, uploads, logs, and artifacts persist in the `metaworks-runtime` named volume.
- The backend is independently usable through `docker run -it ... bash` for CLI-only pipeline work.
