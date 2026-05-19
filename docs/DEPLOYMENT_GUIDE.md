# MetaWorks 2.0 Deployment Guide

This guide covers deploying MetaWorks 2.0 in various environments: local development, production servers, and HPC clusters.

## Recommended 2026 Deployment

MetaWorks supports a split deployment:

- `backend` container runs FastAPI and executes Snakemake pipelines via its own conda runtime
- `frontend` container builds the Vue 3 SPA and serves it via Caddy, proxying `/api/*` to the backend
- `deploy/docker-compose.yml` runs the split stack
- default runtime is `conda` inside the backend container, so users only need Docker/Compose on the host

## Table of Contents

- [Quick Start](#quick-start)
- [Local Development](#local-development)
- [Server Deployment](#server-deployment)
- [HPC Deployment](#hpc-deployment)
- [Building the UI](#building-the-ui)
- [Configuration](#configuration)
- [Troubleshooting](#troubleshooting)

## Quick Start

### Prerequisites

- Python 3.9+
- Node.js 18+ (for building the UI)
- Docker (optional, for containerized runs)
- Conda/Mamba (for local runs)

### Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/Hajibabaei-Lab/MetaWorks-2.0.git
   cd MetaWorks-2.0
   ```

2. **Create conda environment:**
   ```bash
   conda env create -f environment.yml
   conda activate MetaWorks
   ```

3. **Build the UI:**
   ```bash
   cd frontend
   npm install
   npm run build
   cd ..
   ```

4. **Start the API server:**
    ```bash
    uvicorn api.main:app --host 0.0.0.0 --port 8000
    ```

5. **Open the web UI:**
   Navigate to `http://localhost:8000` in your browser.

## Local Development

### Using Conda Runtime

For development and testing on a local machine with conda:

```bash
# Activate the MetaWorks environment
conda activate MetaWorks

# Start the API server
uvicorn api.main:app --host 0.0.0.0 --port 8000 --reload
```

The `--reload` flag enables hot-reloading during development.

### Using Docker Runtime

To use Docker for running the pipeline while running the API server locally:

1. **Build the MetaWorks Docker image:**
   ```bash
   docker build -t metaworks:latest .
   ```

2. **Start the API server:**
   ```bash
   uvicorn api.main:app --host 0.0.0.0 --port 8000
   ```

3. **Configure runs in the web UI:**
   - Set Runtime to "Docker"
   - Container image: `docker://metaworks:latest`
   - Add bind paths as needed (e.g., `/path/to/your/data:/MetaWorks/input`)

### Development Mode

For UI development with hot-reloading:

```bash
# Terminal 1: Start Vite dev server
cd frontend
npm run dev

# Terminal 2: Start API server
uvicorn api.main:app --host 0.0.0.0 --port 8000 --reload
```

The Vite dev server runs on port 5173 and proxies API requests to port 8000.

## Server Deployment

### Using Docker Compose (Recommended)

Docker Compose provides the easiest way to deploy to a server:

1. **Build and start:**
   ```bash
   docker compose -f deploy/docker-compose.yml up -d --build
   ```

2. **Configure environment variables:**
   Adjust environment variables in `deploy/docker-compose.yml` or via a `.env` file next to it.

3. **Access the UI:**
   Navigate to `http://your-server:8080` (Caddy serves the frontend and proxies `/api/*` to the backend on port 8000).

### Manual Deployment

For more control over the deployment:

1. **Install dependencies:**
   ```bash
   conda env create -f environment.yml
   conda activate MetaWorks
   ```

2. **Build the UI:**
   ```bash
   cd frontend
   npm install
   npm run build
   cd ..
   ```

3. **Configure environment variables:**
   See the full environment variable table below for all `METAWORKS_*` variables.

4. **Start the API server:**
   ```bash
   uvicorn api.main:app --host 0.0.0.0 --port 8000
   ```

5. **Serve the frontend:**
   Use a reverse proxy (Caddy, Nginx) to serve the built frontend and proxy `/api/*` to port 8000.

## HPC Deployment

### Control Node on HPC

Many HPC clusters don't allow running web servers on login nodes. You have several options:

#### Option 1: Dedicated Control Node

Deploy MetaWorks on a separate server or virtual machine that has network access to the HPC cluster:

1. **Deploy on a separate server** following the [Server Deployment](#server-deployment) instructions.
2. **Configure shared storage** between the control node and HPC:
   - Use NFS, Lustre, or similar shared filesystem
   - Mount the same paths on both systems
3. **Configure bind paths** to point to shared storage

#### Option 2: SSH Tunneling

Run MetaWorks on the login node and access via SSH tunnel:

1. **Start MetaWorks on the HPC login node:**
   ```bash
   module load conda
   conda activate MetaWorks
   uvicorn api.main:app --host 127.0.0.1 --port 8000
   ```

2. **Create SSH tunnel from your local machine:**
   ```bash
   ssh -L 8000:localhost:8000 user@hpc-cluster.edu
   ```

3. **Access the UI** at `http://localhost:8000` in your browser.

#### Option 3: Interactive Job

Run MetaWorks within an HPC interactive job:

1. **Request an interactive node:**
   ```bash
   srun --pty --nodes=1 --cpus-per-task=4 --mem=16G bash
   ```

2. **Start MetaWorks** as above.

#### Option 4: Reverse Proxy with Apache/Nginx

If the cluster has a web server:

1. **Configure reverse proxy** to forward requests to MetaWorks running on a compute node
2. **Use appropriate authentication** (e.g., .htaccess for Apache, auth_basic for Nginx)

### HPC-Specific Configuration

For HPC deployment, use Apptainer (Singularity) runtime:

```bash
# Build Singularity image
sudo singularity build metaworks.sif docker://metaworks:latest

# Use in MetaWorks UI:
# Runtime: Apptainer
# Container image: /path/to/metaworks.sif
# Bind paths: /shared/storage:/MetaWorks/storage
```

### Scheduler Integration

Scheduler integration (SLURM, PBS) is planned for a future release. The current architecture uses local `subprocess.Popen` for pipeline execution. When implemented, the API surface will remain the same — only the runner backend changes.

## Building the UI

### Prerequisites

- Node.js 18+ and npm

### Build Process

1. **Install dependencies:**
   ```bash
   cd frontend
   npm install
   ```

2. **Build for production:**
   ```bash
   npm run build
   ```
   This creates the `frontend/dist` directory with optimized static files.

3. **Build for development (hot-reload):**
   ```bash
   npm run dev
   ```
   This starts a dev server on port 5173 with hot-reload enabled.

4. **Preview production build:**
   ```bash
   npm run preview
   ```

## Configuration

### Environment Variables

MetaWorks uses environment variables (prefix `METAWORKS_`) for all runtime configuration. These are loaded by `api/config.py` via pydantic-settings.

| Variable | Description | Default |
|----------|-------------|---------|
| `METAWORKS_REPO_ROOT` | Path to MetaWorks repository | Auto-detected from `api/` parent |
| `METAWORKS_RUN_ROOT` | Directory for run outputs | `<repo>/runtime/runs` |
| `METAWORKS_STAGING_ROOT` | Directory for staging uploads | `<repo>/runtime/staging` |
| `METAWORKS_ARTIFACT_ROOT` | Directory for artifacts | `<repo>/runtime/artifacts` |
| `METAWORKS_CLASSIFIER_ROOT` | Directory for classifier files | `<repo>/runtime/classifiers` |
| `METAWORKS_ADAPTER_ROOT` | Directory for adapter FASTA files | `<repo>/runtime/adapters` |
| `METAWORKS_STATE_FILE` | Path to run state JSON | `<repo>/runtime/state/runs.json` |
| `METAWORKS_SNAKEFILES` | Workflow entrypoints (JSON) | `{"esv":"workflow/Snakefile","otu":"workflow/Snakefile"}` |
| `METAWORKS_DEFAULT_RUNTIME` | Default runtime (conda/docker/apptainer) | `docker` |
| `METAWORKS_ALLOWED_RUNTIMES` | Comma-separated allowed runtimes | `docker,apptainer` |
| `METAWORKS_RETENTION_POLICY` | Run data retention mode | `until_download` |
| `METAWORKS_CONTAINER_IMAGE` | Default container image URI | `docker://metaworks:latest` |
| `METAWORKS_BIND_PATHS` | Additional bind mounts (JSON list) | `[]` |
| `METAWORKS_SINGULARITY_CACHE` | Apptainer/Singularity cache directory | `<repo>/runtime/cache` |
| `METAWORKS_DEFAULT_CORES` | Cores requested per job | `32` |
| `METAWORKS_LOG_TAIL_LINES` | Default lines for log tail endpoint | `200` |
| `METAWORKS_SERVE_LEGACY_UI` | Serve legacy static UI from backend | `true` |
| `METAWORKS_CORS_ALLOWED_ORIGINS` | Comma-separated allowed CORS origins | `http://localhost:5173,...` |

## Troubleshooting

### UI Not Loading

**Problem:** Accessing the UI shows old content or 404

**Solution:** Build the Vue UI:
```bash
cd frontend
npm install
npm run build
```

### Port Already in Use

**Problem:** `OSError: [Errno 48] Address already in use`

**Solution:** Use a different port:
```bash
uvicorn api.main:app --host 0.0.0.0 --port 8080
```

### Permission Denied on Runtime Directories

**Problem:** Cannot write to runtime directories

**Solution:** Ensure proper permissions:
```bash
chmod -R 755 runtime/
chmod -R 777 runtime/runs
```

### Docker Runtime Issues

**Problem:** Docker container cannot access files

**Solution:** Check bind paths and ensure:
1. Host paths exist
2. Bind paths use correct format: `host:/container`
3. Container has read permissions

### HPC Scheduler Issues

**Problem:** Jobs not submitting to scheduler

**Solution:** Scheduler integration (SLURM, PBS) is planned for a future release. The current implementation uses local subprocess execution only.

### Build Errors

**Problem:** `npm install` fails

**Solution:** Clear cache and retry:
```bash
cd frontend
rm -rf node_modules package-lock.json
npm install
```

### CORS Errors in Development

**Problem:** Browser shows CORS errors

**Solution:** The FastAPI app allows all origins by default. If you've modified it, ensure CORS middleware is configured:
```python
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

## Additional Resources

- [Configuration Guide](CONFIGURATION_GUIDE.md)
- [Module Standards](MODULE_STANDARDS.md)
- [Remote API Usage](REMOTE_API_UI.md)

## Support

For issues and questions:
- Open an issue on GitHub
- Contact the development team
- Check the documentation in the `docs/` directory
