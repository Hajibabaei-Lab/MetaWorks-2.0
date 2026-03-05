# MetaWorks 2.0 Deployment Guide

This guide covers deploying MetaWorks 2.0 in various environments: local development, production servers, and HPC clusters.

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
   conda activate metaworks
   ```

3. **Build the UI:**
   ```bash
   cd ui
   npm install
   npm run build
   cd ..
   ```

4. **Start the server:**
   ```bash
   uvicorn api.main:app --host 0.0.0.0 --port 8000
   ```

5. **Open the web UI:**
   Navigate to `http://localhost:8000` in your browser.

## Local Development

### Using Conda Runtime

For development and testing on a local machine with conda:

```bash
# Activate the metaworks environment
conda activate metaworks

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
   - Add bind paths as needed (e.g., `/path/to/your/data:/workspace/input`)

### Development Mode

For UI development with hot-reloading:

```bash
# Terminal 1: Start Vite dev server
cd ui
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
   docker-compose up -d --build
   ```

2. **Configure volumes:**
   Edit `docker-compose.yml` to mount your data directories:
   ```yaml
   volumes:
     - ./runtime:/workspace/runtime
     - /path/to/your/data:/workspace/input:ro
     - ./config:/workspace/config:ro
   ```

3. **Configure environment variables:**
   Adjust environment variables in `docker-compose.yml` as needed.

4. **Access the UI:**
   Navigate to `http://your-server:8000`

### Manual Deployment

For more control over the deployment:

1. **Install dependencies:**
   ```bash
   conda env create -f environment.yml
   conda activate metaworks
   ```

2. **Build the UI:**
   ```bash
   cd ui
   npm install
   npm run build
   cd ..
   ```

3. **Configure environment variables:**
   Create a `.env` file:
   ```bash
   METAWORKS_REPO_ROOT=/path/to/MetaWorks-2.0
   METAWORKS_RUN_ROOT=/path/to/runtime/runs
   METAWORKS_ARTIFACT_ROOT=/path/to/runtime/artifacts
   METAWORKS_ADAPTER_ROOT=/path/to/runtime/adapters
   METAWORKS_CLASSIFIER_ROOT=/path/to/runtime/classifiers
   METAWORKS_STAGING_ROOT=/path/to/runtime/staging
   METAWORKS_SINGULARITY_CACHE=/path/to/runtime/cache
   METAWORKS_DEFAULT_RUNTIME=docker
   METAWORKS_CONTAINER_IMAGE=docker://metaworks:latest
   ```

4. **Start with a production server:**
   ```bash
   gunicorn api.main:app -w 4 -k uvicorn.workers.UvicornWorker --bind 0.0.0.0:8000
   ```

5. **Set up a reverse proxy (optional but recommended):**

   **Nginx example:**
   ```nginx
   server {
       listen 80;
       server_name metaworks.example.com;

       location / {
           proxy_pass http://localhost:8000;
           proxy_set_header Host $host;
           proxy_set_header X-Real-IP $remote_addr;
       }
   }
   ```

### Using systemd (Linux)

Create a systemd service for automatic startup:

```ini
# /etc/systemd/system/metaworks.service
[Unit]
Description=MetaWorks Control Node
After=network.target

[Service]
Type=simple
User=metaworks
WorkingDirectory=/path/to/MetaWorks-2.0
Environment="PATH=/path/to/conda/envs/metaworks/bin"
ExecStart=/path/to/conda/envs/metaworks/bin/uvicorn api.main:app --host 0.0.0.0 --port 8000
Restart=always

[Install]
WantedBy=multi-user.target
```

Enable and start the service:
```bash
sudo systemctl enable metaworks
sudo systemctl start metaworks
```

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
   conda activate metaworks
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
# Bind paths: /shared/storage:/workspace/storage
```

### Scheduler Integration

MetaWorks integrates with SLURM scheduler out of the box. Configure in the web UI:

1. **Set Runtime** to "Docker" or "Apptainer"
2. **Configure Scheduler Settings** in system config:
   ```yaml
   scheduler:
     type: slurm
     partition: compute
     account: your-account
     time: "24:00:00"
     cpus-per-task: 8
     memory: 32G
   ```

## Building the UI

### Prerequisites

- Node.js 18+ and npm

### Build Process

1. **Install dependencies:**
   ```bash
   cd ui
   npm install
   ```

2. **Build for production:**
   ```bash
   npm run build
   ```
   This creates the `ui/dist` directory with optimized static files.

3. **Build for development (hot-reload):**
   ```bash
   npm run dev
   ```
   This starts a dev server on port 5173 with hot-reload enabled.

4. **Preview production build:**
   ```bash
   npm run preview
   ```

### Zero-Download Deployment

The Vue UI is built into static files that are served by FastAPI. Users don't need to install Node.js or run build commands - they simply access the web interface.

## Configuration

### Environment Variables

MetaWorks uses environment variables for configuration:

| Variable | Description | Default |
|----------|-------------|---------|
| `METAWORKS_REPO_ROOT` | Path to MetaWorks repository | Auto-detected |
| `METAWORKS_RUN_ROOT` | Directory for run outputs | `runtime/runs` |
| `METAWORKS_ARTIFACT_ROOT` | Directory for artifacts | `runtime/artifacts` |
| `METAWORKS_ADAPTER_ROOT` | Directory for adapter files | `runtime/adapters` |
| `METAWORKS_CLASSIFIER_ROOT` | Directory for classifiers | `runtime/classifiers` |
| `METAWORKS_STAGING_ROOT` | Directory for staging | `runtime/staging` |
| `METAWORKS_SINGULARITY_CACHE` | Singularity cache directory | `runtime/cache` |
| `METAWORKS_DEFAULT_RUNTIME` | Default runtime (conda/docker/apptainer) | `docker` |
| `METAWORKS_CONTAINER_IMAGE` | Default container image | None |

### System Configuration

Edit `config/system_config.yaml` for system-wide settings:

```yaml
# Scheduler configuration
scheduler:
  type: slurm  # slurm, pbs, local
  partition: compute
  account: your-account

# Resource defaults
resources:
  default_cores: 8
  default_memory: 32G
  default_time: "24:00:00"

# API settings
api:
  max_runs: 100
  log_retention_days: 30
```

### User Configuration

Users can override defaults via `config/user_config.yaml` or through the web UI.

## Troubleshooting

### UI Not Loading

**Problem:** Accessing `http://localhost:8000` shows old UI or 404

**Solution:** Build the Vue UI:
```bash
cd ui
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

**Solution:** Check:
1. Scheduler is configured in system config
2. Account and partition are valid
3. You have proper cluster permissions
4. Scheduler commands (sbatch, qsub) are available

### Build Errors

**Problem:** `npm install` fails

**Solution:** Clear cache and retry:
```bash
cd ui
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
