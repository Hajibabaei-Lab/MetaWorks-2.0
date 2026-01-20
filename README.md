# MetaWorks 2.0

A flexible, web-based control node for managing bioinformatics pipeline runs (ESV/OTU analysis) across local machines, servers, and HPC clusters.

## Features

- **Modern Web UI (work in progress) **: Vue 3 + Vite based interface with real-time updates
- **Multi-Environment Support**: Run locally, on servers, or on HPC clusters
- **Flexible Runtimes**: Conda, Docker, and Apptainer (Singularity) support
- **Workflow Presets**: Pre-configured templates for COI, 16S, and custom analyses
- **Smart Configuration**: Interactive config editor with validation and help tooltips
- **Real-Time Monitoring**: Live progress tracking, log streaming, and status updates
- **Asset Management**: Drag-and-drop upload for classifiers and adapters
- **Scheduler Integration**: Native SLURM support with extensible scheduler architecture

## Quick Start

### Prerequisites

- Python 3.9+
- Node.js 18+ (only for building UI - users don't need this)
- Conda/Mamba (optional, for local runs)
- Docker/Apptainer (optional, for containerized runs)

### Installation

1. **Clone and setup:**
   ```bash
   git clone https://github.com/Hajibabaei-Lab/MetaWorks-2.0.git
   cd MetaWorks-2.0
   ```

2. **Create conda environment:**
   ```bash
   conda env create -f environment.yml
   conda activate metaworks
   ```

3. **Build the UI (one-time setup):**
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

5. **Open in browser:**
   Navigate to `http://localhost:8000`

That's it! The web interface is now ready to use.

### Using Docker Compose (Easiest)

For a complete, containerized deployment:

```bash
docker-compose up -d
```

Access the UI at `http://localhost:8000`.

## Architecture

### Control Node Pattern

MetaWorks uses a **control node** architecture:

1. **Web UI**: Provides user interface for submitting and monitoring runs
2. **API Server**: FastAPI backend managing run lifecycle
3. **Job Manager**: Handles scheduler integration and job execution
4. **Runtime Layer**: Supports Conda, Docker, and Apptainer for running pipelines

This separation allows the control node to run anywhere (laptop, server, HPC) while the actual pipeline runs execute on appropriate compute resources.

### Zero-Download Experience

Users simply access the web interface - no frontend build steps required. The Vue UI is pre-built into static files served by FastAPI.

## Documentation

- **[Deployment Guide](docs/DEPLOYMENT_GUIDE.md)** - Deploy to local, server, or HPC environments
- **[Configuration Guide](docs/CONFIGURATION_GUIDE.md)** - System and user configuration
- **[Module Standards](docs/MODULE_STANDARDS.md)** - Creating custom modules
- **[Remote API Usage](docs/REMOTE_API_UI.md)** - Using the API programmatically

## Usage Overview

### Submitting a Run

1. **Choose a workflow preset** (COI Standard, 16S Microbiome, or Custom)
2. **Configure parameters**:
   - Runtime type (Conda, Docker, Apptainer)
   - Input directory and sample source
   - Resource requirements (cores, memory)
3. **Edit config sections** (optional):
   - Click "Load [ESV/OTU] sections" to see available parameters
   - Modify fields as needed - help tooltips explain each option
   - Only changed values are sent with the run
4. **Upload assets** (optional):
   - Drag-and-drop classifier and adapter files
   - Reference them in your config
5. **Submit to scheduler** and monitor progress

### Monitoring Runs

- **Auto-refresh**: Runs automatically refresh every 5 seconds
- **Progress tracking**: Percentage complete, current step, time estimates
- **Log streaming**: Real-time log output in the browser
- **Actions**: Cancel, download logs, download artifacts, delete runs

### Development Mode

For UI development with hot-reload:

```bash
# Terminal 1: Vite dev server
cd ui
npm run dev

# Terminal 2: API server
uvicorn api.main:app --host 0.0.0.0 --port 8000 --reload
```

## Project Structure

```
MetaWorks-2.0/
├── api/                    # FastAPI backend
│   ├── main.py             # Application entry point
│   ├── config.py           # Configuration management
│   ├── job_manager.py      # Job lifecycle management
│   ├── schemas.py         # Pydantic models
│   └── routes/            # API endpoints
│       ├── runs.py         # Run CRUD operations
│       ├── configs.py      # Configuration management
│       └── assets.py      # Asset uploads
├── ui/                    # Vue 3 + Vite frontend
│   ├── src/
│   │   ├── App.vue       # Main application
│   │   ├── main.js       # Entry point
│   │   └── components/   # Vue components
│   ├── dist/             # Build output (generated)
│   ├── package.json       # Node dependencies
│   └── vite.config.js    # Vite configuration
├── modules/               # Snakemake pipeline modules
├── runtime/              # Runtime directories (generated)
├── config/               # Configuration files
├── docs/                # Documentation
└── docker-compose.yml     # Docker deployment
```

## Deployment Options

### Local Development

Run on your laptop for testing and development:
```bash
conda activate metaworks
uvicorn api.main:app --host 0.0.0.0 --port 8000 --reload
```

### Production Server

Deploy to a server for team access:
```bash
docker-compose up -d
```

### HPC Cluster

Deploy on HPC with multiple options:
- **Dedicated control node** with shared storage
- **SSH tunneling** from local machine
- **Interactive job** on compute node
- **Reverse proxy** with authentication

See [Deployment Guide](docs/DEPLOYMENT_GUIDE.md) for detailed instructions.

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

See [Module Standards](docs/MODULE_STANDARDS.md) for guidance on creating new pipeline modules.

## License

[Specify your license here]

## Citation

If you use MetaWorks in your research, please cite:

```
[Add citation information here]
```

## Acknowledgments

- Hajibabaei Lab
- Contributors and community members

## Support

- Open an issue on GitHub
- Check [documentation](docs/)
- Contact the development team
