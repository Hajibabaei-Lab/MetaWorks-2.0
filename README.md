# MetaWorks 2.0

A flexible, web-based control node for managing bioinformatics pipeline runs (ESV/OTU analysis) across local machines, servers, and HPC clusters.

## Features

- **Modern Web UI (work in progress) **: Vue 3 + Vite based interface with real-time updates
- **Multi-Environment Support**: Run locally, on servers, or on HPC clusters
- **Flexible Runtimes**: Conda, Docker, and Apptainer (Singularity) support
- **Workflow Presets**: Pre-configured templates for COI, 16S, and custom analyses
- **Schema-Driven Configuration**: Dynamic config forms generated from the pipeline schema
- **Real-Time Monitoring**: Live progress tracking, log streaming, and status updates
- **Scheduler Integration**: Pluggable scheduler architecture (SLURM planned)

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
   conda activate MetaWorks
   ```

3. **Build the UI (one-time setup):**
   ```bash
   cd frontend
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

### Using Docker Compose

```bash
cd deploy
cp .env.example .env
docker compose up --build
```

Access the UI at `http://localhost:8080`.

## Architecture

### Control Node Pattern

MetaWorks uses a **control node** architecture:

1. **Web UI**: Vue 3 SPA served from `frontend/`
2. **API Server**: FastAPI backend managing run lifecycle
3. **Job Manager**: Handles scheduler integration and job execution
4. **Runtime Layer**: Supports Conda, Docker, and Apptainer for running pipelines

This separation allows the control node to run anywhere (laptop, server, HPC) while the actual pipeline runs execute on appropriate compute resources.

### Split Frontend Deployment

The backend now exposes a stable `/api` surface and optional legacy static UI serving. The
recommended deployment runs the standalone frontend in its own container and proxies `/api/*` to
FastAPI, so the runner remains independently usable without the web app.

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
    - Upload classifier and adapter files via the UI
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
cd frontend
npm run dev

# Terminal 2: API server
uvicorn api.main:app --host 0.0.0.0 --port 8000 --reload
```

## Project Structure

```
MetaWorks-2.0/
├── api/                    # FastAPI backend (routes, services, schemas)
├── config/                 # Pipeline defaults and marker presets
├── deploy/                 # Docker Compose split deployment
├── docs/                   # Documentation
├── frontend/               # Vue 3 + TypeScript SPA
├── lib/                    # Config management, runtime builders, exceptions
├── tests/                  # pytest test suite (162 tests)
├── workflow/               # Snakemake pipeline (rules/, scripts/, profiles/)
├── Makefile                # Dev, test, lint, build commands
└── environment.yml         # Conda environment
```

## Deployment Options

### Local Development

Run on your laptop for testing and development:
```bash
 conda activate MetaWorks
 uvicorn api.main:app --host 0.0.0.0 --port 8000 --reload
```

### Production Server

Deploy the recommended split stack:
```bash
cd deploy
docker compose up -d --build
```

For the quickest smoke test after startup, submit a run against
`/MetaWorks/tests/testing_data`, which is already bundled into the backend image.

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

GNU General Public License v3.0

## Origin

MetaWorks 2.0 is an updated and extended version of [MetaWorks v1](https://github.com/terrimporter/MetaWorks) by Terri Porter, with new features, a modular Snakemake architecture, a FastAPI backend, a Vue 3 frontend, and Perl scripts migrated to Python.

## Citation

If you use MetaWorks in your research, please cite the MetaWorks paper:
Porter, T. M., & Hajibabaei, M. (2022). MetaWorks: A flexible, scalable bioinformatic pipeline for high-throughput multi-marker biodiversity assessments. PLOS ONE, 17(9), e0274260. doi: 10.1371/journal.pone.0274260

The original v1 repository is archived at: Teresita M. Porter. (2020, June 25). MetaWorks: A Multi-Marker Metabarcode Pipeline (Version v1.13.0). Zenodo. http://doi.org/10.5281/zenodo.4741407

If you use this dataflow for making COI taxonomic assignments, please cite the COI classifier publication:
Porter, T. M., & Hajibabaei, M. (2018). Automated high throughput animal CO1 metabarcode classification. Scientific Reports, 8, 4226.

If you use the pseudogene filtering methods, please cite the pseudogene publication: Porter, T.M., & Hajibabaei, M. (2021). Profile hidden Markov model sequence analysis can help remove putative pseudogenes from DNA barcoding and metabarcoding datasets. BMC Bioinformatics, 22: 256.

If you use the RDP classifier, please cite the publication:
Wang, Q., Garrity, G. M., Tiedje, J. M., & Cole, J. R. (2007). Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Applied and Environmental Microbiology, 73(16), 5261–5267. doi:10.1128/AEM.00062-07

Last updated: May 2026

## Acknowledgments

- Hajibabaei Lab
- Terri Porter
- Alex Song
- Contributors and community members

## Support

- Open an issue on GitHub
- Check [documentation](docs/)
- Contact the development team
