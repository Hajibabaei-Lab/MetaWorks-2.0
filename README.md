# MetaWorks 2.0

[![CI](https://github.com/Hajibabaei-Lab/MetaWorks-2.0/actions/workflows/ci.yml/badge.svg)](https://github.com/Hajibabaei-Lab/MetaWorks-2.0/actions/workflows/ci.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-2.0.0-blue.svg)](https://github.com/Hajibabaei-Lab/MetaWorks-2.0)

MetaWorks 2.0 is a modular, API-backed, web-operable metabarcoding platform for ESV and OTU workflows. It updates the original MetaWorks Snakemake pipeline for current research use by adding validated configuration profiles, reusable rule modules, a FastAPI control node, a Vue 3 web interface, and single-node container-oriented deployment.

## Statement of Need

Metabarcoding studies often need to run repeatable ESV/OTU analyses across marker genes, classifiers, runtimes, and computing environments. MetaWorks 2.0 is designed for biodiversity researchers, molecular ecology labs, and bioinformatics teams that need a reusable workflow for high-throughput multi-marker sequence processing without maintaining several separate ad hoc pipelines.

The v2.0.0 release fills the gap between the archived MetaWorks v1 command-line workflow and a maintainable research software platform. It keeps the established MetaWorks analysis concepts while replacing monolithic workflow files with a module registry, marker presets, validated configuration merging, Python helper scripts, API-based run control, and tested single-node Conda/Docker execution. Apptainer support is configuration-driven and intended for compatible deployments, but it has not yet been broadly validated across HPC systems.

## Features

- Modular Snakemake workflow for ESV and OTU metabarcoding analyses.
- Marker presets for COI, 16S, 12S, 18S, 28S, ITS, rbcL, and related use cases.
- Schema-driven configuration system for defaults, marker profiles, and user overrides.
- FastAPI backend for run submission, status, logs, artifacts, configuration rendering, and managed adapters/classifiers.
- Vue 3 + TypeScript web UI for submitting runs, monitoring progress, and managing runtime assets.
- Tested runtime support for host Conda and Docker, with experimental Apptainer/Singularity-oriented execution.
- **HPC integration (Apptainer/Singularity, SLURM, PBS) has not been tested yet.** Apptainer support is configuration-driven but has not been validated across HPC systems. Scheduler integration is planned for a future release.
- Split Docker Compose deployment with backend and frontend containers.
- SINTAX classifier support alongside RDP classifier workflows.
- Parallel RDP classifier wrapper (`parallel_rdp.py`) that splits input FASTA into chunks and runs multiple concurrent JVM processes for faster taxonomic assignment.
- Automated backend and frontend tests with GitHub Actions CI.

## Requirements

For local source installs:

- Python 3.9+
- Conda or Mamba
- Node.js 18+ and npm, only when building or developing the frontend
- Docker, optional for containerized runs
- Apptainer/Singularity, optional and experimental

For release-image deployment:

- Docker Desktop or Docker Engine with Compose support

## Installation

### Option A: Docker Compose Release Images

Use this path for the easiest end-user deployment once the v2.0.0 release images are available on GitHub Container Registry.

```bash
git clone https://github.com/Hajibabaei-Lab/MetaWorks-2.0.git
cd MetaWorks-2.0/deploy
cp .env.example .env
docker compose -f docker-compose.ghcr.yml pull
docker compose -f docker-compose.ghcr.yml up
```

Open `http://127.0.0.1:8080`. The bundled smoke-test data is available inside the backend image at `/MetaWorks/tests/testing_data`.

### Option B: Build from Source

```bash
git clone https://github.com/Hajibabaei-Lab/MetaWorks-2.0.git
cd MetaWorks-2.0
conda env create -f environment.yml
conda activate MetaWorks
cd frontend
npm install
npm run build
cd ..
uvicorn api.main:app --host 0.0.0.0 --port 8000
```

Open `http://localhost:8000`.

### Option C: Developer Split Stack

```bash
git clone https://github.com/Hajibabaei-Lab/MetaWorks-2.0.git
cd MetaWorks-2.0/deploy
cp .env.example .env
docker compose up --build
```

Open `http://127.0.0.1:8080`.

## Usage

### Web UI Run Submission

1. Choose a workflow: `ESV` or `OTU`.
2. Choose a marker profile, such as `coi`, `16s`, or `its`.
3. Select a runtime: Conda or Docker. Apptainer is available for experimental deployments.
4. Set the input directory and sample source.
5. Optionally upload adapter FASTA files or classifier assets.
6. Submit the run and monitor status, logs, and artifacts from the Runs view.

For the bundled test data, use:

```text
/MetaWorks/tests/testing_data
```

If needed, set adapter overrides to:

```json
{"trimming": {"adapters": "/MetaWorks/tests/adapters_anchored.fasta"}}
```

### Programmatic Configuration

```python
from lib.config.config_manager import ConfigManager

config = ConfigManager.load_from_dict(
    profile="coi",
    workflow="esv",
    user_overrides={"input": {"fastq_dir": "/data/my_samples"}},
    repo_root="/path/to/MetaWorks-2.0",
)

workflow_config = config.export_for_workflow()
```

### Direct Snakemake Use

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile config/defaults.yaml \
  --configfile config/presets/coi.yaml \
  --configfile my_run.yaml \
  --cores 4
```

## Documentation

- [Deployment Guide](docs/DEPLOYMENT_GUIDE.md): local, server, Docker Compose, and experimental HPC deployment notes.
- [Configuration Guide](docs/CONFIGURATION_GUIDE.md): profiles, configuration layers, and pipeline data flow.
- [Module Parameters](docs/MODULE_PARAMETERS.md): detailed parameter reference, complete examples, and common use cases.
- [Module Standards](docs/MODULE_STANDARDS.md): workflow module registry, rule conventions, and extension guidance.
- [Remote API Usage](docs/REMOTE_API_UI.md): FastAPI control node, runtime-aware endpoints, and UI/API behavior.
- [Docker Deployment README](deploy/README.md): end-user and developer Compose workflows.

For the API reference, interactive docs are available at `http://localhost:8000/api/docs` when the backend is running.

## Development and Tests

For development setup and local quality checks, see [CONTRIBUTING.md](CONTRIBUTING.md).

GitHub Actions runs backend lint/tests and frontend typecheck/tests/build on pushes and pull requests to `main` and `development`.

## Contributing

Contributions are welcome. Please read [CONTRIBUTING.md](CONTRIBUTING.md) before opening an issue or pull request. New workflow modules should follow [Module Standards](docs/MODULE_STANDARDS.md), and changes that affect behavior should include tests or documented verification steps.

## Support

- Open an issue on GitHub for bugs, feature requests, and support questions.
- Include the run ID, runtime mode, relevant logs, input/config summary, and MetaWorks version when reporting run failures.
- For deployment questions, start with [docs/DEPLOYMENT_GUIDE.md](docs/DEPLOYMENT_GUIDE.md) and [deploy/README.md](deploy/README.md).

## License

MetaWorks 2.0 is released under the GNU General Public License v3.0. See [LICENSE](LICENSE).

## Citation

If you use MetaWorks 2.0 in research, cite the archived software release:

Song, Y., & Hajibabaei, M. (2026). MetaWorks 2.0 (Version v2.0.0). Zenodo. https://doi.org/10.5281/zenodo.22210765

Citation metadata is also available in [CITATION.cff](CITATION.cff). Please also cite the MetaWorks paper:

Porter, T. M., & Hajibabaei, M. (2022). MetaWorks: A flexible, scalable bioinformatic pipeline for high-throughput multi-marker biodiversity assessments. PLOS ONE, 17(9), e0274260. https://doi.org/10.1371/journal.pone.0274260

The original v1 repository is archived as:

Porter, T. M. (2020). MetaWorks: A Multi-Marker Metabarcode Pipeline (Version v1.13.0). Zenodo. https://doi.org/10.5281/zenodo.4741407

Also cite the relevant method dependencies for your analysis:

- COI classifier: Porter, T. M., & Hajibabaei, M. (2018). Automated high throughput animal CO1 metabarcode classification. Scientific Reports, 8, 4226. https://doi.org/10.1038/s41598-018-22505-4
- Pseudogene filtering: Porter, T. M., & Hajibabaei, M. (2021). Profile hidden Markov model sequence analysis can help remove putative pseudogenes from DNA barcoding and metabarcoding datasets. BMC Bioinformatics, 22, 256. https://doi.org/10.1186/s12859-021-04180-x
- RDP classifier: Wang, Q., Garrity, G. M., Tiedje, J. M., & Cole, J. R. (2007). Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Applied and Environmental Microbiology, 73(16), 5261-5267. https://doi.org/10.1128/AEM.00062-07

Also cite the tools used in your analysis:

- Snakemake: Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., & Köster, J. (2021). Sustainable data analysis with Snakemake. F1000Research, 10, 33. https://doi.org/10.12688/f1000research.29032.2
- VSEARCH: Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4, e2584. https://doi.org/10.7717/peerj.2584
- Cutadapt: Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10-12. https://doi.org/10.14806/ej.17.1.200
- SeqPrep: St. John, J. (2016). SeqPrep. https://github.com/jstjohn/SeqPrep

## Acknowledgments

MetaWorks 2.0 is an updated and extended version of [MetaWorks v1](https://github.com/terrimporter/MetaWorks) by Terri Porter. Development is led by the Hajibabaei Lab with contributions from Yaye Song, Terri Porter and Mehrdad Hajibabaei.
