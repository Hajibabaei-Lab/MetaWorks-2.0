# Changelog

This changelog starts with the v2.0.0 release. Earlier MetaWorks v1 releases were
tracked through Git history, GitHub releases, and the original v1 Zenodo archive.

## v2.0.0

MetaWorks v2.0.0 is a major migration from the v1.13.0 command-line pipeline to a
modular, API-backed, web-operable platform for ESV and OTU metabarcoding workflows.

### Added

- Added a single modular Snakemake entry point at `workflow/Snakefile` with
  dynamically loaded rule modules for trimming, denoising, clustering, ITSx
  extraction, classification, pseudogene filtering, statistics, results
  consolidation, and global ESV/OTU aggregation.
- Added a module registry with dependency validation, topological include ordering,
  activation modes, terminal target resolution, and module metadata.
- Added a three-layer configuration system that merges defaults, marker presets,
  and user overrides into a validated frozen configuration object.
- Added schema generation for frontend-driven configuration forms.
- Added marker presets under `config/presets/` for supported marker workflows.
- Added a FastAPI backend with endpoints for health checks, settings, profiles,
  configuration schema/rendering, run lifecycle operations, logs, artifacts,
  adapters, and classifiers.
- Added a Vue 3, TypeScript, Vite, Element Plus, and TanStack Query frontend for
  submitting runs, monitoring runs, configuring pipeline settings, and managing
  runtime assets.
- Added JSON-backed run state management with file locking and process monitoring.
- Added Docker Compose deployment with separate backend and frontend services,
  plus development-oriented compose support.
- Added release automation for backend and frontend container images published to
  GitHub Container Registry.
- Added SINTAX classifier support alongside RDP.
- Added a parallel RDP classifier wrapper (`workflow/scripts/parallel_rdp.py`) that splits input FASTA files into chunks and runs concurrent `rdp_classifier` processes using a thread pool, reducing classification wall-clock time on multi-core systems.
- Added 166 automated tests covering API surface, configuration behavior, schema
  generation, and pipeline scripts.

### Changed

- Replaced six monolithic v1 Snakemake workflows with a composable rule-module
  architecture.
- Migrated Perl helper scripts to Python and consolidated duplicated abundance
  merging logic into a single implementation.
- Replaced flat, manually edited workflow YAMLs with validated module-oriented
  configuration.
- Standardized script command-line interfaces with `argparse`.
- Reworked output routing so ESV, OTU, ITSx, pseudogene, and global workflows are
  selected by module configuration instead of parse-time branching.
- Expanded runtime support from Conda-only execution to Conda, Docker, and
  Apptainer/Singularity-oriented operation.
- Moved web-submitted run outputs, staged assets, logs, and state into runtime
  directories managed by the backend.
- Marked trimming and denoising intermediates as temporary so Snakemake removes
  them once consumed, reducing peak disk usage during runs; per-sample trimmed
  FASTA files (`trimmed/{sample}.fasta.gz`) remain durable artifacts. Because
  intermediates are deleted, reruns after late-stage failures re-execute
  trimming from raw reads.

### Fixed

- Fixed nullable optional configuration schema metadata so the Configure page can
  load backend schema responses that include `null` values.
- Fixed append-mode output bugs that could duplicate partial results on reruns.
- Fixed marker coverage gaps by centralizing marker definitions used by taxonomy
  and abundance scripts.
- Added shell strict-mode handling across rule modules while preserving valid
  empty-result filtering behavior.
- Fixed FASTA header normalization so sequence contents are not modified.
- Fixed configuration routing issues affecting report type, output directory
  defaults, and marker-aware result headers.
- Fixed multiple script logic issues affecting HMM pseudogene filtering, minimum
  abundance thresholds, sample ordering, and global OTU result mapping.

### Known Issues

- `lib/exceptions.py` still defines a custom `RuntimeError` class that shadows the
  Python built-in name. This is deferred because current pipeline code does not
  trigger the conflicting import path.
