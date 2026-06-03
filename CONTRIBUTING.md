# Contributing to MetaWorks 2.0

Thank you for helping improve MetaWorks. This guide describes how to report issues, propose changes, and prepare pull requests for the v2.x codebase.

## Support and Issue Reports

Use GitHub issues for bug reports, feature requests, documentation gaps, and support questions. For run failures, include:

- MetaWorks version or commit SHA.
- Runtime mode: Conda, Docker, or Apptainer.
- Operating system or server/HPC context.
- Workflow and profile, such as `ESV` with `coi`.
- Input layout summary, without private sample data.
- Relevant run ID, log excerpts, and rendered config if available.
- Steps needed to reproduce the problem.

## Development Setup

```bash
git clone https://github.com/Hajibabaei-Lab/MetaWorks-2.0.git
cd MetaWorks-2.0
conda env create -f environment.yml
conda activate MetaWorks
cd frontend
npm install
cd ..
```

For frontend development, run the backend and Vite dev server in separate terminals:

```bash
uvicorn api.main:app --host 0.0.0.0 --port 8000 --reload
cd frontend && npm run dev
```

## Quality Checks

Run the full local check set before opening a pull request:

```bash
make lint
make test
make build
```

Equivalent individual commands:

```bash
ruff check .
PYTHONPATH=$(pwd) pytest -v
cd frontend && npm run typecheck
cd frontend && npm run test
cd frontend && npm run build
```

If a check cannot be run locally, state that in the pull request and explain why.

## Pull Requests

Pull requests should:

- Describe the user-facing behavior or research workflow affected.
- Keep changes focused to the relevant API, frontend, workflow, docs, or config area.
- Include tests for behavior changes, or document an objective manual verification path.
- Update docs when changing configuration, deployment, APIs, or workflow modules.
- Avoid committing local runtime outputs, artifacts, caches, or generated run data.

## Workflow Module Changes

New or modified Snakemake modules should follow [docs/MODULE_STANDARDS.md](docs/MODULE_STANDARDS.md). In particular:

- Add or update the module registry entry in `lib/config/module_registry.py`.
- Add default configuration and validation for new parameters.
- Keep rule names lowercase with underscores.
- Include log directives for rules.
- Use write-only outputs rather than append-mode outputs.
- Add focused tests for module-selection or script behavior where practical.

## Release Contributions

Release-facing changes should update the relevant metadata:

- `CHANGELOG.md` for notable user-facing changes.
- `CITATION.cff` and `.zenodo.json` for citation/archive metadata changes.
- Version fields in Python, API, frontend, and config files when preparing a tagged release.

## AI Assistance Disclosure

If generative AI tools are used to draft code, docs, tests, or release materials, disclose the tool and scope of assistance in the pull request. Human contributors remain responsible for reviewing, validating, and licensing all submitted material.
