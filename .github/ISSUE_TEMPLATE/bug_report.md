---
name: Bug report
about: Report something that is not working as expected
title: "[bug] "
labels: bug
body:
  - type: textarea
    id: summary
    attributes:
      label: Summary
      description: What happened, and what did you expect to happen?
    validations:
      required: true
  - type: textarea
    id: repro
    attributes:
      label: Steps to reproduce
      placeholder: |
        1. Config used (marker preset, user_config.yaml, or Web UI settings)
        2. Command or UI action
        3. See error
    validations:
      required: true
  - type: textarea
    id: logs
    attributes:
      label: Relevant log output
      description: Snakemake error output or log file contents (e.g. output_dir/logs/...)
      render: shell
  - type: textarea
    id: environment
    attributes:
      label: Environment
      placeholder: |
        - MetaWorks version / commit:
        - Execution mode: conda / Docker / Apptainer, CLI or Web UI
        - OS:
    validations:
      required: true
