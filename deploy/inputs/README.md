Place user FASTQ input data here for the self-contained Docker deployment.

- Files in this folder are mounted into the backend container at `/data`
- In the web UI, use input paths like `/data/my-run`
- For the bundled smoke test data, use `/MetaWorks/tests/testing_data`
- A smoke test config is provided at `tests/testing_data/user_config.yaml` — it runs the full pipeline against the bundled COI test data using the built-in RDP classifier (no external database needed):
  ```
  snakemake --snakefile workflow/Snakefile --configfile tests/testing_data/user_config.yaml --cores 2
  ```

