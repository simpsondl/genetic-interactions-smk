# SMK_PIPELINE

This repository contains a Snakemake-based pipeline for processing genetic interaction (GI) screens. It coordinates a set of R scripts (under `workflow/scripts/`) and Snakemake rules (under `workflow/rules/`) to produce per-screen and per-score GI outputs and hit-calls.

## Layout

- `config/` - configuration files (primary: `config.yaml`).
- `workflow/Snakefile` - main Snakemake workflow entrypoint.
- `workflow/rules/` - Snakemake rule files.
- `workflow/scripts/` - R scripts executed by rules. These expect to be run via Snakemake and access inputs/outputs/params through the `snakemake` object.
- `workflow/envs/` - optional conda/environment YAMLs used by rules (if present).
- `outputs/` - pipeline outputs (organized by screen/score) and logs.

## Requirements

- Snakemake (recommended installed in a conda environment).
- R (with tidyverse packages such as `readr`, `dplyr`, `tidyr`) available to the R scripts.
- Conda (recommended) or other environment manager to create reproducible environments. Environment YAMLs (if included) are under `workflow/envs/`.

Note: many R scripts use the `snakemake` object provided by Snakemake and redirect R stdout/messages to the rule `log` file when present. Static linters (lintr/R CMD check) may warn about `snakemake` being an undefined global; this is expected for scripts designed to run under Snakemake.

## Quick start (PowerShell)

1. Prepare your environment. Example using conda (adjust to the specific env YAML in `workflow/envs/` if present):

```powershell
# create and activate an environment (example file name)
# conda env create -f workflow/envs/smk-env.yaml -n smk-pipeline
# conda activate smk-pipeline

# verify snakemake is available
snakemake --version
```

2. Inspect or edit `config/config.yaml` to set screens, thresholds and other parameters used by rules.

3. Dry-run to see which jobs would be executed:

```powershell
snakemake -n --snakefile workflow/Snakefile --configfile config/config.yaml
```

4. Run the pipeline (example using 4 cores):

```powershell
snakemake --cores 4 --snakefile workflow/Snakefile --configfile config/config.yaml
```

5. Run a single target (for debugging). Example: run the GI score file for a particular screen/score:

```powershell
snakemake outputs/gi_scores/<screen>/all_scores/all_gis_<score>.tsv --cores 1 --snakefile workflow/Snakefile --configfile config/config.yaml
```

Replace `<screen>` and `<score>` with concrete names present in your `config/config.yaml`.

## Logs and troubleshooting

- Rule-level logs are written when a rule declares a `log:` path. By convention this pipeline writes logs to `outputs/logs/` (check the specific rule for the exact file). R scripts will redirect stdout and messages to the Snakemake-provided `log` file where available.
- If R scripts fail with an error referencing `snakemake`, ensure you run them through Snakemake (they expect the `snakemake` object). To test an R script interactively, you can construct a small wrapper that sets a `snakemake` list with `input`, `output`, `params` and `log` entries.

## Configuration

Edit `config/config.yaml` to control which screens and scores are processed and to set thresholds used by hit-calling rules (for example `HIT_THRESHOLD` and `DIFFERENTIAL_HIT_THRESHOLD`). Many rules use lambda params which select config entries based on rule wildcards (for example using the `{screen}` wildcard to pick screen-specific parameters).

## Common commands

PowerShell examples:

```powershell
# dry-run entire workflow
snakemake -n --snakefile workflow/Snakefile --configfile config/config.yaml

# run specific rule target (one file)
snakemake outputs/gi_scores/SCREEN2022/all_scores/all_gis_gi_score.tsv --cores 1 --snakefile workflow/Snakefile --configfile config/config.yaml

# run the full pipeline using available cores
snakemake --cores 0 --snakefile workflow/Snakefile --configfile config/config.yaml
```

Note: `--cores 0` lets Snakemake decide or use all available cores.