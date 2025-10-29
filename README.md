# Calling Differential Genetic Interactions

This repository contains a Snakemake-based pipeline for processing genetic interaction (GI) screens. It coordinates a set of R scripts (under `workflow/scripts/`) and Snakemake rules (under `workflow/rules/`) to produce per-screen, per-replicate, and per-phenotype genetic interaction scores and hit calls. Additionally calculates differential interaction scores across two screens, determines hits, and reproduces key analyses from our accompanying manuscript.

This pipeline can be adapted to work on other datasets. However, in order to use the workflow as is, input data must closely resemble the structure of the provided counts tables. In particular, column names and order are expected to match those in `workflow/data/counts`. If this is not the case, some steps may fail or present with unexpected behavior. 

The provided configuration file in `config/` will reproduce the analyses as presented in our manuscript; by adjusting these parameters, researchers can easily re-analyze our data as desired (e.g., with different filtering thresholds, or thresholds for determining hits).  

## Layout

- `config/` - configuration files (primary: `config.yaml`).
- `workflow/data` - counts files and metadata from Screen 2022 and Screen 2023
- `workflow/Snakefile` - main Snakemake workflow entrypoint.
- `workflow/rules/` - Snakemake rule files.
- `workflow/scripts/` - R scripts executed by rules. These expect to be run via Snakemake and access inputs/outputs/params through the `snakemake` object.
- `workflow/envs/` - optional conda/environment YAMLs used by rules (if present).
- `outputs/` - pipeline outputs (organized by screen/score) and logs, created when pipeline is run.

## Requirements

- Snakemake (recommended installed in a conda environment).
- R (with tidyverse packages such as `readr`, `dplyr`, `tidyr`) available to the R scripts.
- Conda (recommended) or other environment manager to create reproducible environments. Environment YAML is under `workflow/envs/`.

Note: many R scripts use the `snakemake` object provided by Snakemake and redirect R stdout/messages to the rule `log` file when present. Static linters (lintr/R CMD check) may warn about `snakemake` being an undefined global; this is expected for scripts designed to run under Snakemake.

## Quick start (Unix / Linux / macOS)

1. Get the code. If you don't yet have the repository locally, clone it. If you already have it, update to the latest changes:

```bash
# clone (first time)
git clone https://github.com/simpsondl/genetic-interactions-smk.git
cd genetic-interactions-smk

# or update an existing local clone
#git pull origin main
```

2. Prepare your environment. Example using conda:

```bash
# create and activate an environment (example file name)
# conda env create -f workflow/envs/smk-env.yaml
# conda activate differential_gi_smk

# verify snakemake is available
snakemake --version
```

3. Inspect or edit `config/config.yaml` to set screens, thresholds and other parameters used by rules.

4. Dry-run to see which jobs would be executed:

```bash
snakemake -n --snakefile workflow/Snakefile --configfile config/config.yaml
```

5. Run the complete pipeline (example using 4 cores):

```bash
snakemake --cores 4 --snakefile workflow/Snakefile --configfile config/config.yaml
```

6. Run a single target (for debugging). Example: run the GI score file for a particular screen/score:

```bash
snakemake outputs/gi_scores/<screen>/all_scores/all_gis_<score>.tsv --cores 1 --snakefile workflow/Snakefile --configfile config/config.yaml
```

Replace `<screen>` and `<score>` with concrete names present in your `config/config.yaml`.

## Expected outputs

If you run the complete pipeline using the parameters in the provided config file, then the pipeline will create a new `outputs` directory and place all associated output files, including logs, there. The complete pipeline takes approximately 4 hours to run on a desktop computer with 16GB RAM and 6 cores, and generates ~12 Gb of data. Ensure that enough space is available before starting the pipeline.

## Logs and troubleshooting

- If you have trouble establishing the conda environment, please ensure that your institution has access to the requested channels and change them as needed.
- If you have trouble establishing the conda environment using the full provided YAML, an alternate YAML is available (`workflow/envs/smk-env-no-r.yaml`) which does not specify installation of any R packages, just a base R version. Once that environment has been made, you can activate it, open R, and install depedencies:

```bash
conda env create -f workflow/envs/smk-env-no-r.yaml
conda activate differential_gi_smk
R
```

```R
install.packages(c("readr", "dplyr", "data.table", "broom", "ggplot2", "ggrepel"))
q() #no need to save workspace
```

- Rule-level logs are written when a rule declares a `log:` path. By convention this pipeline writes logs to `outputs/logs/` (check the specific rule for the exact file). R scripts will redirect stdout and messages to the Snakemake-provided `log` file where available.
- If R scripts fail with an error referencing `snakemake`, ensure you run them through Snakemake (they expect the `snakemake` object). To test an R script interactively, you can construct a small wrapper that sets a `snakemake` list with `input`, `output`, `params` and `log` entries.

## Configuration

Edit `config/config.yaml` to control which screens and scores are processed and to set thresholds used by hit-calling rules (for example `HIT_THRESHOLD` and `DIFFERENTIAL_HIT_THRESHOLD`). Many rules use lambda params which select config entries based on rule wildcards (for example using the `{screen}` wildcard to pick screen-specific parameters).

## Common commands

Shell examples (bash), run from within the repository directory `genetic-interactions-smk`:

```bash
# dry-run entire workflow
snakemake -n --snakefile workflow/Snakefile --configfile config/config.yaml

# run specific rule target (one file)
snakemake outputs/gi_scores/SCREEN2022/all_scores/all_gis_gi_score.tsv --cores 1 --snakefile workflow/Snakefile --configfile config/config.yaml

# run the full pipeline using available cores (adjust --cores as available)
snakemake --cores 8 --snakefile workflow/Snakefile --configfile config/config.yaml
```

## Support for zipped counts files

This pipeline accepts compressed counts files in ZIP or TSV format for the initial counts input. This was done to aid distribution of the raw counts files, which are large (>100 Mb). The R scripts that read counts (`workflow/scripts/apply_filters.R` and `workflow/scripts/calculate_phenotypes.R`) have been made zip-aware and will accept either:

- A plain TSV: `data/counts/{screen}_raw_counts.tsv`, or
- A ZIP file: `data/counts/{screen}_raw_counts.zip` which contains one or more files; the scripts will look for the first entry matching `_raw_counts.tsv` (case-insensitive) and read it.

## Notes and troubleshooting

- The R scripts search inside the ZIP for the first file whose name ends with `_raw_counts.tsv` (case-insensitive). If your ZIP uses a different inner filename, either rename the inner file to match that pattern or create a ZIP that contains the expected filename.
- The R scripts still expect to be run via Snakemake (they rely on the `snakemake` object for inputs/outputs/params). To test the scripts directly, construct a small wrapper that sets a `snakemake` list with the required keys (or run via Snakemake with a small test config).

## Citation

If you use these datasets or workflow in your work, please cite the following reference: