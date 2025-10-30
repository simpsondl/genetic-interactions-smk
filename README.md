# Calling Differential Genetic Interactions

This repository contains a Snakemake-based pipeline for processing genetic interaction (GI) screens. It coordinates a set of R scripts (under `workflow/scripts/`) and Snakemake rules (under `workflow/rules/`) to produce per-screen, per-replicate, and per-phenotype genetic interaction scores and hit calls. Additionally calculates differential interaction scores across two screens, determines hits, and reproduces key analyses from our accompanying manuscript.

This pipeline can be adapted to work on other datasets. However, in order to use the workflow as is, input data must closely resemble the structure of the provided counts tables. In particular, column names and order are expected to match those in `workflow/data/counts`. If this is not the case, some steps may fail or present with unexpected behavior. 

The provided configuration file in `config/` will reproduce the analyses as presented in our manuscript; by adjusting these parameters, researchers can easily re-analyze our data as desired (e.g., with different filtering thresholds, or thresholds for determining hits).  

Data presented in our manuscript, along with data from two other manuscripts, are also available through our interactive data portal: https://parpi.princeton.edu/map

## Requirements

- Snakemake (recommended installed in a conda environment).
- R (with tidyverse packages such as `readr`, `dplyr`, `tidyr`) available to the R scripts.
- Conda (recommended) or other environment manager to create reproducible environments. Environment YAML for establishing the same conda environment used in our manuscript is under `workflow/envs/`.


## Layout

- `config/` - configuration file (primary: `config.yaml`).
- `workflow/data` - counts files and metadata from Screen 2022 and Screen 2023
- `workflow/Snakefile` - main Snakemake workflow entrypoint.
- `workflow/rules/` - Snakemake rule files.
- `workflow/scripts/` - R scripts executed by rules. These expect to be run via Snakemake and access inputs/outputs/params through the `snakemake` object.
- `workflow/envs/` - conda/environment YAMLs to produce same environment used in original data analysis.
- `outputs/` - pipeline outputs (organized by screen/score) and logs, created when pipeline is run.

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


## Configuration

Edit `config/config.yaml` to control which screens and scores are processed and to set thresholds used by hit-calling rules (for example `HIT_THRESHOLD` and `DIFFERENTIAL_HIT_THRESHOLD`). Many rules use lambda params which select config entries based on rule wildcards (for example using the `{screen}` wildcard to pick screen-specific parameters).


## Expected outputs

If you run the complete pipeline using the parameters in the provided config file, then the pipeline will create a new `outputs` directory and place all associated output files, including logs, there. The complete pipeline takes approximately 3 hours to run on a desktop computer with 16GB RAM and 6 cores, and generates ~10 Gb of data. Ensure that enough space is available before starting the pipeline.

## Support for zipped counts files

This pipeline accepts compressed counts files in ZIP or TSV format for the initial counts input. This was done to aid distribution of the raw counts files, which are large (>100 Mb). The R scripts that read counts (`workflow/scripts/apply_filters.R` and `workflow/scripts/calculate_phenotypes.R`) have been made zip-aware and will accept either:

- A plain TSV: `data/counts/{screen}_raw_counts.tsv`, or
- A ZIP file: `data/counts/{screen}_raw_counts.zip` which contains one or more files; the scripts will look for the first entry matching `_raw_counts.tsv` (case-insensitive) and read it.


## Logs and troubleshooting

- Rule-level logs are written when a rule declares a `log:` path. This pipeline writes logs to `outputs/logs/` with a name matching the rule. R scripts will redirect stdout and messages to the Snakemake-provided `log` file where available. Logging files are descriptive and contain some debugging information as well as summaries.
- This pipeline supports dual logging, where messages are simultaneously sent to stdout as well as the indicated log file. This is helpful for monitoring workflow progress; however, if snakemake was invoked with multiple cores, messages from all cores print to stdout without indication of which job each message belongs to and messages can sometimes interrupt each other. Individual log files are deconvoluted and allow progress of parallel jobs to be monitored independently. 
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

- To test an R script interactively, you can construct a small wrapper that sets a `snakemake` list with `input`, `output`, `params` and `log` entries.


## Citation

If you use these datasets or workflow in your work, please cite the following reference:

> Simpson D, Ling J, Jing Y, Adamson B. Mapping the Genetic Interaction Network of PARP inhibitor Response. bioRxiv [Preprint]. 2023 Aug 20:2023.08.19.553986. doi: 10.1101/2023.08.19.553986. PMID: 37645833; PMCID: PMC10462155.