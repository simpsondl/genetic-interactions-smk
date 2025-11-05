# Genetic Interaction Analysis Pipeline

A Snakemake-based pipeline for end-to-end processing of genetic interaction (GI) screens, from raw counts through quality filtering, phenotype calculation, GI score computation, differential analysis, and hit identification.

**📊 Interactive Data Portal**: https://parpi.princeton.edu/map

**📖 Full Documentation**: See [PIPELINE_DOCUMENTATION.md](PIPELINE_DOCUMENTATION.md) for comprehensive setup guide, configuration options, and pipeline details.

## Overview

This pipeline:
- Processes multiple screens in parallel with flexible, per-screen configuration
- Calculates customizable phenotypes (Gamma, Tau, Rho, etc.)
- Computes genetic interaction scores at construct and gene levels
- Performs differential interaction analysis across conditions
- Applies sophisticated quality control and correlation filtering
- Supports compressed (ZIP) input files and high-precision numerical calculations

## Quick Start

### 1. Clone the repository
```bash
git clone https://github.com/simpsondl/genetic-interactions-smk.git
cd genetic-interactions-smk
```

### 2. Set up environment
```bash
# Create and activate conda environment
conda env create -f workflow/envs/smk-env.yaml
conda activate differential_gi_smk
```

### 3. Configure the pipeline
The provided config file (`config/config.yaml`) contains the parameters that were used to process the provided screens. By creating a new config file and passing its path as a command line argument, or by editing `config/config.yaml` directly, users can specify:
- Screen names and input file locations
- Phenotype definitions and treatment names
- Population doublings (if normalizing)
- Filter thresholds and hit-calling parameters

See [Configuration Guide](PIPELINE_DOCUMENTATION.md#configuration-guide) for detailed instructions.

### 4. Run the pipeline
```bash
# Dry run to preview execution plan
snakemake -n --snakefile workflow/Snakefile --configfile config/config.yaml

# Run complete pipeline (adjust cores as needed)
snakemake --cores 4 --snakefile workflow/Snakefile --configfile config/config.yaml
```

## Requirements

- **Snakemake** (≥9.12.0)
- **R** with tidyverse packages (readr, dplyr, broom, data.table)
- **Conda** (recommended for reproducible environments)

Environment specifications are provided in `workflow/envs/`.

## Key Features

- ✅ **Automatic metadata generation** for simplified input requirements
- ✅ **ZIP file support** for compressed count files (>100 MB)
- ✅ **Per-screen parameter overrides** for maximum flexibility
- ✅ **Multiple phenotype processing** in a single run
- ✅ **High-precision numerical I/O** to prevent rounding errors
- ✅ **Sequential correlation filtering** for quality control
- ✅ **Parallel execution** across screens, replicates, and phenotypes
- ✅ **Comprehensive logging** for debugging and transparency

## Output

Results are organized in `outputs/` directory:
- `counts/` - Validated counts with metadata
- `phenotypes/` - Calculated and filtered phenotypes
- `gi_scores/` - Genetic interaction scores (construct, gene, discriminant levels)
- `annotations/` - ID mapping files
- `misc_results/` - QC metrics and filter flags

The complete pipeline takes ~1 hour to process the two provided screens on a desktop (16GB RAM, 6 cores) and generates ~11 GB of final outputs.

## Documentation

For comprehensive documentation including:
- Detailed configuration guide with examples
- Pipeline stage descriptions and workflows
- Feature explanations and adaptability options
- Troubleshooting guide and validation checklist

See **[PIPELINE_DOCUMENTATION.md](PIPELINE_DOCUMENTATION.md)**

## Citation

If you use these datasets or workflow in your work, please cite the following reference:

> Simpson D, Ling J, Jing Y, Adamson B. Mapping the Genetic Interaction Network of PARP inhibitor Response. bioRxiv [Preprint]. 2023 Aug 20:2023.08.19.553986. doi: 10.1101/2023.08.19.553986. PMID: 37645833; PMCID: PMC10462155.