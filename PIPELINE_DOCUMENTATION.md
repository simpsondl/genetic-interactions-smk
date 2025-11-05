# Genetic Interaction Pipeline Documentation

## Table of Contents

1. [Overview](#overview)
2. [Pipeline Architecture](#pipeline-architecture)
3. [Configuration Guide](#configuration-guide)
   - [Essential Settings](#essential-settings)
   - [Screen-Specific Parameters](#screen-specific-parameters)
   - [Filter Parameters](#filter-parameters)
   - [Advanced Settings](#advanced-settings)
4. [Pipeline Stages](#pipeline-stages)
   - [Stage 1: Metadata Generation](#stage-1-metadata-generation)
   - [Stage 2: Preprocessing](#stage-2-preprocessing)
   - [Stage 3: GI Score Calculation](#stage-3-gi-score-calculation)
   - [Stage 4: Differential Analysis](#stage-4-differential-analysis)
   - [Stage 5: Hit Calling](#stage-5-hit-calling)
   - [Stage 6: Clustering & Visualization](#stage-6-clustering--visualization)
5. [Key Features & Adaptability](#key-features--adaptability)
6. [Output Structure](#output-structure)
7. [Troubleshooting](#troubleshooting)

---

## Overview

This Snakemake-based pipeline processes genetic interaction (GI) screening data end-to-end, from raw counts through quality filtering, phenotype calculation, genetic interaction score computation, differential analysis, and hit identification. The pipeline is designed to be **highly configurable** and **adaptable** to different experimental designs while maintaining reproducibility and precision.

**Key Capabilities:**
- Process multiple screens simultaneously
- Calculate multiple phenotype types (Gamma, Tau, Rho, etc.)
- Compute genetic interaction scores at construct and gene levels
- Perform differential interaction analysis across conditions
- Apply sophisticated quality control filters
- Support both replicate-level and averaged analyses
- Generate reproducible results with high numerical precision

---

## Pipeline Architecture

The pipeline follows a **modular design** with six major stages executed in sequence:

```
Raw Counts (ZIP/TSV)
    ↓
[1] Metadata Generation → Validates/generates construct annotations
    ↓
[2] Preprocessing → Filters, phenotypes, normalization, correlation QC
    ↓
[3] GI Score Calculation → Construct-level and gene-level scores
    ↓
[4] Differential Analysis → Compares conditions (e.g., Tau - Gamma = Nu)
    ↓
[5] Hit Calling → Identifies significant interactions using discriminant scores
    ↓
[6] Clustering & Visualization → Optional diagnostic plots and figures
```

**Key Design Principles:**
- **Temporary file cleanup**: Intermediate RDS workspaces are automatically cleaned up
- **Parallel processing**: Multiple screens and replicates can be processed simultaneously
- **Fail-fast validation**: Early validation catches configuration errors before long computations
- **High precision**: Custom I/O functions preserve numerical precision across stages

---

## Configuration Guide

The pipeline is controlled by `config/config.yaml` by default. You can either edit that file directly, or create your own config file, and pass that file as a command line option: 
```bash
snakemake -c6 --snakefile workflow/Snakefile --configfile path/to/your/config.yaml process_all
```

This section explains how to set up your configuration.

### Global parameters

These essential parameters define the scope of the pipeline run, including where input and output files are, and what files (screens) to process.

#### 1. Directory Configuration
```yaml
OUTPUTS_DIR: outputs  # Where all results are written
LOGS_DIR: logs        # Where log files are stored
```
💡 **Tip**: Use relative paths from the working directory (where you submitted the job) or absolute paths.

#### 2. Screen Selection
```yaml
SCREENS: [screen2022, screen2023]
```
- List all screens to process
- Names must match filename prefixes in your counts files
- Pipeline processes all listed screens in parallel (assuming multiple cores are used)

#### 3. Input Counts Location
```yaml
COUNTS_FILEPATH_PATTERN: manuscript_data/counts/{screen}_raw_counts.zip
```
- Provide a full path or one relative to the working directory. Use `{screen}` as a placeholder (snakemake automatically replaces this with the screen names defined in SCREENS), but otherwise the file path should match your input counts file name.
- Supports both `.tsv` and `.zip` formats
- **Cool Feature**: ZIP files are automatically extracted - the pipeline finds and uses the first file ending in `_counts.tsv`, `_counts_only.tsv`, or `_counts_with_metadata.tsv` for each screen.

Example patterns:
```yaml
# For TSV files directly:
COUNTS_FILEPATH_PATTERN: manuscript_data/counts/{screen}_raw_counts.tsv

# For ZIP files (recommended for large files):
COUNTS_FILEPATH_PATTERN: manuscript_data/counts/{screen}_raw_counts.zip

# For different directory structure:
COUNTS_FILEPATH_PATTERN: /path/to/data/{screen}/{screen}_counts.zip
```

### Phenotype definition parameters

These settings control the behavior of specific steps in the pipeline. Importantly, each of these settings can be overridden in a screen-specific manner, so that different settings can be applied to different screens in the same run (see `Enabling screen-specific settings`).

#### 4. Count Column Structure
```yaml
COUNT_PREFIXES: [T0, DMSO, TREATED]  # Treatment/timepoint identifiers
REPLICATES: [R1, R2]                  # Replicate suffixes
INITIAL_CONDITION_ID: T0              # Which prefix is the starting condition
```
These define your column naming scheme. For example:
- `T0_R1`, `T0_R2` (initial condition replicates)
- `DMSO_R1`, `DMSO_R2` (control treatment replicates)
- `TREATED_R1`, `TREATED_R2` (experimental treatment replicates)

The pipeline supports `_`, `-`, or `.` as separators between COUNT_PREFIXES and REPLICATES in column names, as well as no separator.

#### 5. Phenotype Definitions
```yaml
PHENOTYPE_PREFIX_MAP:
  GAMMA: [DMSO, T0]      # Control phenotype: DMSO / T0
  TAU: [TREATED, T0]     # Treatment phenotype: TREATED / T0
  RHO: [TREATED, DMSO]   # Relative phenotype: TREATED / DMSO
```
- **Format**: `Phenotype_Name: [Numerator, Denominator]`
- Prefixes must match those in `COUNT_PREFIXES`
- You can define additional phenotypes as long as both `Numerator` and `Denominator` appears in COUNT_PREFIXES 
- Names can be customized (not limited to Gamma/Tau/Rho)
- Any number of phenotypes can be defined (1+)

#### 6. GI Phenotype Selection
```yaml
GI_PHENOTYPES: [GAMMA, TAU]
```
- Specifies which phenotypes to use for genetic interaction scoring
- Must be a subset of keys in `PHENOTYPE_PREFIX_MAP`
- Reduces computation time if you don't need all genetic interaction scores calculated for all phenotypes

### Filter parameters

#### 7. Count Quality Filters
```yaml
INDIVIDUAL_SGRNA_MEDIAN_THRESHOLD: 35
COMBINATION_SGRNA_COUNT_THRESHOLD: 50
```
- `INDIVIDUAL_SGRNA_MEDIAN_THRESHOLD`: Minimum median count across replicates for single guides
- `COMBINATION_SGRNA_COUNT_THRESHOLD`: Minimum count threshold for guide combinations
- Constructs failing these filters are flagged and removed from genetic interaction score calculations

#### 8. Phenotype Processing
```yaml
PSEUDOCOUNT: 10              # Added before log transformation
NORMALIZE: TRUE              # Normalize by population doublings
AVERAGE_REPLICATES: TRUE     # Generate averaged phenotype columns
```
- NORMALIZE uses population doubling information to adjust phenotype calculations, enabling better comparison of phenotypes across screen arms (information must be collected during experiment)
- AVERAGE_REPLICATES will average replicates together and the pipeline will additionally process these averaged phenotypes

#### 9. Correlation Filtering
```yaml
NO_CORRELATION_FILTER:
  - GAMMA
  - TAU

NO_CORRELATION_THRESHOLD: 0.25
CORRELATION_FILTER_MODE: avg_only  # Options: avg_only, any, all
```
**🎯 Smart Feature**: Sequential correlation filtering
- Applied in the order specified in `NO_CORRELATION_FILTER`
- Checks correlation between dual guide combinations and their constituent single guides for each sgRNA
- Three modes:
  - `avg_only`: Filter based on average correlation (used in manuscript)
  - `any`: Filter if ANY single replicate (or average, if present) shows poor correlation
  - `all`: Filter only if ALL replicates (and average, if present) show poor correlation

#### 10. Differential Comparisons
```yaml
DIFFERENTIAL_COMPARISONS:
  NU: [TAU, GAMMA]  # NU = TAU - GAMMA (differential interaction)
```
- Define which phenotypes to compare
- **Format**: `Name: [Phenotype1, Phenotype2]` calculates Phenotype1 - Phenotype2
- Can define multiple comparisons

**Example - Multiple comparisons:**
```yaml
DIFFERENTIAL_COMPARISONS:
  NU: [TAU, GAMMA]
  DELTA: [RHO, GAMMA]
  EPSILON: [TAU, RHO]
```

#### 11. Hit Calling Thresholds
```yaml
HIT_THRESHOLDS:
  GAMMA: 0.995    # 99.5th percentile for synergy, 0.5th for suppression
  TAU: 0.995
  RHO: 0.995
  NU: 1           # 100th percentile (most extreme) for differential
```
- Thresholds represent percentile cutoffs of discriminant scores based on the distribution of discriminant scores from pairs involving at least one non-targeting pseudogene (i.e., controls)
- Values between 0 and 1
- Higher values = more stringent (fewer hits)

#### 12. Clustering Parameters
```yaml
PHENOTYPES_TO_CLUSTER: [GAMMA.OI.Avg, TAU.OI.Avg, NU.OI.Avg]
```
- Generates diagnostic clustering plots
- Use orientation-independent averaged scores for best results


### Screen-Specific Parameters

**🌟 Adaptability Feature**: Override global settings for individual screens by prefixing parameters with the screen name. Use all capital letters in config file variable names.

```yaml
# Global default:
COUNT_PREFIXES: [T0, DMSO, TREATED]

# Screen-specific override:
SCREEN2023_COUNT_PREFIXES: [T0, CONTROL, EXPERIMENTAL]
SCREEN2023_PHENOTYPE_PREFIX_MAP:
  GAMMA: [CONTROL, T0]
  TAU: [EXPERIMENTAL, T0]
```

**Example use cases:**
- Different treatment names across screens
- Different replicate structures
- Screen-specific thresholds
- Varied phenotype definitions

The following configuration variables can be set in a screen-specific way:

- `COUNT_PREFIXES`
- `REPLICATES`
- `INITIAL_CONDITION_ID`
- `PHENOTYPE_PREFIX_MAP`
- `GI_PHENOTYPES`
- `INDIVIDUAL_SGRNA_MEDIAN_THRESHOLD`
- `COMBINATION_SGRNA_COUNT_THRESHOLD`
- `PSEUDOCOUNT`
- `NORMALIZE`
- `AVERAGE_REPLICATES`
- `NO_CORRELATION_FILTER`
- `NO_CORRELATION_THRESHOLD`
- `CORRELATION_FILTER_MODE`
- `HIT_THRESHOLDS`

The _DOUBLINGS variable(s) *must* be set in a screen-specific way:

```yaml
SCREEN2022_DOUBLINGS:
  DMSO:
    R1: 7.19
    R2: 7.73
  NIRAP:
    R1: 5.32
    R2: 5.36

SCREEN2023_DOUBLINGS:
  DMSO:
    R1: 8.04
    R2: 7.98
  NIRAP:
    R1: 5.40
    R2: 4.98
```
- Only required with `NORMALIZE: TRUE` 
- Define per-treatment, per-replicate values
- Treatment names must match your `COUNT_PREFIXES`



#### 13. Internal Metadata Columns (⚠️ DANGER ZONE)
```yaml
EXPECTED_META_COLUMNS:
  - FirstPosition
  - FirstPseudogene
  - SecondPosition
  - SecondPseudogene
  - ConstructID
  - GuideCombinationID
  - PseudogeneCombinationID
  - PseudogeneCombinationName
  - Category
  - Control
  - Identical
  - Orientation
```
**⚠️ Warning**: Only modify if you intend to re-work downstream scripts. These column names are expected and are called by name throughout the pipeline. They are auto-generated if missing from the provided counts file.

---

## Pipeline Stages

### Stage 1: Metadata Generation

**Purpose**: Validate, enrich, and standardize input counts files so downstream rules receive a consistent, annotated counts table. This stage accepts both fully annotated inputs and raw counts, and will generate the minimal metadata required for downstream GI scoring when needed.

**Rules:**
1. `prepare_counts_with_metadata` — read input (ZIP or TSV), detect/repair column names, and produce a canonical counts-with-metadata table.
2. `generate_id_maps` — create mapping tables that translate internal construct/guide identifiers to canonical IDs used elsewhere in the pipeline.
3. `validate_counts` — run structural and content checks and emit a validation marker used by downstream rules.

**What happens:**
- Input detection: `prepare_counts_with_metadata` accepts either a plain TSV or a ZIP file (the ZIP is inspected for the first matching `_counts.tsv`, `_counts_only.tsv`, or `_counts_with_metadata.tsv` entry).
- Column harmonization: column names are checked against the expected naming conventions (prefix + replicate). The code tolerates common separators (`.`, `_`, `-`) and normalizes column names so downstream scripts can rely on predictable names.
- Metadata verification and partial generation: if required metadata columns are missing, the script attempts to infer or generate them deterministically (see "Partial metadata generation" below). If inference is impossible (for example, missing all guide identifiers), the pipeline exits with a clear error message.
- ID maps creation: `generate_id_maps` writes several per-screen mapping files (construct, guide combination, gene combination, pseudogene mapping) which are used to annotate and aggregate results downstream.
- Validation: `validate_counts` performs sanity checks (types, duplicates, control presence) and writes `outputs/misc_results/{screen}_counts_validated.txt`. Downstream rules depend on this marker so validation failures stop the workflow early.

**Partial metadata generation (use cases & behavior):**
- Many distributed count files omit derived columns such as `ConstructID`, `GuideCombinationID`, or aggregated gene-pair IDs. `prepare_counts_with_metadata` will generate these when possible using heuristics based on available guide/position columns and consistent naming.
- Example generated fields: `ConstructID`, `GuideCombinationID`, `PseudogeneCombinationID`, `PseudogeneCombinationName`, `Category` (e.g., NT+NT, control), and orientation markers.
- The script preserves any valid supplied metadata, generates missing columns deterministically (so runs are reproducible), and logs all assumptions. When a required provenance field is entirely missing and cannot be inferred, the script exits with an actionable error describing what the user must supply.

**ID maps produced (what they contain & why):**
- `outputs/annotations/{screen}_construct_id_map.tsv` — maps raw construct identifiers (barcode or raw string) to canonical `ConstructID` used throughout the pipeline. Helps trace results back to raw inputs.
- `outputs/annotations/{screen}_guide_combination_id_map.tsv` — maps dual-guide constructs to combined guide IDs (e.g., `sgc_...`) and records the two constituent guide IDs and their positions; used for single-guide extraction and correlation checks.
- `outputs/annotations/{screen}_gene_combination_id_map.tsv` — aggregates constructs by targeted gene pair and provides gene-level identifiers used for gene-level scoring and reporting.
- `outputs/annotations/{screen}_pseudogene_id_map.tsv` — provides mappings for pseudogene or legacy identifiers when datasets include two-level mappings (construct → pseudogene → gene), preserving compatibility with older naming schemes.

**Data integrity checks and safeguards:**
- Column existence and type checks: required count columns and metadata columns are verified; non-numeric count columns raise errors.
- Replicate matching: ensures numerator/denominator pairs exist for configured `COUNT_PREFIXES` and `REPLICATES`; missing pairs are reported as errors.
- Non-targeting control checks: verifies presence of NT controls used for centering and QC; if none are present the pipeline warns or aborts depending on downstream needs.
- Duplicate detection: duplicate `ConstructID` or conflicting mappings are detected and reported.
- Validation marker: `validate_counts` writes `outputs/misc_results/{screen}_counts_validated.txt` (a temporary marker). Downstream rules require this marker, so failures halt processing early.

**Inputs:**
- Raw counts files (ZIP or TSV) or an already-annotated counts TSV.

**Outputs:**
- `outputs/counts/{screen}_counts_with_metadata.tsv` — canonical counts table (original or augmented with generated metadata).
- `outputs/annotations/{screen}_construct_id_map.tsv`
- `outputs/annotations/{screen}_guide_combination_id_map.tsv`
- `outputs/annotations/{screen}_gene_combination_id_map.tsv`
- `outputs/annotations/{screen}_pseudogene_id_map.tsv`
- `outputs/misc_results/{screen}_counts_validated.txt` (validation marker)

**Logging & transparency:**
- All generated metadata, assumptions, and non-fatal warnings are recorded in the rule logs, which can be found in `logs/{screen}/metadata_generation/` by default.

For the canonical internal column names the pipeline expects, see the "Internal Metadata Columns (⚠️ DANGER ZONE)" section above — these names are referenced directly across the R scripts and should only be changed with corresponding code updates.

---

### Stage 2: Preprocessing

**Purpose**: Filter low-quality constructs, compute per-replicate phenotypes, center to non-targeting controls, optionally normalize by population doublings, and prepare orientation-independent and single-sgRNA phenotype tables for GI scoring.

**Rules covered:**
1. `apply_count_filters` — flag low-count constructs and create filter flags used downstream
2. `calculate_phenotypes` — compute per-replicate phenotype columns (log2 ratios) and center to NT medians
3. `normalize_phenotypes` — divide phenotypes by population-doubling values (if enabled)
4. `create_averaged_phenotypes` — compute orientation-independent and averaged replicate phenotype tables
5. `apply_correlation_filter` — run no-correlation filters producing filtered phenotype tables and QC reports

What each rule does (details):
- `apply_count_filters` reads the canonical counts-with-metadata, applies `INDIVIDUAL_SGRNA_MEDIAN_THRESHOLD` and `COMBINATION_SGRNA_COUNT_THRESHOLD` (screen-specific overrides allowed), and writes `outputs/misc_results/{screen}_filter_flags.tsv`. Filter flags are preserved and passed to phenotype calculation so filtered constructs can be excluded from scoring but remain documented.
- `calculate_phenotypes` builds fractional abundances (adds pseudocount), computes per-phenotype per-replicate log2 ratios using `PHENOTYPE_PREFIX_MAP`, and then centers each phenotype column by subtracting the median across non-targeting (`Category == 'NT+NT'`) constructs. The centering step is deterministic and logged. This rule creates `outputs/phenotypes/{screen}_phenotypes.raw.tsv` (temporary).
- `normalize_phenotypes` (if `NORMALIZE: TRUE`) divides the already-centered phenotype values by the appropriate `DOUBLINGS` values (screen-specific mapping required for normalization). For Rho-type phenotypes the code uses (d_treated - d_reference) as the divisor. This rule writes a normalized phenotype TSV (temporary).
- `create_averaged_phenotypes` produces three main outputs per screen: processed phenotypes, orientation-independent phenotypes (per GuideCombinationID), and single-sgRNA phenotypes derived from NT-containing constructs. It uses `AVERAGE_REPLICATES` and `PHENOTYPE_PREFIX_MAP` to determine which averaged columns to make.
- `apply_correlation_filter` evaluates correlation between dual-guide constructs and their constituent single-guide phenotypes for phenotypes listed in `NO_CORRELATION_FILTER`. Modes (`avg_only`, `any`, `all`) control stringency. It outputs filtered phenotype tables per phenotype that fail the correlation checks, a full filter flags table, and correlation results/summary TSVs.

Key inputs, outputs and artifacts:
- Inputs: `outputs/counts/{screen}_counts_with_metadata.tsv`, `outputs/misc_results/{screen}_filter_flags.tsv` (from apply_count_filters)
- Temp outputs: `outputs/phenotypes/{screen}_phenotypes.raw.tsv`, `outputs/phenotypes/{screen}_phenotypes.normalized.tsv`
- Persistent outputs: `outputs/phenotypes/{screen}_processed_phenotypes.tsv`, `outputs/phenotypes/{screen}_orientation_independent_phenotypes.tsv`, `outputs/phenotypes/{screen}_single_sgRNA_phenotypes.tsv`, `outputs/phenotypes/{screen}_filtered_{phenotype}_phenotypes.tsv`, and `outputs/misc_results/{screen}_*_correlation_*.tsv`

Checks & safeguards:
- Ensures count columns required by `PHENOTYPE_PREFIX_MAP` are present and numeric
- Verifies NT controls exist before centering; if absent, logs and may abort depending on downstream needs
- Confirms `DOUBLINGS` mapping exists for each used prefix/suffix when normalization is enabled and raises a clear error if not
- All filter decisions are logged and stored in flag files for reproducibility

Performance and practical notes:
- Normalization is optional; if you lack doubling measurements, set `NORMALIZE: FALSE` and proceed without doubling-based scaling
- Use `snakemake -n` to preview which phenotype files will be produced for your configured replicates and averages

---

### Stage 3: GI Score Calculation

**Purpose**: Fit models that compare observed dual-guide phenotypes to expectation from single guides, derive construct-level GI scores, aggregate to gene level, and compute discriminants used for hit calling.

**Rules covered:**
1. `compute_genetic_interaction_scores` — model fitting and per-construct GI scores
2. `calculate_gene_level_scores` — aggregate construct scores to gene-pair level
3. `calculate_discriminant_scores` — normalize scores relative to control distributions for hit-calling

Detailed flow and behavior:
- Inputs: orientation-independent and single-sgRNA phenotype tables produced in Stage 2. The `compute_genetic_interaction_scores` rule resolves score wildcards (e.g., `GAMMA.OI.R1`, `GAMMA.OI.Avg`) based on `GI_PHENOTYPES`, `REPLICATES`, and `AVERAGE_REPLICATES`.
- Modeling: for each sgRNA the pipeline fits a quadratic model between combinations involving that sgRNA and the other genes' individual phenotypes. The GI score is the residual (observed − expected). Model estimates and statistics are written to `outputs/gi_scores/{screen}/models/` for diagnostic use.
- Workspaces & precision: intermediate R workspaces are saved with high-precision IO (`r_precise_io.R`) as `.rds` files. These are marked as `temp()` so Snakemake cleans them up after downstream rules finish.
- Aggregation: `calculate_gene_level_scores` aggregates construct-level GI scores into gene-level scores by taking the mean.
- Discriminant: `calculate_discriminant_scores` computes a discriminant (`-log10(Mann-Whitney p-value) * |Interaction Score|`) used for percentile-based hit calling.

Outputs and artifacts:
- Construct-level: `outputs/gi_scores/{screen}/construct_scores/all_gis_{score}.tsv` and per-construct directories under `individual_scores/{score}`
- Gene-level: `outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv`
- Discriminant and hits-ready tables: `outputs/gi_scores/{screen}/discriminant_scores/discriminant_scores_{score}.tsv`
- Model diagnostics: `outputs/gi_scores/{screen}/models/model_estimates_{score}.tsv`, `model_stats_{score}.tsv`

Checks & safeguards:
- Ensures sufficient non-control constructs exist to estimate control distributions; if not, discriminant estimation will fail with clear logs
- Detects and reports constructs with insufficient single-guide coverage
- Logs model convergence and unusual parameter estimates in model stat files

Practical notes:
- The rule supports per-replicate and averaged (`.Avg`) score generation; use `GI_PHENOTYPES` and `REPLICATES` in `config.yaml` to control which scores are built

---

### Stage 4: Differential Analysis

**Purpose**: Produce differential GI scores that capture how interactions change between two phenotypes (for example, treatment vs control).

**Rule covered:**
- `calculate_differential_scores`

Detailed behavior:
- The rule reads construct-level GI scores for the two phenotypes specified in `DIFFERENTIAL_COMPARISONS` (mapping keys → `[treated, reference]`) and computes differential scores as (treated − reference) per construct.
- Supports the same replicate suffixes and averaged `.Avg` outputs as construct GI scoring; the helper functions build replicate name lists dynamically from `REPLICATES` and `AVERAGE_REPLICATES`.
- Differential gene-level aggregation and discriminant calculation follow the same approach as Stage 3, producing gene-level differential tables and discriminant-normalized differential scores.

Inputs & outputs:
- Inputs: two construct-level GI score files and associated `.rds` workspaces for the phenotypes being compared
- Outputs: `outputs/gi_scores/{screen}/construct_scores/all_gis_{comparison}.{rep}.tsv`, `outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{comparison}.{rep}.tsv`, and discriminant files for the differential comparison

Checks & safeguards:
- Validates that both phenotypes in a comparison exist in the configured `GI_PHENOTYPES`
- Reports missing replicate-specific construct files if a replicate is absent for either phenotype

Practical notes:
- Define any number of differential comparisons in `DIFFERENTIAL_COMPARISONS` in `config.yaml` to create custom contrasts (Nu, Delta, etc.)

---

### Stage 5: Hit Calling

**Purpose**: Convert discriminant scores into categorical hit calls (Interaction (synergistic, buffering) / No Interaction) using percentile-based thresholds, supporting screen- and phenotype-specific thresholds.

**Rule covered:**
- `call_hits`

How thresholds are resolved:
- The rule first attempts to find a screen-specific mapping in `_HIT_THRESHOLDS` (a mapping of phenotype → percentile). If not present it falls back to global `HIT_THRESHOLDS`.
- Percentiles are applied to the discriminant score distribution derived from control constructs (pairs containing non-targeting pseudogenes) to compute cutoffs.

Classification & outputs:
- Output file: `outputs/gi_scores/{screen}/discriminant_scores/discriminant_hits_{score}.tsv` which adds a `hit` column.

Checks & safeguards:
- Validates that the control-derived distribution exists and contains enough values to estimate percentile cutoffs. If insufficient controls are present, the rule fails with a descriptive error.
- Logs the thresholds used and the number of hits called in the rule log for auditing.

Practical notes:
- To increase/decrease stringency adjust `HIT_THRESHOLDS` in `config.yaml`. For experimental tuning, run `call_hits` on a specific discriminant file to quickly iterate on thresholds.

---

### Stage 6: Clustering & Visualization

**Purpose**: Generate diagnostic plots and manuscript figures (optional).

**Rules:**
- `diagnostic_plot` - Generate clustering diagnostic plots
- `generate_manuscript_figures` - Create publication-quality figures

**Diagnostic Plots:**
- Visualize score distributions
- Identify outliers
- QC clustering behavior
- Can be generated for selected phenotypes via `PHENOTYPES_TO_CLUSTER`

**Manuscript Figures:**
- Reproduces specific figures from publication
- Includes correlation plots, score distributions, example interactions
- Optional target (not run by default)

**Outputs:**
- `outputs/gi_scores/{screen}/clusters/diagnostic_plot_{score}.svg`
- `outputs/manuscript_figures/figure_*.{svg,png}` (if requested)

---

## Key Features & Adaptability

### 1. 🗜️ **Compressed File Support**
- **Feature**: Automatic ZIP file detection and extraction
- **Why**: Large count files (>100 MB) are easier to distribute compressed
- **How**: Pipeline searches for `*_counts.tsv` inside ZIP files
- **Benefit**: No manual extraction needed

### 2. 🔁 **Flexible Screen Configuration**
- **Feature**: Process any number of screens simultaneously
- **Why**: Batch analysis of multiple experiments
- **How**: Simply list screens in `SCREENS` parameter
- **Benefit**: Consistent processing, parallel execution

### 3. 🎛️ **Per-Screen Parameter Overrides**
- **Feature**: Override any global parameter for specific screens
- **Why**: Screens often have different experimental conditions
- **How**: Prefix parameter name with screen name (e.g., `SCREEN2023_DOUBLINGS`)
- **Examples**:
  - Different treatment names
  - Different replicate numbers
  - Screen-specific thresholds
  - Unique phenotype definitions
- **Benefit**: Maximum flexibility without separate config files

### 4. 🎨 **Custom Phenotype Definitions**
- **Feature**: Define unlimited phenotypes with any names
- **Why**: Different experimental designs require different comparisons
- **How**: Add entries to `PHENOTYPE_PREFIX_MAP`
- **Examples**:
  ```yaml
  PHENOTYPE_PREFIX_MAP:
    EarlyGrowth: [T4, T0]
    LateGrowth: [T8, T4]
    DrugEffect: [TREATED, DMSO]
    DoseResponse: [HIGH_DOSE, LOW_DOSE]
  ```
- **Benefit**: Not limited to standard Gamma/Tau/Rho nomenclature

### 5. 🔬 **Multi-Phenotype Processing**
- **Feature**: Calculate GI scores for multiple phenotypes in one run
- **Why**: Compare interactions across different conditions
- **How**: List all desired phenotypes in `GI_PHENOTYPES`
- **Benefit**: Comprehensive analysis without re-running pipeline

### 6. ⚖️ **Differential Interaction Analysis**
- **Feature**: Compare GI scores between any two phenotypes
- **Why**: Identify condition-specific interactions
- **How**: Define comparisons in `DIFFERENTIAL_COMPARISONS`
- **Examples**:
  ```yaml
  DIFFERENTIAL_COMPARISONS:
    DrugSpecific: [Drug, Control]
    DoseDependent: [HighDose, LowDose]
    Temporal: [Late, Early]
  ```
- **Benefit**: Unlimited comparison types

### 7. 🎯 **Sequential Correlation Filtering**
- **Feature**: Apply correlation filter in specified order
- **Why**: Quality control for dual-guide experiments
- **How**: List phenotypes in `NO_CORRELATION_FILTER`
- **Detail**: 
  - Checks if dual-guide phenotype correlates with single guides
  - Filters cascade (later filters see results of earlier ones)
  - Three modes for different stringency levels
- **Benefit**: Removes low-quality constructs systematically

### 8. 🔢 **High-Precision Numerical I/O**
- **Feature**: Custom RDS workspace saving/loading with extended precision
- **Why**: Prevent precision loss in intermediate steps
- **How**: Automatic use of `r_precise_io.R` functions
- **Details**:
  - Sets `digits = 17` and `scipen = 999`
  - Uses uncompressed RDS (version 3)
  - Preserves numerical accuracy across pipeline stages
- **Benefit**: Reproducible results to full floating-point precision

### 9. 🧹 **Automatic Temporary File Cleanup**
- **Feature**: Intermediate RDS workspaces marked as temporary
- **Why**: Reduce disk space usage
- **How**: Snakemake's `temp()` directive
- **Benefit**: Final results retained, intermediate files cleaned automatically

### 10. 📦 **Integrated Package Management**
- **Feature**: Conda environments specified per rule
- **Why**: Ensures reproducible software environments
- **How**: Each rule has `conda: "../envs/smk-env.yaml"`
- **Usage**: Run with `--software-deployment-method conda`
- **Benefit**: Automatic environment creation and management

### 11. 🔍 **Automatic Metadata Generation**
- **Feature**: Generate construct annotations if not provided
- **Why**: Simplifies input requirements
- **How**: `prepare_counts_with_metadata.R` detects and generates
- **Benefit**: Works with minimal input files OR fully annotated files

### 12. 📊 **Multi-Level Score Reporting**
- **Feature**: Scores at construct, gene, and discriminant levels
- **Why**: Different analyses require different granularities
- **How**: Automatic generation at all levels
- **Benefit**: 
  - Construct level: Most detailed, for specific guide analysis
  - Gene level: Aggregated, more robust for biological interpretation
  - Discriminant: Normalized, for statistical hit calling

### 13. ⚙️ **Orientation Independence**
- **Feature**: Automatic averaging of forward/reverse orientations
- **Why**: Reduces technical variation in dual-guide experiments
- **How**: `create_averaged_phenotypes.R`
- **Benefit**: More robust GI score estimates

### 14. 📋 **Comprehensive Logging**
- **Feature**: Detailed logs for every rule, per screen
- **Why**: Debugging and transparency
- **How**: Dual logging to console and file (`dual_logging.R`)
- **Structure**: `logs/{screen}/{stage}/{screen}_{rule}.log`
- **Benefit**: Easy troubleshooting and auditing

### 15. 🚀 **Parallel Execution**
- **Feature**: Process multiple screens, replicates, and phenotypes in parallel
- **Why**: Reduce total runtime
- **How**: Snakemake's DAG-based execution with `--cores`
- **Benefit**: Efficient use of computational resources

---

## Output Structure

The pipeline generates organized outputs in `outputs/` directory:

```
outputs/
├── annotations/              # ID mapping files
│   ├── {screen}_construct_id_map.tsv
│   ├── {screen}_gene_combination_id_map.tsv
│   ├── {screen}_guide_combination_id_map.tsv
│   └── {screen}_pseudogene_id_map.tsv
│
├── counts/                   # Validated counts with metadata
│   └── {screen}_counts_with_metadata.tsv
│
├── phenotypes/               # Calculated phenotypes at various stages
│   ├── {screen}_processed_phenotypes.tsv
│   ├── {screen}_orientation_independent_phenotypes.tsv
│   ├── {screen}_single_sgRNA_phenotypes.tsv
│   └── {screen}_filtered_{phenotype}_phenotypes.tsv
│
├── gi_scores/                # Genetic interaction scores
│   └── {screen}/
│       ├── construct_scores/     # Per-construct GI scores
│       │   └── all_gis_{score}.tsv
│       ├── gene_combination_scores/  # Aggregated gene-level scores
│       │   └── gene_combination_scores_{score}.tsv
│       ├── discriminant_scores/  # Normalized scores + hits
│       │   ├── discriminant_scores_{score}.tsv
│       │   └── discriminant_hits_{score}.tsv
│       ├── individual_scores/    # Per-gene-pair directories
│       │   └── {score}/
│       ├── models/               # Model statistics
│       │   ├── model_estimates_{score}.tsv
│       │   └── model_stats_{score}.tsv
│       └── clusters/             # Diagnostic plots
│           └── diagnostic_plot_{score}.svg
│
├── misc_results/             # QC and filtering results
│   ├── {screen}_filter_flags.tsv
│   ├── {screen}_full_filter_flags.tsv
│   ├── {screen}_correlation_results.tsv
│   └── {screen}_correlation_summary.tsv
│
└── manuscript_figures/       # Publication figures (optional)
    └── figure_*.{svg,png}
```

**Key Files to Examine:**

1. **`{screen}_counts_with_metadata.tsv`**: Starting point with all metadata
2. **`{screen}_full_filter_flags.tsv`**: Which constructs passed/failed each filter
3. **`{screen}_correlation_results.tsv`**: Detailed correlation QC metrics
4. **`all_gis_{score}.tsv`**: Complete construct-level GI scores
5. **`gene_combination_scores_{score}.tsv`**: Gene-level aggregated scores
6. **`discriminant_hits_{score}.tsv`**: Final hit calls (synergy/suppression)

---

## Troubleshooting

### Common Issues

#### 1. **ZIP file not found or empty**
**Symptom**: Error reading counts file
**Solutions**:
- Check `COUNTS_FILEPATH_PATTERN` in config
- Ensure `{screen}` placeholder is used correctly
- Verify ZIP contains `*_counts.tsv` file
- Try with uncompressed TSV for testing

#### 2. **Missing metadata columns**
**Symptom**: Validation error about missing columns
**Solutions**:
- Let pipeline auto-generate metadata (recommended)
- If providing metadata, ensure all `EXPECTED_META_COLUMNS` are present
- Check column name spelling and capitalization

#### 3. **Doublings not defined**
**Symptom**: Error in normalization step
**Solutions**:
- Define `{SCREEN}_DOUBLINGS` for each screen
- Treatment names in doublings must match `COUNT_PREFIXES`
- Set `NORMALIZE: FALSE` if you don't have doubling information

#### 4. **No outputs generated**
**Symptom**: Pipeline completes but no files
**Solutions**:
- Check `SCREENS` parameter matches your file prefixes
- Verify `GI_PHENOTYPES` is not empty
- Run with `snakemake -n` to see what would be executed

#### 5. **Phenotype calculation fails**
**Symptom**: Error in `calculate_phenotypes` step
**Solutions**:
- Verify column naming matches `COUNT_PREFIXES` + `REPLICATES`
- Check `PHENOTYPE_PREFIX_MAP` uses valid prefix names
- Ensure count columns exist for all prefix/replicate combinations

#### 6. **Memory issues**
**Symptom**: R process killed or out-of-memory error
**Solutions**:
- Reduce number of cores used (lower `--cores` value)
- Process fewer screens simultaneously
- Consider processing screens separately
- Increase system memory allocation

#### 7. **Hit calling produces unexpected results**
**Symptom**: Too many or too few hits
**Solutions**:
- Adjust `HIT_THRESHOLDS` (lower = more hits, higher = fewer hits)
- Check discriminant score distributions in TSV files
- Verify control constructs are properly labeled
- Review correlation filtering results

### Validation Checklist

Before running the pipeline, verify:

- [ ] All screens listed in `SCREENS` have corresponding count files
- [ ] `COUNT_PREFIXES` matches your column naming scheme
- [ ] `REPLICATES` matches your experimental design
- [ ] `PHENOTYPE_PREFIX_MAP` uses valid `COUNT_PREFIXES`
- [ ] `GI_PHENOTYPES` is a subset of `PHENOTYPE_PREFIX_MAP` keys
- [ ] `DOUBLINGS` defined for each screen (if `NORMALIZE: TRUE`)
- [ ] `DIFFERENTIAL_COMPARISONS` uses phenotypes from `GI_PHENOTYPES`

### Testing the Pipeline

**Quick test run:**
```bash
# Dry run to see execution plan
snakemake -n --snakefile workflow/Snakefile --configfile config/config.yaml process_all

# Run just metadata generation for one screen
snakemake outputs/counts/screen2022_counts_with_metadata.tsv --cores 1

# Test complete preprocessing for one screen
snakemake outputs/phenotypes/screen2022_processed_phenotypes.tsv --cores 2

# Full pipeline (all screens)
snakemake --cores 6 --snakefile workflow/Snakefile --configfile config/config.yaml process_all
```

### Getting Help

1. **Check log files**: `logs/{screen}/{stage}/{rule}.log`
2. **Examine intermediate outputs**: Look at TSV files at each stage
3. **Run with verbose mode**: Add `-p` flag to Snakemake command
4. **Test individual rules**: Target specific output files
5. **Review R script**: Check corresponding script in `workflow/scripts/`

---

## Summary

This pipeline provides a **complete, flexible, and reproducible** workflow for genetic interaction analysis. Key strengths:

✅ **Automated**: Minimal manual intervention required
✅ **Adaptable**: Extensive configuration options for different experimental designs  
✅ **Robust**: Comprehensive quality control and filtering
✅ **Precise**: High-precision numerical calculations maintained throughout
✅ **Scalable**: Process multiple screens in parallel
✅ **Transparent**: Detailed logging and intermediate outputs
✅ **Reproducible**: Conda integration and version-controlled workflows

---

*Last updated: 2025-11-04*
