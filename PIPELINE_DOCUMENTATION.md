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

The pipeline is controlled by `config/config.yaml`. This section explains how to set up your configuration.

### Essential Settings

#### 1. Directory Configuration
```yaml
OUTPUTS_DIR: ../outputs  # Where all results are written
LOGS_DIR: ../logs        # Where log files are stored
```
💡 **Tip**: Use relative paths from the `workflow/` directory or absolute paths.

#### 2. Screen Selection
```yaml
SCREENS: [screen2022, screen2023]
```
- List all screens to process
- Names must match filename prefixes in your counts files
- Pipeline processes all listed screens in parallel

#### 3. Input Counts Location
```yaml
COUNTS_FILEPATH_PATTERN: data/counts/{screen}_raw_counts.zip
```
- Use `{screen}` as placeholder (automatically replaced with each screen name)
- Supports both `.tsv` and `.zip` formats
- **Cool Feature**: ZIP files are automatically extracted - the pipeline finds `*_counts.tsv` files inside

Example patterns:
```yaml
# For TSV files directly:
COUNTS_FILEPATH_PATTERN: data/counts/{screen}_raw_counts.tsv

# For ZIP files (recommended for large files):
COUNTS_FILEPATH_PATTERN: data/counts/{screen}_raw_counts.zip

# For different directory structure:
COUNTS_FILEPATH_PATTERN: /path/to/data/{screen}/{screen}_counts.zip
```

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

#### 5. Phenotype Definitions
```yaml
PHENOTYPE_PREFIX_MAP:
  Gamma: [DMSO, T0]      # Control phenotype: DMSO / T0
  Tau: [TREATED, T0]     # Treatment phenotype: TREATED / T0
  Rho: [TREATED, DMSO]   # Relative phenotype: TREATED / DMSO
```
- **Format**: `Phenotype_Name: [Numerator, Denominator]`
- Prefixes must match those in `COUNT_PREFIXES`
- You can define as many phenotypes as needed
- Names can be customized (not limited to Gamma/Tau/Rho)

#### 6. GI Phenotype Selection
```yaml
GI_PHENOTYPES: [Gamma, Tau]
```
- Specifies which phenotypes to use for genetic interaction scoring
- Must be a subset of keys in `PHENOTYPE_PREFIX_MAP`
- Reduces computation time if you don't need all phenotypes

### Screen-Specific Parameters

**🌟 Adaptability Feature**: Override global settings for individual screens by prefixing parameters with the screen name.

#### Population Doublings
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
- Used for normalization when `NORMALIZE: TRUE`
- Define per-treatment, per-replicate
- Treatment names must match your `COUNT_PREFIXES`

#### Screen-Specific Overrides
Any parameter can be made screen-specific:

```yaml
# Global default:
COUNT_PREFIXES: [T0, DMSO, TREATED]

# Screen-specific override:
SCREEN2023_COUNT_PREFIXES: [T0, CONTROL, EXPERIMENTAL]
SCREEN2023_PHENOTYPE_PREFIX_MAP:
  Gamma: [CONTROL, T0]
  Tau: [EXPERIMENTAL, T0]
```

**Example use cases:**
- Different treatment names across screens
- Different replicate structures
- Screen-specific thresholds
- Varied phenotype definitions

### Filter Parameters

#### Count Quality Filters
```yaml
INDIVIDUAL_SGRNA_MEDIAN_THRESHOLD: 35
COMBINATION_SGRNA_COUNT_THRESHOLD: 50
```
- `INDIVIDUAL_SGRNA_MEDIAN_THRESHOLD`: Minimum median count across replicates for single guides
- `COMBINATION_SGRNA_COUNT_THRESHOLD`: Minimum count threshold for guide combinations
- Constructs failing these filters are flagged and removed for genetic interaction score calculations

#### Phenotype Processing
```yaml
PSEUDOCOUNT: 10              # Added before log transformation
NORMALIZE: TRUE              # Normalize by population doublings
AVERAGE_REPLICATES: TRUE     # Generate averaged phenotype columns
```

#### Correlation Filtering
```yaml
NO_CORRELATION_FILTER:
  - Gamma
  - Tau

NO_CORRELATION_THRESHOLD: 0.25
CORRELATION_FILTER_MODE: avg_only  # Options: avg_only, any, all
```
**🎯 Smart Feature**: Sequential correlation filtering
- Applied in the order specified in `NO_CORRELATION_FILTER`
- Checks correlation between dual guide combinations and their constituent single guides for each sgRNA
- Three modes:
  - `avg_only`: Filter based on average correlation
  - `any`: Filter if ANY single replicate (or average, if present) shows poor correlation
  - `all`: Filter only if ALL replicates (and average, if present) show poor correlation

### Advanced Settings

#### Differential Comparisons
```yaml
DIFFERENTIAL_COMPARISONS:
  Nu: [Tau, Gamma]  # Nu = Tau - Gamma (differential interaction)
```
- Define which phenotypes to compare
- **Format**: `Name: [Phenotype1, Phenotype2]` calculates Phenotype1 - Phenotype2
- Can define multiple comparisons

**Example - Multiple comparisons:**
```yaml
DIFFERENTIAL_COMPARISONS:
  Nu: [Tau, Gamma]
  Delta: [Rho, Gamma]
  Epsilon: [Tau, Rho]
```

#### Hit Calling Thresholds
```yaml
HIT_THRESHOLDS:
  Gamma: 0.995    # 99.5th percentile for synergy, 0.5th for suppression
  Tau: 0.995
  Rho: 0.995
  Nu: 1           # 100th percentile (most extreme) for differential
```
- Thresholds represent percentile cutoffs of discriminant scores based on control distributions
- Values between 0 and 1
- Higher values = more stringent (fewer hits)
- Can be screen-specific: `SCREEN2022_HIT_THRESHOLDS`

#### Clustering Parameters
```yaml
PHENOTYPES_TO_CLUSTER: [Gamma.OI.Avg, Tau.OI.Avg, Nu.OI.Avg]
```
- Generates diagnostic clustering plots
- Use orientation-independent averaged scores for best results

#### Internal Metadata Columns (⚠️ DANGER ZONE)
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

**Purpose**: Validate and enrich counts files with construct metadata.

**Rules:**
1. `prepare_counts_with_metadata`
2. `generate_id_maps`
3. `validate_counts`

**What Happens:**
- Reads raw counts (auto-detects ZIP vs TSV)
- Checks for required metadata columns
- **Auto-generates** metadata if not present (construct IDs, gene combinations, etc.)
- Creates ID mapping files for downstream use
- Validates column structure and data integrity

**Inputs:**
- Raw counts files (ZIP or TSV)

**Outputs:**
- `outputs/counts/{screen}_counts_with_metadata.tsv`
- `outputs/annotations/{screen}_*_id_map.tsv` (4 types)

**🔧 Adaptability**: If your counts already have metadata, the pipeline detects and uses them. If not, it generates them automatically.

---

### Stage 2: Preprocessing

**Purpose**: Filter low-quality constructs and calculate normalized phenotypes.

**Rules:**
1. `apply_count_filters`
2. `calculate_phenotypes`
3. `normalize_phenotypes`
4. `create_averaged_phenotypes`
5. `apply_correlation_filter`

**Workflow:**

```
Counts + Metadata
    ↓
Apply Count Filters → Flag low-count constructs
    ↓
Calculate Phenotypes → Compute log2 ratios with pseudocounts
    ↓
Normalize Phenotypes → Adjust for population doublings (optional)
    ↓
Create Averaged Phenotypes → Generate replicate averages + orientation-independent
    ↓
Apply Correlation Filter → QC based on single-guide correlation
```

**Key Features:**

**Count Filtering:**
- Individual sgRNA median counts across replicates
- Combination sgRNA minimum counts
- Flags stored but data retained for transparency

**Phenotype Calculation:**
- Log2 transformation: `log2((count + pseudocount))`
- Ratio calculation: `log2(numerator) - log2(denominator)`
- Per-replicate and averaged

**Normalization:**
- Adjusts for different growth rates between treatments
- Formula: `phenotype / population_doublings`
- Optional (controlled by `NORMALIZE` parameter)

**Orientation Independence:**
- Averages forward and reverse orientations
- Reduces technical variation
- Key for symmetric genetic interactions

**Correlation Filtering:**
- Compares dual-guide phenotypes to constituent single guides
- Three correlation modes (avg_only, any, all)
- Sequential application (filters compound)
- Generates detailed correlation reports

**Outputs:**
- `outputs/phenotypes/{screen}_processed_phenotypes.tsv`
- `outputs/phenotypes/{screen}_orientation_independent_phenotypes.tsv`
- `outputs/phenotypes/{screen}_single_sgRNA_phenotypes.tsv`
- `outputs/phenotypes/{screen}_filtered_{phenotype}_phenotypes.tsv` (per phenotype in filter list)
- `outputs/misc_results/{screen}_filter_flags.tsv`
- `outputs/misc_results/{screen}_full_filter_flags.tsv`
- `outputs/misc_results/{screen}_correlation_results.tsv`
- `outputs/misc_results/{screen}_correlation_summary.tsv`

---

### Stage 3: GI Score Calculation

**Purpose**: Compute genetic interaction scores comparing observed vs expected phenotypes.

**Rules:**
1. `compute_genetic_interaction_scores`
2. `calculate_gene_level_scores`
3. `calculate_discriminant_scores`

**Method:**

For each construct, the pipeline:
1. **Fits a linear model**: Compares dual-guide phenotype to sum of single-guide phenotypes
2. **Calculates GI score**: Deviation from additivity (residual from model)
3. **Aggregates to gene level**: Median across all constructs targeting same gene pair
4. **Computes discriminant**: Normalizes GI scores for hit calling

**Mathematical Framework:**

```
Expected (additive) = sgRNA_A_phenotype + sgRNA_B_phenotype
Observed = dual_guide_phenotype
GI Score = Observed - Expected

Discriminant = (GI Score - median(negative controls)) / MAD(negative controls)
```

**Scoring Hierarchy:**

```
Construct Level → Individual dual-guide combinations
    ↓
Gene Level → Aggregated across all constructs per gene pair
    ↓
Discriminant → Normalized scores for statistical inference
```

**🎯 Cool Feature**: The pipeline automatically handles:
- Per-replicate scores (e.g., `Gamma.OI.R1`, `Gamma.OI.R2`)
- Averaged scores (e.g., `Gamma.OI.Avg`)
- Multiple phenotypes in parallel

**Outputs:**
- `outputs/gi_scores/{screen}/construct_scores/all_gis_{score}.tsv`
- `outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv`
- `outputs/gi_scores/{screen}/discriminant_scores/discriminant_scores_{score}.tsv`
- `outputs/gi_scores/{screen}/individual_scores/{score}/` (directory with per-gene-pair files)
- `outputs/gi_scores/{screen}/models/model_estimates_{score}.tsv`
- `outputs/gi_scores/{screen}/models/model_stats_{score}.tsv`

---

### Stage 4: Differential Analysis

**Purpose**: Identify condition-specific interactions by comparing GI scores across phenotypes.

**Rule:**
- `calculate_differential_scores`

**Method:**

Calculates differential interactions (e.g., Nu = Tau - Gamma):
- Compares treatment-specific GI scores (Tau) to control GI scores (Gamma)
- Identifies interactions that change between conditions
- Produces construct-level, gene-level, and discriminant scores

**Use Cases:**
- Drug-specific interactions (treatment vs control)
- Context-dependent interactions
- Genetic background effects

**Example:**
```
Gamma GI (DMSO):    Gene A + Gene B = -0.5 (suppression)
Tau GI (Treatment): Gene A + Gene B = +1.2 (synergy)
Nu Differential:    1.2 - (-0.5) = 1.7 (condition-specific synergy)
```

**🌟 Adaptability**: Define unlimited differential comparisons in config:
```yaml
DIFFERENTIAL_COMPARISONS:
  Nu: [Tau, Gamma]           # Treatment - Control
  Delta: [Tau2, Tau1]        # Dose comparison
  CustomName: [PhenoA, PhenoB]  # Any comparison you define
```

**Outputs:**
- `outputs/gi_scores/{screen}/construct_scores/all_gis_{comparison}.{rep}.tsv`
- `outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{comparison}.{rep}.tsv`
- `outputs/gi_scores/{screen}/discriminant_scores/discriminant_scores_{comparison}.{rep}.tsv`

---

### Stage 5: Hit Calling

**Purpose**: Identify statistically significant genetic interactions.

**Rule:**
- `call_hits`

**Method:**

Uses discriminant scores and percentile thresholds:
1. Ranks all constructs/gene pairs by discriminant score
2. Applies threshold to identify hits
3. Classifies interactions:
   - **Synergy**: Discriminant > upper threshold (more negative than expected)
   - **Suppression**: Discriminant < lower threshold (more positive than expected)
   - **No Interaction**: Within threshold bounds

**Threshold Interpretation:**
```yaml
HIT_THRESHOLDS:
  Gamma: 0.995  # Top and bottom 0.5% called as hits
  Nu: 1.0       # Only most extreme values (100th percentile)
```

**Classification Example:**
- Threshold = 0.995
- Upper cutoff = 99.5th percentile of discriminant
- Lower cutoff = 0.5th percentile of discriminant
- Hits are outside these bounds

**🎯 Smart Feature**: Screen-specific and phenotype-specific thresholds allow fine-tuned control over hit stringency.

**Outputs:**
- `outputs/gi_scores/{screen}/discriminant_scores/discriminant_hits_{score}.tsv`

**Hit File Contents:**
- All discriminant score columns
- `hit` column: "Synergy", "Suppression", or "No Interaction"
- `threshold_used`: Which threshold was applied

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
- [ ] Output directories exist and have write permissions

### Testing the Pipeline

**Quick test run:**
```bash
# Dry run to see execution plan
snakemake -n --snakefile workflow/Snakefile --configfile config/config.yaml

# Run just metadata generation for one screen
snakemake outputs/counts/screen2022_counts_with_metadata.tsv --cores 1

# Test complete preprocessing for one screen
snakemake outputs/phenotypes/screen2022_processed_phenotypes.tsv --cores 2

# Full pipeline (all screens)
snakemake --cores 4 --snakefile workflow/Snakefile --configfile config/config.yaml
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

**Getting Started:**
1. Edit `config/config.yaml` following the [Configuration Guide](#configuration-guide)
2. Place count files in `workflow/data/counts/`
3. Run: `snakemake --cores 4 --snakefile workflow/Snakefile --configfile config/config.yaml`
4. Find results in `outputs/` directory

For questions or issues, consult the [Troubleshooting](#troubleshooting) section or examine pipeline logs.

---

*Last updated: 2025-11-03*
