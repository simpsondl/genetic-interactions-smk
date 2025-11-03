import os


def _choose_counts(wildcards):
    """Choose the counts input according to config.

    Configuration key `COUNTS_SOURCE` controls behavior and accepts:
      - "manuscript" (default): use raw input counts under data/counts/{screen}_raw_counts.tsv or .zip
      - "counts_only": use the generated counts with metadata produced by the pipeline
        (../outputs/counts/{screen}_counts_with_metadata.tsv). This creates a DAG
        dependency on the `generate_id_maps` rule.
      - "counts_metadata": use a pre-provided counts-with-metadata file. The pattern to
        use may be supplied via `COUNTS_FILEPATH_PATTERN` in config.yaml;
        if not supplied, it defaults to ../outputs/counts/{screen}_counts_with_metadata.tsv
    """
    mode = config.get("COUNTS_SOURCE", "manuscript")
    mode = str(mode).lower()
    # Provided pattern (optional)
    provided_pattern = config.get("COUNTS_FILEPATH_PATTERN", "data/counts/{screen}_counts_with_metadata.tsv")

    if mode == "counts_only":
        # Depend on the pipeline's generated file (generate_id_maps rule)
        return f"{OUTPUTS_DIR}/counts/{wildcards.screen}_counts_with_metadata.tsv"
    elif mode == "counts_metadata":
        # Allow the user to provide a pattern; substitute {screen}
        try:
            return provided_pattern.format(screen=wildcards.screen)
        except Exception:
            # fallback to the conventional path
            return f"{OUTPUTS_DIR}/counts/{wildcards.screen}_counts_with_metadata.tsv"
    else:
        # manuscript mode (default) -- existing behaviour: prefer TSV then ZIP
        tsv = f"data/counts/{wildcards.screen}_raw_counts.tsv"
        zipf = f"data/counts/{wildcards.screen}_raw_counts.zip"
        if os.path.exists(tsv):
            return tsv
        if os.path.exists(zipf):
            return zipf
        raise FileNotFoundError(f"Neither {{tsv}} nor {{zipf}} found for screen {wildcards.screen}")


def _expand_scores_for_screen(wc=None):
    screens = config.get("SCREENS")
    if screens is None:
        # fallback: infer from config keys
        screens = [k.lower().split("_GI_SCORES")[0].lower() for k in config if k.endswith("_GI_SCORES")]

    targets = []
    # global phenotype/replicate settings
    global_phens = config.get("GI_PHENOTYPES")
    replicates = config.get("REPLICATES", []) or []
    average_reps = bool(config.get("AVERAGE_REPLICATES", False))

    for sc in screens:
        # Prefer explicit per-screen list if present (backwards compatible)
        per_screen_key = f"{sc.upper()}_GI_SCORES"
        per_screen = config.get(per_screen_key)
        if per_screen:
            scores = per_screen
        else:
            # Build from GI_PHENOTYPES + REPLICATES + AVERAGE_REPLICATES
            if not global_phens:
                raise Exception(
                    f"Please define either {per_screen_key} in config.yaml or a global GI_PHENOTYPES list"
                )
            scores = []
            for ph in global_phens:
                # orientation-independent (OI) is the convention used across the scripts
                for r in replicates:
                    scores.append(f"{ph}.OI.{r}")
                if average_reps:
                    scores.append(f"{ph}.OI.Avg")

        # Use expand to produce the filesystem targets; expand accepts a list for the 'score' wildcard
        targets += expand(f"{OUTPUTS_DIR}/gi_scores/{{screen}}/individual_scores/{{score}}",
                          screen=sc, score=scores)

    return targets


def _scores_for_screen(sc):
    """Return the list of score names for a given screen.

    Preference order:
    1. Per-screen config key (e.g., SCREEN2023_GI_SCORES)
    2. Built from global GI_PHENOTYPES, REPLICATES, and AVERAGE_REPLICATES
    """
    per_screen_key = f"{sc.upper()}_GI_SCORES"
    per_screen = config.get(per_screen_key)
    if per_screen:
        return per_screen

    global_phens = config.get("GI_PHENOTYPES")
    if not global_phens:
        raise Exception(f"Please define either {per_screen_key} in config.yaml or a global GI_PHENOTYPES list")

    replicates = config.get("REPLICATES", []) or []
    average_reps = bool(config.get("AVERAGE_REPLICATES", False))
    scores = []
    for ph in global_phens:
        for r in replicates:
            scores.append(f"{ph}.OI.{r}")
        if average_reps:
            scores.append(f"{ph}.OI.Avg")
    return scores


def _expand_gene_level_score_targets(wc=None):
    screens = config.get("SCREENS")
    if screens is None:
        screens = [k.lower().split("_GI_SCORES")[0].lower() for k in config if k.endswith("_GI_SCORES")]
    targets = []
    for sc in screens:
        scores = _scores_for_screen(sc)
        targets += expand(f"{OUTPUTS_DIR}/gi_scores/{{screen}}/gene_combination_scores/gene_combination_scores_{{score}}.tsv",
                          screen=sc, score=scores)
    return targets


def _expand_discriminant_score_targets(wc=None):
    screens = config.get("SCREENS")
    if screens is None:
        screens = [k.lower().split("_GI_SCORES")[0].lower() for k in config if k.endswith("_GI_SCORES")]
    targets = []
    for sc in screens:
        scores = _scores_for_screen(sc)
        targets += expand(f"{OUTPUTS_DIR}/gi_scores/{{screen}}/discriminant_scores/discriminant_scores_{{score}}.tsv",
                          screen=sc, score=scores)
    return targets
    

def _expand_differential_score_targets(wc=None):
    screens = config.get("SCREENS")
    if screens is None:
        # same fallback as above
        screens = [k.lower().split("_GI_SCORES")[0].lower() for k in config if k.endswith("_GI_SCORES")]
    # Allow configurable differential comparisons mapping in config.yaml. Default keeps Nu = [Tau, Gamma]
    differential_map = config.get("DIFFERENTIAL_COMPARISONS", {"Nu": ["Tau", "Gamma"]})
    comps = list(differential_map.keys())

    # Build differential replicate suffixes
    replicates = config.get("REPLICATES", []) or []
    average_reps = bool(config.get("AVERAGE_REPLICATES", False))
    oi_prefix = "OI"
    reps = [f"{oi_prefix}.{r}" for r in replicates]
    if average_reps:
        reps.append(f"{oi_prefix}.Avg")
    
    targets = []
    for sc in screens:
        # For each comparison key (e.g., Nu) expand the differential construct score outputs
        for comp in comps:
            targets += expand(f"{OUTPUTS_DIR}/gi_scores/{{screen}}/construct_scores/all_gis_{{comp}}.{{rep}}.tsv",
                              screen=sc, comp=comp, rep=reps)
    return targets


def _expand_hit_targets(wc=None):
    screens = config.get("SCREENS")
    if screens is None:
        screens = [k.lower().split("_GI_SCORES")[0].lower() for k in config if k.endswith("_GI_SCORES")]
    # determine differential reps (backwards compatible)
    diff_reps_cfg = config.get("DIFFERENTIAL_SCORES")
    if diff_reps_cfg is None:
        replicates = config.get("REPLICATES", []) or []
        average_reps = bool(config.get("AVERAGE_REPLICATES", False))
        oi_prefix = "OI"
        reps = [f"{oi_prefix}.{r}" for r in replicates]
        if average_reps:
            reps.append(f"{oi_prefix}.Avg")
    else:
        reps = diff_reps_cfg
    targets = []
    for sc in screens:
        # Gamma/Tau (or other phenotype-derived scores)
        scores = _scores_for_screen(sc)
        targets += expand(f"{OUTPUTS_DIR}/gi_scores/{{screen}}/discriminant_scores/discriminant_hits_{{score}}.tsv",
                          screen=sc, score=scores)
        # Differential comparisons: produce discriminant hit files for each configured
        # differential comparison key (previously hard-coded to 'Nu')
        differential_map = config.get("DIFFERENTIAL_COMPARISONS", {"Nu": ["Tau", "Gamma"]})
        comps = list(differential_map.keys())
        if reps:
            for comp in comps:
                targets += expand(f"{OUTPUTS_DIR}/gi_scores/{{screen}}/discriminant_scores/discriminant_hits_{{comp}}.{{rep}}.tsv",
                                  screen=sc, comp=comp, rep=reps)
    return targets


def _expand_diagnostic_plots(wc=None):
    screens = config.get("SCREENS")
    if screens is None:
        screens = [k.lower().split("_GI_SCORES")[0].lower() for k in config if k.endswith("_GI_SCORES")]
    targets = []
    for sc in screens:
        targets += expand(f"{OUTPUTS_DIR}/gi_scores/{{screen}}/clusters/diagnostic_plot_{{score}}.svg",
                          screen=sc, score=config.get("PHENOTYPES_TO_CLUSTER"))
    return targets


def _gi_orientation_indep_phenotype(wildcards):
    """Return the orientation-independent phenotype file path for any configured phenotype.
    
    Extracts the phenotype prefix from the score wildcard (e.g., "Gamma" from "Gamma.OI.R1")
    and constructs the path. For phenotypes in NO_CORRELATION_FILTER, returns the filtered
    file; otherwise returns the last filtered phenotype file based on NO_CORRELATION_FILTER order.
    Supports screen-specific NO_CORRELATION_FILTER configurations.
    """
    sc = wildcards.screen
    score = str(wildcards.score)
    # Extract the phenotype prefix (first part before the dot)
    pheno_prefix = score.split(".")[0]
    pheno_lower = pheno_prefix.lower()
    
    # Check if this phenotype has correlation filtering applied (screen-specific or global)
    no_corr_filter = config.get(f"{sc.upper()}_NO_CORRELATION_FILTER",
                                 config.get("NO_CORRELATION_FILTER", []))
    if pheno_prefix in no_corr_filter:
        # Use filtered file for this phenotype
        result = f"{OUTPUTS_DIR}/phenotypes/{sc}_filtered_{pheno_lower}_phenotypes.tsv"
    else:
        # Use the last filtered phenotype file (last item in NO_CORRELATION_FILTER list)
        if no_corr_filter:
            last_filtered_pheno = no_corr_filter[-1].lower()
            result = f"{OUTPUTS_DIR}/phenotypes/{sc}_filtered_{last_filtered_pheno}_phenotypes.tsv"
        else:
            # Fallback if NO_CORRELATION_FILTER is empty
            result = f"{OUTPUTS_DIR}/phenotypes/{sc}_orientation_independent_phenotypes.tsv"
    
    return result


def _gi_single_sgRNA_phenotype(wildcards):
    """Return the single-sgRNA phenotype file path for any configured phenotype.
    
    Extracts the phenotype prefix from the score wildcard (e.g., "Gamma" from "Gamma.OI.R1")
    and constructs the path. For phenotypes in NO_CORRELATION_FILTER, returns the filtered
    file; otherwise returns the last filtered single-sgRNA phenotype file based on NO_CORRELATION_FILTER order.
    Supports screen-specific NO_CORRELATION_FILTER configurations.
    """
    sc = wildcards.screen
    score = str(wildcards.score)
    # Extract the phenotype prefix (first part before the dot)
    pheno_prefix = score.split(".")[0]
    pheno_lower = pheno_prefix.lower()
    
    # Check if this phenotype has correlation filtering applied (screen-specific or global)
    no_corr_filter = config.get(f"{sc.upper()}_NO_CORRELATION_FILTER",
                                 config.get("NO_CORRELATION_FILTER", []))
    if pheno_prefix in no_corr_filter:
        # Use filtered file for this phenotype
        result = f"{OUTPUTS_DIR}/phenotypes/{sc}_filtered_{pheno_lower}_single_sgRNA_phenotypes.tsv"
    else:
        # Use the last filtered phenotype file (last item in NO_CORRELATION_FILTER list)
        if no_corr_filter:
            last_filtered_pheno = no_corr_filter[-1].lower()
            result = f"{OUTPUTS_DIR}/phenotypes/{sc}_filtered_{last_filtered_pheno}_single_sgRNA_phenotypes.tsv"
        else:
            # Fallback if NO_CORRELATION_FILTER is empty
            result = f"{OUTPUTS_DIR}/phenotypes/{sc}_single_sgRNA_phenotypes.tsv"
    
    return result