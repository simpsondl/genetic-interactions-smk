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
        return f"../outputs/counts/{wildcards.screen}_counts_with_metadata.tsv"
    elif mode == "counts_metadata":
        # Allow the user to provide a pattern; substitute {screen}
        try:
            return provided_pattern.format(screen=wildcards.screen)
        except Exception:
            # fallback to the conventional path
            return f"../outputs/counts/{wildcards.screen}_counts_with_metadata.tsv"
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
    for sc in screens:
        targets += expand("../outputs/gi_scores/{screen}/individual_scores/{score}",
                          screen=sc, score=config[f"{sc.upper()}_GI_SCORES"])
    return targets


def _expand_gene_level_score_targets(wc=None):
    screens = config.get("SCREENS")
    if screens is None:
        screens = [k.lower().split("_GI_SCORES")[0].lower() for k in config if k.endswith("_GI_SCORES")]
    targets = []
    for sc in screens:
        targets += expand("../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv",
                          screen=sc, score=config[f"{sc.upper()}_GI_SCORES"])
    return targets


def _expand_discriminant_score_targets(wc=None):
    screens = config.get("SCREENS")
    if screens is None:
        screens = [k.lower().split("_GI_SCORES")[0].lower() for k in config if k.endswith("_GI_SCORES")]
    targets = []
    for sc in screens:
        targets += expand("../outputs/gi_scores/{screen}/discriminant_scores/discriminant_scores_{score}.tsv",
                          screen=sc, score=config[f"{sc.upper()}_GI_SCORES"])
    return targets
    

def _expand_differential_score_targets(wc=None):
    screens = config.get("SCREENS")
    if screens is None:
        # same fallback as above
        screens = [k.lower().split("_GI_SCORES")[0].lower() for k in config if k.endswith("_GI_SCORES")]
    reps = config.get("DIFFERENTIAL_SCORES")
    if reps is None:
        raise Exception("Please add a DIFFERENTIAL_SCORES list to config.yaml, e.g. DIFFERENTIAL_SCORES: [OI.R1, OI.R2]")
    targets = []
    for sc in screens:
        # differential construct scores (Nu.*) were moved into construct_scores
        targets += expand("../outputs/gi_scores/{screen}/construct_scores/all_gis_Nu.{rep}.tsv",
                          screen=sc, rep=reps)
    return targets


def _expand_hit_targets(wc=None):
    screens = config.get("SCREENS")
    if screens is None:
        screens = [k.lower().split("_GI_SCORES")[0].lower() for k in config if k.endswith("_GI_SCORES")]
    reps = config.get("DIFFERENTIAL_SCORES", [])
    targets = []
    for sc in screens:
        # Gamma/Tau
        targets += expand("../outputs/gi_scores/{screen}/discriminant_scores/discriminant_hits_{score}.tsv",
                          screen=sc, score=config[f"{sc.upper()}_GI_SCORES"])
        # Nu
        if reps:
            targets += expand("../outputs/gi_scores/{screen}/discriminant_scores/discriminant_hits_Nu.{rep}.tsv",
                              screen=sc, rep=reps)
    return targets


def _expand_differential_hit_targets(wc=None):
    screens = config.get("SCREENS")
    if screens is None:
        # same fallback as above
        screens = [k.lower().split("_GI_SCORES")[0].lower() for k in config if k.endswith("_GI_SCORES")]
    reps = config.get("DIFFERENTIAL_SCORES")
    if reps is None:
        raise Exception("Please add a DIFFERENTIAL_SCORES list to config.yaml, e.g. DIFFERENTIAL_SCORES: [OI.R1, OI.R2]")
    targets = []
    for sc in screens:
        # differential hit outputs are now produced as discriminant_hits_Nu.{rep}.tsv
        targets += expand("../outputs/gi_scores/{screen}/discriminant_scores/discriminant_hits_Nu.{rep}.tsv",
                          screen=sc, rep=reps)
    return targets


def _expand_diagnostic_plots(wc=None):
    screens = config.get("SCREENS")
    if screens is None:
        screens = [k.lower().split("_GI_SCORES")[0].lower() for k in config if k.endswith("_GI_SCORES")]
    targets = []
    for sc in screens:
        targets += expand("../outputs/gi_scores/{screen}/clusters/diagnostic_plot_{score}.svg",
                          screen=sc, score=config.get("PHENOTYPES_TO_CLUSTER"))
    return targets


def _gi_orientation_indep_phenotype(wildcards):
    """Return the orientation-independent phenotype file path for Gamma/Tau."""
    sc = wildcards.screen
    score = str(wildcards.score)
    if score.startswith("Gamma"):
        return f"../outputs/phenotypes/{sc}_filtered_gamma_phenotypes.tsv"
    if score.startswith("Tau"):
        return f"../outputs/phenotypes/{sc}_filtered_tau_phenotypes.tsv"
    raise ValueError(f"Unsupported score wildcard for GI phenotype: {score}")


def _gi_single_sgRNA_phenotype(wildcards):
    """Return the single-sgRNA phenotype file path for Gamma/Tau."""
    sc = wildcards.screen
    score = str(wildcards.score)
    if score.startswith("Gamma"):
        return f"../outputs/phenotypes/{sc}_filtered_gamma_single_sgRNA_phenotypes.tsv"
    if score.startswith("Tau"):
        return f"../outputs/phenotypes/{sc}_filtered_tau_single_sgRNA_phenotypes.tsv"
    raise ValueError(f"Unsupported score wildcard for GI single-sgRNA phenotype: {score}")