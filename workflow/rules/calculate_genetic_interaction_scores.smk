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
        targets += expand("../outputs/gi_scores/{screen}/differential_scores/differential_scores_{rep}.tsv",
                          screen=sc, rep=reps)
    return targets

rule compute_genetic_interaction_scores:
    input:
        input_orientation_indep_phenotypes="../outputs/phenotypes/{screen}_filtered_phenotypes.tsv",
        input_single_sgRNA_phenotypes="../outputs/phenotypes/{screen}_filtered_single_sgRNA_phenotypes.tsv"
    params:
        # the single score for this job
        score=lambda wildcards: wildcards.score
    # Use a separate job per score (avoid using functions in outputs). Snakemake will create
    # one job for each score by expanding the {score} wildcard elsewhere in the workflow.
    output:
        output_dir=directory("../outputs/gi_scores/{screen}/individual_scores/{score}"),
        output_all_scores="../outputs/gi_scores/{screen}/construct_scores/all_gis_{score}.tsv",
        output_model_estimates="../outputs/gi_scores/{screen}/models/model_estimates_{score}.tsv",
        output_model_stats="../outputs/gi_scores/{screen}/models/model_stats_{score}.tsv"
    script:
        "../scripts/calculate_gi_scores.R"

rule calculate_gene_level_scores:
    input:
        input_gi_scores="../outputs/gi_scores/{screen}/construct_scores/all_gis_{score}.tsv",
        input_idmap="data/annotations/{screen}_id_to_name_mapping.tsv"
    params:
        score=lambda wildcards: wildcards.score,
        screen=lambda wildcards: wildcards.screen
    output:
        output_gene_level_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv"
    script:
        "../scripts/calculate_gene_level_scores.R"

rule calculate_discriminant_scores:
    input:
        input_gi_scores="../outputs/gi_scores/{screen}/construct_scores/all_gis_{score}.tsv",
        input_gene_level_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv"
    params:
        score=lambda wildcards: wildcards.score,
        screen=lambda wildcards: wildcards.screen
    output:
        output_discriminant_scores="../outputs/gi_scores/{screen}/discriminant_scores/discriminant_scores_{score}.tsv"
    script:
        "../scripts/calculate_discriminant_scores.R"

rule calculate_differential_scores:
    input:
        input_gamma_gi_scores="../outputs/gi_scores/{screen}/construct_scores/all_gis_Gamma.{rep}.tsv",
        input_tau_gi_scores="../outputs/gi_scores/{screen}/construct_scores/all_gis_Tau.{rep}.tsv",
        input_idmap="data/annotations/{screen}_id_to_name_mapping.tsv"
    output:
        output_differential_scores="../outputs/gi_scores/{screen}/differential_scores/differential_scores_{rep}.tsv",
        output_gene_differential_scores="../outputs/gi_scores/{screen}/differential_scores/gene_differential_scores_{rep}.tsv",
        output_discriminant_differential_scores="../outputs/gi_scores/{screen}/differential_scores/discriminant_differential_scores_{rep}.tsv"
    script:
        "../scripts/calculate_differential_scores.R"

rule compute_all_genetic_interaction_scores:
    # wrapper rule to build all scores for a given screen; call with --config screen=<screen>
    input:
        lambda wildcards: _expand_scores_for_screen(wildcards)

rule compute_all_gene_level_scores:
    input:
        lambda wildcards: _expand_gene_level_score_targets(wildcards)

rule compute_all_discriminant_scores:
    input:
        lambda wildcards: _expand_discriminant_score_targets(wildcards)

rule compute_all_differential_scores:
    input:
        lambda wildcards: _expand_differential_score_targets(wildcards)

#rule call_hits: