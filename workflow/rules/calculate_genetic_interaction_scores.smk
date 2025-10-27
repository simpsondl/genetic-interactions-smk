def _expand_scores_for_screen(wc=None):
    sc = config.get("screen") if isinstance(config, dict) else None
    if sc is None:
        raise Exception("Please provide the target screen via --config screen=<screen>")
    return expand("../outputs/gi_scores/{screen}/individual_scores/{score}", 
                  screen=sc, score=config[f"{sc.upper()}_GI_SCORES"])


rule compute_genetic_interaction_scores:
    input:
        input_orientation_indep_phenotypes="../outputs/phenotypes/{screen}_filtered_phenotypes.tsv",
        input_single_sgRNA_phenotypes="../outputs/phenotypes/{screen}_filtered_single_sgRNA_phenotypes.tsv"
    # Use a separate job per score (avoid using functions in outputs). Snakemake will create
    # one job for each score by expanding the {score} wildcard elsewhere in the workflow.
    output:
        out_dir=directory("../outputs/gi_scores/{screen}/individual_scores/{score}")
    params:
        # the single score for this job
        score=lambda wildcards: wildcards.score,
    script:
        "../scripts/calculate_gi_scores.R"


rule compute_all_genetic_interaction_scores:
    # wrapper rule to build all scores for a given screen; call with --config screen=<screen>
    input:
        lambda wildcards: _expand_scores_for_screen(wildcards)

rule calculate_gene_level_scores:
    input:
        input_gi_scores="../outputs/gi_scores/{screen}/all_scores/all_gis_{score}.tsv",
        input_idmap="../data/annotation/{screen}_id_to_name_mapping.tsv"
    output:
        output_gene_level_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv"
    script:
        "../scripts/calculate_gene_level_scores.R"

rule calculate_discriminant_scores:
    input:
        input_gi_scores="../outputs/gi_scores/{screen}/all_scores/all_gis_{score}.tsv",
        input_gene_level_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv"
    output:
        output_discriminant_scores="../outputs/gi_scores/{screen}/discriminant_scores/discriminant_scores_{score}.tsv"
    script:
        "../scripts/calculate_discriminant_scores.R"

rule calculate_differential_scores:
    input:
        input_gamma_gi_scores="../outputs/gi_scores/{screen}/all_scores/all_gis_Gamma.OI.{rep}.tsv",
        input_tau_gi_scores="../outputs/gi_scores/{screen}/all_scores/all_gis_Tau.OI.{rep}.tsv",
        input_idmap="../data/annotation/{screen}_id_to_name_mapping.tsv"
    output:
        output_differential_scores="../outputs/gi_scores/{screen}/differential_scores/differential_scores_{rep}.tsv",
        output_gene_differential_scores="../outputs/gi_scores/{screen}/differential_scores/gene_differential_scores_{rep}.tsv"
    script:
        "../scripts/calculate_differential_scores.R"

#rule call_hits: