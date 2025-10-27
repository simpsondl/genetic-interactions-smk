def _expand_scores_for_screen(wc=None):
    sc = config.get("screen") if isinstance(config, dict) else None
    if sc is None:
        raise Exception("Please provide the target screen via --config screen=<screen>")
    return expand("../outputs/gi_scores/{screen}/individual_scores/{score}", screen=sc, score=config[f"{sc.upper()}_GI_SCORES"])


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