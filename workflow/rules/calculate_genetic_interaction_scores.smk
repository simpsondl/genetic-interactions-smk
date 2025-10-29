rule compute_genetic_interaction_scores:
    input:
        input_orientation_indep_phenotypes="../outputs/phenotypes/{screen}_filtered_phenotypes.tsv",
        input_single_sgRNA_phenotypes="../outputs/phenotypes/{screen}_filtered_single_sgRNA_phenotypes.tsv"
    output:
        output_dir=directory("../outputs/gi_scores/{screen}/individual_scores/{score}"),
        output_all_scores="../outputs/gi_scores/{screen}/construct_scores/all_gis_{score}.tsv",
        output_model_estimates="../outputs/gi_scores/{screen}/models/model_estimates_{score}.tsv",
        output_model_stats="../outputs/gi_scores/{screen}/models/model_stats_{score}.tsv"
    log:
        "../outputs/logs/{screen}/{screen}_{score}_compute_genetic_interaction_scores.log"
    params:
        score=lambda wildcards: wildcards.score
    script:
        "../scripts/calculate_gi_scores.R"

rule calculate_gene_level_scores:
    input:
        input_gi_scores="../outputs/gi_scores/{screen}/construct_scores/all_gis_{score}.tsv",
        input_idmap="data/annotations/{screen}_id_to_name_mapping.tsv"
    output:
        output_gene_level_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv"
    log:
        "../outputs/logs/{screen}/{screen}_{score}_calculate_gene_level_scores.log"
    params:
        score=lambda wildcards: wildcards.score,
        screen=lambda wildcards: wildcards.screen
    script:
        "../scripts/calculate_gene_level_scores.R"

rule calculate_discriminant_scores:
    input:
        input_gi_scores="../outputs/gi_scores/{screen}/construct_scores/all_gis_{score}.tsv",
        input_gene_level_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv"
    output:
        output_discriminant_scores="../outputs/gi_scores/{screen}/discriminant_scores/discriminant_scores_{score}.tsv"
    log:
        "../outputs/logs/{screen}/{screen}_{score}_calculate_discriminant_scores.log"
    params:
        score=lambda wildcards: wildcards.score,
        screen=lambda wildcards: wildcards.screen
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
    log:
        "../outputs/logs/{screen}/{screen}_{rep}_calculate_differential_scores.log"
    params:
        rep=lambda wildcards: wildcards.rep,
        screen=lambda wildcards: wildcards.screen
    script:
        "../scripts/calculate_differential_scores.R"

rule call_hits:
    input:
        input_scores="../outputs/gi_scores/{screen}/discriminant_scores/discriminant_scores_{score}.tsv"
    output:
        output_hits="../outputs/gi_scores/{screen}/discriminant_scores/discriminant_hits_{score}.tsv"
    log:
        "../outputs/logs/{screen}/{screen}_{score}_call_hits.log"
    params:
        score=lambda wildcards: wildcards.score,
        screen=lambda wildcards: wildcards.screen,
        threshold=config["HIT_THRESHOLD"]
    script:
        "../scripts/call_hits.R"

rule call_differential_hits:
    input:
        input_scores="../outputs/gi_scores/{screen}/differential_scores/discriminant_differential_scores_{rep}.tsv"
    output:
        output_hits="../outputs/gi_scores/{screen}/differential_scores/differential_hits_{rep}.tsv"
    log:
        "../outputs/logs/{screen}/{screen}_{rep}_call_differential_hits.log"
    params:
        rep=lambda wildcards: wildcards.rep,
        screen=lambda wildcards: wildcards.screen,
        threshold=config["DIFFERENTIAL_HIT_THRESHOLD"]
    script:
        "../scripts/call_differential_hits.R"

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

rule identify_all_hits:
    input:
        lambda wildcards: _expand_hit_targets(wildcards)

rule identify_differential_hits:
    input:
        lambda wildcards: _expand_differential_hit_targets(wildcards)
