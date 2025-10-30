rule compute_gamma_genetic_interaction_scores:
    input:
        input_orientation_indep_phenotypes="../outputs/phenotypes/{screen}_filtered_gamma_phenotypes.tsv",
        input_single_sgRNA_phenotypes="../outputs/phenotypes/{screen}_filtered_gamma_single_sgRNA_phenotypes.tsv"
    output:
        # constrain score wildcard to values beginning with 'Gamma'
        output_dir=directory("../outputs/gi_scores/{screen}/individual_scores/{score,Gamma.*}"),
        output_all_scores="../outputs/gi_scores/{screen}/construct_scores/all_gis_{score,Gamma.*}.tsv",
        output_model_estimates="../outputs/gi_scores/{screen}/models/model_estimates_{score,Gamma.*}.tsv",
        output_model_stats="../outputs/gi_scores/{screen}/models/model_stats_{score,Gamma.*}.tsv",
        output_workspace=temp("../outputs/gi_scores/{screen}/construct_scores/gi_workspace_{score,Gamma.*}.rds")
    log:
        "../outputs/logs/{screen}/{screen}_{score,Gamma.*}_compute_genetic_interaction_scores.log"
    params:
        screen=lambda wildcards: wildcards.screen,
        score=lambda wildcards: wildcards.score
    script:
        "../scripts/calculate_gi_scores.R"

rule compute_tau_genetic_interaction_scores:
    input:
        input_orientation_indep_phenotypes="../outputs/phenotypes/{screen}_filtered_tau_phenotypes.tsv",
        input_single_sgRNA_phenotypes="../outputs/phenotypes/{screen}_filtered_tau_single_sgRNA_phenotypes.tsv"
    output:
        # constrain score wildcard to values beginning with 'Tau' 
        output_dir=directory("../outputs/gi_scores/{screen}/individual_scores/{score,Tau.*}"),
        output_all_scores="../outputs/gi_scores/{screen}/construct_scores/all_gis_{score,Tau.*}.tsv",
        output_model_estimates="../outputs/gi_scores/{screen}/models/model_estimates_{score,Tau.*}.tsv",
        output_model_stats="../outputs/gi_scores/{screen}/models/model_stats_{score,Tau.*}.tsv",
        output_workspace=temp("../outputs/gi_scores/{screen}/construct_scores/gi_workspace_{score,Tau.*}.rds")
    log:
        "../outputs/logs/{screen}/{screen}_{score,Tau.*}_compute_genetic_interaction_scores.log"
    params:
        score=lambda wildcards: wildcards.score,
        screen=lambda wildcards: wildcards.screen
    script:
        "../scripts/calculate_gi_scores.R"

rule calculate_gene_level_scores:
    input:
        input_gi_scores="../outputs/gi_scores/{screen}/construct_scores/all_gis_{score}.tsv",
        input_gi_workspace="../outputs/gi_scores/{screen}/construct_scores/gi_workspace_{score}.rds",
        input_idmap="data/annotations/{screen}_id_to_name_mapping.tsv"
    output:
        output_gene_level_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv",
        output_gene_level_workspace=temp("../outputs/gi_scores/{screen}/gene_combination_scores/gene_level_workspace_{score}.rds")
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
        input_gene_level_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv",
        input_gene_level_workspace="../outputs/gi_scores/{screen}/gene_combination_scores/gene_level_workspace_{score}.rds"
    output:
        output_discriminant_scores="../outputs/gi_scores/{screen}/discriminant_scores/discriminant_scores_{score}.tsv",
        output_discriminant_workspace=temp("../outputs/gi_scores/{screen}/discriminant_scores/discriminant_workspace_{score}.rds")
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
        input_gamma_workspace="../outputs/gi_scores/{screen}/construct_scores/gi_workspace_Gamma.{rep}.rds",
        input_tau_workspace="../outputs/gi_scores/{screen}/construct_scores/gi_workspace_Tau.{rep}.rds",
        input_idmap="data/annotations/{screen}_id_to_name_mapping.tsv"
    output:
        output_differential_scores="../outputs/gi_scores/{screen}/differential_scores/differential_scores_{rep}.tsv",
        output_gene_differential_scores="../outputs/gi_scores/{screen}/differential_scores/gene_differential_scores_{rep}.tsv",
        output_discriminant_differential_scores="../outputs/gi_scores/{screen}/differential_scores/discriminant_differential_scores_{rep}.tsv",
        output_diff_workspace=temp("../outputs/gi_scores/{screen}/differential_scores/diff_scores_workspace_{rep}.rds")
    log:
        "../outputs/logs/{screen}/{screen}_{rep}_calculate_differential_scores.log"
    params:
        rep=lambda wildcards: wildcards.rep,
        screen=lambda wildcards: wildcards.screen
    script:
        "../scripts/calculate_differential_scores.R"

rule call_hits:
    input:
        input_scores="../outputs/gi_scores/{screen}/discriminant_scores/discriminant_scores_{score}.tsv",
        input_discriminant_workspace="../outputs/gi_scores/{screen}/discriminant_scores/discriminant_workspace_{score}.rds"
    output:
        output_hits="../outputs/gi_scores/{screen}/discriminant_scores/discriminant_hits_{score}.tsv",
        output_hits_workspace=temp("../outputs/gi_scores/{screen}/discriminant_scores/hits_workspace_{score}.rds")
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
        input_scores="../outputs/gi_scores/{screen}/differential_scores/discriminant_differential_scores_{rep}.tsv",
        input_diff_workspace="../outputs/gi_scores/{screen}/differential_scores/diff_scores_workspace_{rep}.rds"
    output:
        output_hits="../outputs/gi_scores/{screen}/differential_scores/differential_hits_{rep}.tsv",
        output_diff_hits_workspace=temp("../outputs/gi_scores/{screen}/differential_scores/diff_hits_workspace_{rep}.rds")
    log:
        "../outputs/logs/{screen}/{screen}_{rep}_call_differential_hits.log"
    params:
        rep=lambda wildcards: wildcards.rep,
        screen=lambda wildcards: wildcards.screen,
        threshold=config["DIFFERENTIAL_HIT_THRESHOLD"]
    script:
        "../scripts/call_differential_hits.R"

##############################################
# Wrapper rules to build all scores and hits #
##############################################

rule compute_all_genetic_interaction_scores:
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
